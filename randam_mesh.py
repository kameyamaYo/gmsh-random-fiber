# -*- coding: utf-8 -*-
import gmsh
import math
import sys


# ============================================================
# 入力読み込み
# ============================================================

def read_centers(fname_center):
    """
    coord_center.dat を読む
    想定形式:
    1   x1   y1
    2   x2   y2
    ...
    """
    centers = []
    with open(fname_center, "r") as f:
        header = f.readline()  # "model number ..." を飛ばす
        for line in f:
            parts = line.split()
            if len(parts) < 3:
                continue
            # 1列目はID, 2列目x, 3列目y
            idx = int(parts[0])
            x = float(parts[1])
            y = float(parts[2])
            centers.append((idx, x, y))
    return centers


def read_info(fname_info):
    """
    input_info.dat を読む
    """
    keys = [
        "time(sec)",
        "Vf",
        "a",
        "b",
        "R",
        "C(r_ratio)",
        "max_x",
        "max_y",
        "volume"
    ]
    params = {}
    with open(fname_info, "r") as f:
        for key in keys:
            line = f.readline()
            if not line:
                raise ValueError("input_info.dat の行数が不足しています")
            parts = line.strip().split()
            params[key] = float(parts[-1])
    return params


# ============================================================
# 幾何処理
# ============================================================

def classify_y_crossing(centers, Ly, r):
    """
    各円がY境界にかかるか判定する
    return:
        lower_cross : 下端にかかる円IDのリスト
        upper_cross : 上端にかかる円IDのリスト
    """
    lower_cross = []
    upper_cross = []

    for idx, x, y in centers:
        if y - r < 0.0:
            lower_cross.append(idx)
        if y + r > Ly:
            upper_cross.append(idx)

    return lower_cross, upper_cross


def build_geometry_yperiodic(centers, Lx, Ly, r):
    """
    Y周期コピー付きで幾何を作る
    """
    occ = gmsh.model.occ

    # 母材矩形
    box = occ.addRectangle(0.0, 0.0, 0.0, Lx, Ly)

    circles = []          # (2,tag) のリスト
    circle_origin = {}    # tag -> (original fiber id, shift)

    # shift = 0 : 元の円
    # shift = +1: y+Ly にコピー
    # shift = -1: y-Ly にコピー
    for idx, x, y in centers:
        # 元の円
        c0 = occ.addDisk(x, y, 0.0, r, r)
        circles.append((2, c0))
        circle_origin[c0] = (idx, 0)

        # 下端にかかるなら、上側へ周期コピー
        if y - r < 0.0:
            c_up = occ.addDisk(x, y + Ly, 0.0, r, r)
            circles.append((2, c_up))
            circle_origin[c_up] = (idx, +1)

        # 上端にかかるなら、下側へ周期コピー
        if y + r > Ly:
            c_dn = occ.addDisk(x, y - Ly, 0.0, r, r)
            circles.append((2, c_dn))
            circle_origin[c_dn] = (idx, -1)

    occ.synchronize()

    # 長方形との交差をとって、セル内の円部分だけ残す
    intersected_circles, _ = occ.intersect(
        [(2, box)], circles,
        removeObject=False,
        removeTool=True
    )
    occ.synchronize()

    # 母材から円を引く
    cut_result, _ = occ.cut(
        [(2, box)], intersected_circles,
        removeObject=True,
        removeTool=False
    )
    occ.synchronize()

    return cut_result, intersected_circles


# ============================================================
# 境界抽出
# ============================================================

def find_boundary_lines(Lx, Ly, matrix_surfaces=None, tol=1e-6):
    """
    外周境界線を抽出する。
    可能なら matrix surface の境界だけを見て、そこから
    left/right/top/bottom を分類する。
    """
    left_lines = []
    right_lines = []
    bottom_lines = []
    top_lines = []

    # どの曲線集合を対象にするか
    if matrix_surfaces is not None and len(matrix_surfaces) > 0:
        # matrix surface の境界だけを見る
        boundary = gmsh.model.getBoundary(
            [(2, s) for s in matrix_surfaces],
            combined=False,
            oriented=False
        )
        curve_entities = sorted(set([b for b in boundary if b[0] == 1]))
    else:
        # fallback: 全曲線を見る
        curve_entities = gmsh.model.getEntities(1)

    for dim, tag in curve_entities:
        xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(dim, tag)

        # 左右
        if abs(xmin - 0.0) < tol and abs(xmax - 0.0) < tol:
            left_lines.append(tag)
        elif abs(xmin - Lx) < tol and abs(xmax - Lx) < tol:
            right_lines.append(tag)

        # 上下
        elif abs(ymin - 0.0) < tol and abs(ymax - 0.0) < tol:
            bottom_lines.append(tag)
        elif abs(ymin - Ly) < tol and abs(ymax - Ly) < tol:
            top_lines.append(tag)

    return (
        sorted(set(left_lines)),
        sorted(set(right_lines)),
        sorted(set(bottom_lines)),
        sorted(set(top_lines)),
    )


def curve_midpoint(tag):
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(1, tag)
    return 0.5 * (xmin + xmax), 0.5 * (ymin + ymax)


def sort_boundary_lines_for_periodic(bottom_lines, top_lines):
    """
    bottom/top を x 座標順に並べる
    """
    bottom_lines = sorted(bottom_lines, key=lambda t: curve_midpoint(t)[0])
    top_lines    = sorted(top_lines,    key=lambda t: curve_midpoint(t)[0])
    return bottom_lines, top_lines


# ============================================================
# メッシュ制御
# ============================================================

def get_circle_boundary_curves(circle_surfaces):
    """
    円表面群の境界曲線を集める
    """
    curves = []
    for dim, s in circle_surfaces:
        bd = gmsh.model.getBoundary([(2, s)], oriented=False)
        curves.extend([c[1] for c in bd if c[0] == 1])
    return sorted(list(set(curves)))


def setup_mesh_field(circle_surfaces, Lx, Ly, r):
    """
    円境界と外周を細かくする
    """
    left_lines, right_lines, bottom_lines, top_lines = find_boundary_lines(Lx, Ly)

    refine_curves = []
    refine_curves.extend(get_circle_boundary_curves(circle_surfaces))
    refine_curves.extend(left_lines + right_lines + bottom_lines + top_lines)
    refine_curves = sorted(list(set(refine_curves)))

    # 要素サイズ設定
    lc_fiber = r / 6.0
    lc_bulk  = r / 2.0

    field_dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(field_dist, "CurvesList", refine_curves)

    field_thr = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(field_thr, "IField", field_dist)
    gmsh.model.mesh.field.setNumber(field_thr, "LcMin", lc_fiber)
    gmsh.model.mesh.field.setNumber(field_thr, "LcMax", lc_bulk)
    gmsh.model.mesh.field.setNumber(field_thr, "DistMin", 0.0)
    gmsh.model.mesh.field.setNumber(field_thr, "DistMax", 2.5 * r)

    gmsh.model.mesh.field.setAsBackgroundMesh(field_thr)


# ============================================================
# Physical group
# ============================================================

def add_physical_groups(cut_result, intersected_circles, Lx, Ly):
    matrix_surfaces = [s[1] for s in cut_result]
    fiber_surfaces  = [s[1] for s in intersected_circles]

    # 2D physical groups
    pg = gmsh.model.addPhysicalGroup(2, matrix_surfaces, tag=1)
    gmsh.model.setPhysicalName(2, pg, "matrix")

    pg = gmsh.model.addPhysicalGroup(2, fiber_surfaces, tag=2)
    gmsh.model.setPhysicalName(2, pg, "fiber")

    # 1D boundary groups
    left_lines, right_lines, bottom_lines, top_lines = find_boundary_lines(Lx, Ly)

    if left_lines:
        pg = gmsh.model.addPhysicalGroup(1, left_lines, tag=11)
        gmsh.model.setPhysicalName(1, pg, "left")
    if right_lines:
        pg = gmsh.model.addPhysicalGroup(1, right_lines, tag=12)
        gmsh.model.setPhysicalName(1, pg, "right")
    if bottom_lines:
        pg = gmsh.model.addPhysicalGroup(1, bottom_lines, tag=13)
        gmsh.model.setPhysicalName(1, pg, "bottom")
    if top_lines:
        pg = gmsh.model.addPhysicalGroup(1, top_lines, tag=14)
        gmsh.model.setPhysicalName(1, pg, "top")


# ============================================================
# 周期境界
# ============================================================

def set_y_periodic(Ly, Lx):
    left_lines, right_lines, bottom_lines, top_lines = find_boundary_lines(Lx, Ly)

    if not bottom_lines or not top_lines:
        raise RuntimeError("bottom/top boundary lines が見つかりません")

    bottom_lines, top_lines = sort_boundary_lines_for_periodic(bottom_lines, top_lines)

    if len(bottom_lines) != len(top_lines):
        raise RuntimeError(
            f"周期境界線の本数が一致しません: "
            f"bottom={len(bottom_lines)}, top={len(top_lines)}"
        )

    # top -> bottom へ y方向 -Ly の平行移動
    affine = [
        1, 0, 0, 0,
        0, 1, 0, -Ly,
        0, 0, 1, 0,
        0, 0, 0, 1
    ]

    gmsh.model.mesh.setPeriodic(1, top_lines, bottom_lines, affine)


# ============================================================
# main
# ============================================================

def main():
    fname_center = sys.argv[1] if len(sys.argv) > 1 else "input/coord_center.dat"
    fname_info   = sys.argv[2] if len(sys.argv) > 2 else "input/info.dat"
    fname_out    = sys.argv[3] if len(sys.argv) > 3 else "output/mesh_yperiodic.msh"

    centers = read_centers(fname_center)
    params  = read_info(fname_info)

    Lx = params["max_x"]
    Ly = params["max_y"]
    r  = params["R"]

    print(f"Lx = {Lx}")
    print(f"Ly = {Ly}")
    print(f"r  = {r}")
    print(f"number of fibers = {len(centers)}")

    lower_cross, upper_cross = classify_y_crossing(centers, Ly, r)
    print("lower crossing fibers:", lower_cross)
    print("upper crossing fibers:", upper_cross)

    gmsh.initialize()
    gmsh.model.add("RandomFiber_YPeriodic")

    try:
        cut_result, intersected_circles = build_geometry_yperiodic(centers, Lx, Ly, r)

        add_physical_groups(cut_result, intersected_circles, Lx, Ly)
        setup_mesh_field(intersected_circles, Lx, Ly, r)

        # まず periodic を設定
        set_y_periodic(Ly, Lx)

        # メッシュオプション
        gmsh.option.setNumber("Mesh.Algorithm", 6)   # Frontal-Delaunay
        gmsh.option.setNumber("Mesh.RecombineAll", 0)  # まずは三角形で安定に
        gmsh.option.setNumber("Mesh.Optimize", 1)
        gmsh.option.setNumber("Mesh.OptimizeNetgen", 1)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)

        gmsh.model.mesh.generate(2)
        gmsh.write(fname_out)

        print(f"[OK] mesh written to {fname_out}")

        # GUI確認したいならコメントアウトを外す
        # gmsh.fltk.run()

    finally:
        gmsh.finalize()


if __name__ == "__main__":
    main()