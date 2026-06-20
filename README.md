# Random fiber Y-periodic mesh generator

ランダム繊維配列に対して、Y方向周期境界を持つ2Dメッシュを生成するスクリプトです。  
Gmsh Python API を用いて、矩形セル中に円形繊維を配置し、母材領域と繊維領域を分割します。

## Files

- `randam_mesh.py` : メッシュ生成コード
- `input/coord_center.dat` : 繊維中心座標
- `input/info.dat` : 半径、セルサイズなどの情報
- `output/` : 生成メッシュ出力先

## Input format

### coord_center.dat
```txt
model number 1.000
1 0.2766842014 0.3943136912
2 1.0236093498 0.3307884051
...
