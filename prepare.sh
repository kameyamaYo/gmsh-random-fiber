sudo apt update
sudo apt install -y gmsh python3-venv python3-pip

# --- Python仮想環境を作り直す ---
python3 -m venv .venv
source .venv/bin/activate

# --- Pythonライブラリ ---
python -m pip install --upgrade pip
python -m pip install gmsh numpy meshio

# --- 確認 ---
gmsh --version
python -c "import gmsh; print('gmsh python =', gmsh.__file__)"