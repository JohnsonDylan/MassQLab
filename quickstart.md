# Quickstart

### 1. Clone
```bash
git clone https://github.com/JohnsonDylan/MassQLab.git
cd MassQLab
```

### 2. virtualenv OR conda
if using virtualenv:
```bash
py -3.9 -m venv env
env\Scripts\activate

pip install jupyterlab massql reportlab openpyxl matplotlib
pip install -r requirements.txt
```

if using conda:
```bash
conda create -n env python=3.9
conda activate env
```

### 3. Dependendencies
```bash
pip install jupyterlab massql reportlab openpyxl matplotlib
```
OR
```bash
pip install -r requirements.txt
```

### 4. massqlab_config.json
```json
{
  "data_directory": "data/",
  "queryfile": "MassQL_Queries.csv"
}
```

### 5. Run
```bash
py src\MassQLab_console.py      # Run full pipeline via console
```
OR
```bash
py src\MassQLab_GUI.py          # Launch GUI (if supported)
```
OR
```bash
jupyter lab                     # Launch Jupyter Lab, navigate to notebooks directory and launch MassQLab_notebook.ipynb
```

### 6. Results
Results are saved in "MassQLab_Output" directory within data_directory

____
