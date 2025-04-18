# Quickstart

### 1. Install
```bash
git clone https://github.com/JohnsonDylan/MassQLab.git
cd MassQLab

py -3.9 -m venv env
env\Scripts\activate

pip install -r requirements.txt
```

### 2. massqlab_config.json
```json
{
  "data_directory": "data/",
  "queryfile": "MassQL_Queries.csv"
}
```

### 3. Run
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

### 4. Results
Results are saved in "MassQLab_Output" directory within data_directory

____
