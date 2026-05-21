# MassQLab Install and Run Instructions

These instructions install MassQLab from GitHub into a fresh Python environment and run the command-line workflow.

The command-line format is:

    massqlab <query_file_or_query_directory> <mzML_file_directory>

Example:

    massqlab MassQL_Queries.csv data/

1. Clone the repository
```
git clone https://github.com/JohnsonDylan/MassQLab.git
cd MassQLab
git checkout dev
```
2. Create and activate a Python environment

Option A: conda

    conda create -n massqlab python=3.9
    conda activate massqlab

Option B: Python venv on macOS/Linux

    python -m venv .venv
    source .venv/bin/activate

3. Update installation tools
```
python -m pip install --upgrade pip setuptools wheel
```
4. Install MassQLab

From the root of the cloned repository:

    pip install -e .

5. Confirm the command is available
```
massqlab --help
```
You should see help text for the massqlab command.

6. Prepare input files

You need:

1. A MassQL query file, or a directory containing exactly one query file.
2. A directory containing one or more .mzML files.

Supported query file types:

    .csv
    .xlsx
    .xls
    .json


7. Run MassQLab

Basic usage:

    massqlab <query_file_or_query_directory> <mzML_file_directory>

Example with a direct query file:

    massqlab MassQL_Queries.csv data/

Example with a custom output directory:

    massqlab MassQL_Queries.csv data/ --output-directory results/


8. Expected outputs

MassQLab creates a timestamped output folder.

By default, the output folder is created inside the .mzML data directory:

    data/
      MassQLab_Output/
        massqlab_<timestamp>/

Expected output files include:

    massqlab_export_bundle.xlsx
    ms1_analysis_export_bundle.xlsx
    peak_gaussian_overlays.pdf
    area_heatmap.pdf
    peak_area_heatmap.pdf


Minimal full example

    git clone https://github.com/JohnsonDylan/MassQLab.git
    cd MassQLab
    git checkout dev

    conda create -n massqlab python=3.9
    conda activate massqlab

    python -m pip install --upgrade pip setuptools wheel
    pip install -e .

    massqlab MassQL_Queries.csv data/
