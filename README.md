# MassQLab
An Implementation of MassQL for Streamlining Common Mass Spectrometry Analysis Workflows

MassQLab applies a series of queries (written in the language of MassQL) to a directory containing mass spectrometry data in mzML format. Results are tabulated and saved as images and csv files.

MassQL: https://github.com/mwang87/MassQueryLanguage


## Usage

See [Setup Guide](setup_guide.md) for Initial Setup Instructions

### Entry Points

**MassQLab** currently supports multiple usage modes, depending on your workflow and preference:

#### A) Command Line

- `src\MassQLab_console.py` — Executes the full analysis pipeline from the command line [^]
- `src\MassQLab_GUI.py` — Launches the graphical interface [^^]

#### B) Jupyter Notebooks

- `notebooks/MassQLab_notebook.ipynb` — Walks through the workflow with manual control
- `notebooks/MassQLab_console.ipynb` — Executes the full analysis pipeline within a notebook [^]
- `notebooks/MassQLab_GUI.ipynb` — Launches the GUI from within a notebook [^^]

#### C) Standalone GUI *(in development)*

- A full desktop interface is planned to support non-technical users and exploratory workflows, with minimal code and no need for Python or environment setup.

[^] Pulls configuration settings from massqlab_config.json

[^^] If GUI hangs or if closed while in progress, may continue running in background and may require force close via task manager

## Configuration

MassQLab configuration is determined by the point of entry (e.g., notebook or command line), with settings automatically populated from the `massqlab_config.json` file located in the project root.

### Example `massqlab_config.json`

```json
{
  "data_directory": "path/to/data",
  "queryfile": "path/to/queryfile",
  "analysis": true,
  "cache_setting": true,
  "convert_raw": false,
  "msconvertexe": false
}
```

### Configuration Fields

- **`data_directory`**  
  Path to the folder containing mass spectrometry files (`.mzML` or raw format).  
  MassQLab currently processes **all files** in this directory.

- **`queryfile`**  
  Path to the query file (JSON, CSV, or XLSX) that defines the MassQL queries to apply to each file in the data directory.

- **`analysis`**  
  Whether to run downstream analysis on results returned from the MassQL queries (`true` or `false`).

- **`cache_setting`**  
  If `true`, saves an intermediate [Feather](https://arrow.apache.org/docs/python/feather.html) file to disk and uses it for faster re-querying.

- **`convert_raw`**  
  If `true`, attempts to convert raw files to `.mzML`.  
  *Note:* This feature is experimental. It is recommended to convert files using [ProteoWizard msconvert](https://proteowizard.sourceforge.io/) before using MassQLab.

- **`msconvertexe`**  
  Path to the `msconvert.exe` executable used for raw-to-mzML conversion (only used if `convert_raw` is `true`). [MSConvert](#msconvert)

### MassQLab_notebook and GUI

These options still ingest massqlab_config.json but allow manual override and control of configuration


## Query File

The **query file** (`queryfile`) defines the MassQL queries to be executed, along with a user-defined name for each query. MassQLab supports query files in the following formats:

- `.json`
- `.csv`
- `.xlsx`

### Example Files

Example query files are provided in the root of the MassQLab repository:

- `MassQL_Queries.json`
- `MassQL_Queries.csv`

### Required Fields

Each query in the file must include the following:

- **`name`**:  
  A unique identifier for the query used in result files and visualizations.  
  - Should be unique within the same query type (MS1 or MS2).  
  - Duplicate names are allowed across types (e.g., same name used for an MS1 and an MS2 query).  
  - Note: Very long query names may cause issues in some output files.

- **`query`**:  
  The MassQL query string itself.  
  - Queries fall into two categories: **MS1** and **MS2**, which are processed independently.  
  - See the **"MassQL Queries"** section for detailed query syntax.

### Additional Parameters

For optional fields and advanced query configuration, see the [Query File Advanced](#query-file-advanced) section.


## MassQL Queries

See https://mwang87.github.io/MassQueryLanguage_Documentation/ for full documentation

### MS1 Query Examples

- Get MS1 scans matching m/z 207.1418 with a tolerance of 2.5 ppm, retention time between 1.0 and 1.2 minutes, and filter for 207.1418 peak intensity  
  `QUERY scaninfo(MS1DATA) FILTER MS1MZ=207.1418:TOLERANCEPPM=2.5 AND RTMIN=1.0 AND RTMAX=1.2`

- Get MS1 scans where an MS2 scan with a product ion is present  
  `QUERY scaninfo(MS1DATA) WHERE MS2PROD=226.18`

---

### MS2 Query Examples

- Get MS2 scans with a precursor ion matching m/z 429.3765, retention time between 9.0 and 9.5 minutes, and return total intensity of each scan  
  `QUERY scaninfo(MS2DATA) WHERE MS2PREC=429.3765:TOLERANCEPPM=2.5 AND RTMIN=9.0 AND RTMAX=9.5`

- Get MS2 scans with a precursor ion matching m/z 429.3765, retention time between 9.0 and 9.5 minutes, and return intensity of the peak with m/z 85.0281  
  `QUERY scaninfo(MS2DATA) WHERE MS2PREC=429.3765:TOLERANCEPPM=2.5 AND RTMIN=9.0 AND RTMAX=9.5 FILTER MS2PROD=85.0281:TOLERANCEPPM=10`

- Get MS2 scans where both a product ion and a neutral loss are present  
  `QUERY scaninfo(MS2DATA) WHERE MS2NL=176.0321 AND MS2PROD=85.02915`

- Get MS2 scans where a product ion matches an arithmetic expression  
  `QUERY scaninfo(MS2DATA) WHERE MS2PROD=144+formula(CH2)`

## Launch

### Command Line

From the project root, ensure your virtual environment is activated and all requirements are installed [(see Setup Guide)](setup_guide.md). Then run one of the following:

```bash
py src\MassQLab_console.py      # Launches the console workflow
```
```bash
py src\MassQLab_GUI.py          # Launches the GUI (if available)
```

### Notebook

From the project root, ensure your virtual environment is activated and all requirements are installed. Then launch Jupyter Lab:

```bash
jupyter lab
```

Navigate to the `notebooks/` directory and open the desired notebook to begin working with MassQLab.


## Results

Results can be found in MassQLab_Output directory located within the defined data_directory. MassQLab_Output directory contains timestamped subdirectories.

- ms1_raw_df.csv: Raw output returned from applying MS1 MassQL queries and merged with query metadata
- ms1_analysis_df.csv: Output after peak fitting analysis of ms1_raw_df
- ms1_RT_analysis_df.csv: Statistical analysis of MS1 peak center retention time relative to retention defined in query
- ms1_traces.pdf: Document containing plots showing gaussian fit of each query applied to each file
- ms1_summary_traces.pdf: Document containing plots of all peaks for each query
- ms1_summary_areas.pdf: Document containing plots of areas of all peaks for each query
- ms1_summary_areas_inverse.pdf: Document containing plots of areas of all peaks for each file

- ms2_raw_df.csv: Raw output returned from applying MS2 MassQL queries and merged with query metadata
- ms2_analysis_df.csv: Output after peak picking analysis of ms2_raw_df
- ms2_plots.pdf: Document containing plots of intensity of each scan returned for each file and query
  - Note: lines are shown between points that share the same MS1 scan
  - Note: highest intensity scan for each condition is used for downstream analysis
- ms2_summary_plots.pdf: Document containing plots of highest intensity scan for each query
- ms2_cluster_plots_group.pdf: Document containing summary plots of each group of queries as defined by "group" parameter in queryfile if present


## Query File Advanced
This section describes parameters that may be used in the queryfile in addition to "name" and "query"

Any additional parameters in queryfile will also be carried over into excel/csv outputs for ad-hoc usage

- "abundance":
    - an expected or anticipated abundance of the peak area (ms1) or intensity (MS2)
    - abundance validity will be measured with a 10% (ms1) or 20% (ms2) threshold of this value
        - abundances not within this threshold will be flagged as invalid and will not be considered in downstream analysis
        - customizing threshold will be added in future release 

- "group":
    - A designation for queries that have some relationship
    - Use the same argument for the "group" parameter in multiple queries to associate the queries with each other
    - Use case 1: Assign acyl-carnitine query and each acyl-carnite subfragment query a "group" argument of "acyl-carnitine"
        - In "ms2_cluster_plots_group" result, each group will be consolidated and visualized together
    - Use case 2: Assign all phosphatidylcholine lipids a "group" argument of "PC"
        - Function will be expanded in future release to automate quantifying total group (ie. PC) content

- "KEGG":
    - Argument does not currently carry any special meaning, but this and any other parameter in the query file will be carried over into excel/csv outputs
    - The query file is a good place to define unique compound identifiers for ad-hoc usage
        - metabolites: KEGG, InChI string/key, CASRN, SMILES, IUPAC, PubChem ID
        - proteins: PDB ID, Gene name, UniProt, Ensembl, RefSeq

## MassQLab Query Details  
1. MS1 queries will return a dataframe of total intensity of MS1 scan vs retention time for each MS1 scan that matches query
    - Peak will be fit with a gaussian and area will be determined from gaussian.
    - Peak must meet detection criteria to be considered valid
      - peak prominence > Intensity_Max / 10
      - peak height > Intensity_Average x 1.1
      - peak center > (RTMIN - (RTMAX - RTMIN))
      - peak center < (RTMAX + (RTMAX - RTMIN))
      - fwhm < (RTMAX - RTMIN)


2. MS2 queries will return a dataframe of total intensity of MS2 scan vs retention time for each MS2 scan that matches query 
   - Downstream analysis will split MS2 query results for each file based on any collision energy grouping  as defined in scan metadata
     - ie. HCD and CID will be split
     - ie. Collision energy 20 vs 30 vs 40 will be split
     - An MS2 group will be each combination, ie. HCD_20, CID_30, etc. 
   - When more than one valid scan is returned from the MassQL query for an MS2 group, the highest intensity scan will be used for downstream analysis

## MSConvert
mzML files can be created from raw files on the fly if MSConvert (part of ProteoWizard software) is installed separately

https://proteowizard.sourceforge.io/download.html

You will need to identify the path to "msconvert.exe" and select it in the "msconvertexe" field in the MassQLab GUI window. You will also need to check the "convert_raw" button.

On my system, the path to mscovert.exe looks like this:
`"C:\Users\username\AppData\Local\Apps\ProteoWizard 3.0.23118.b2ed96f 64-bit\msconvert.exe"`

