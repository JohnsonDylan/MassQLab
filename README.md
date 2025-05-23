# MassQLab
An Implementation of MassQL for Streamlining Common Mass Spectrometry Analysis Workflows

MassQLab applies a series of queries (written in the language of MassQL) to a directory containing mass spectrometry data in mzML format. Results are tabulated and saved as images and csv files.

MassQL: https://github.com/mwang87/MassQueryLanguage


# Installation

See [Setup Guide](setup_guide.md) or [Quickstart](quickstart.md) for Initial Setup Instructions 

## Usage - Entry Points

**MassQLab** currently supports multiple usage modes, depending on your workflow and preference:

### A) Command Line

- `py src\MassQLab_console.py` — Executes the full analysis pipeline from the command line [^]
- `py src\MassQLab_GUI.py` — Launches the graphical interface [^^]

### B) Jupyter Notebooks

- `notebooks\MassQLab_notebook.ipynb` — Walks through the workflow with manual control
- `notebooks\MassQLab_console.ipynb` — Executes the full analysis pipeline within a notebook [^]
- `notebooks\MassQLab_GUI.ipynb` — Launches the GUI from within a notebook [^^]

### C) Standalone GUI *(in development)*

- A full desktop interface is planned to support non-technical users and exploratory workflows, with minimal code and no need for Python or environment setup.

[^] Requires massqlab_config.json to be configured upfront

[^^] If GUI hangs or if closed while in progress, may continue running in background and may require force close via task manager

# Configuration

Each entry point loads configurations populated from the `massqlab_config.json` file located in the project root. 

Some entry points (e.g., MassQLab_console.py, MassQLab_console.ipynb) run the full workflow based only on this configuration. 

Other entry points (e.g., MassQLab_notebook.ipynb, MassQLab_GUI.py, MassQLab_GUI.ipynb) permit manual modification of configuration at runtime.

### Example `massqlab_config.json`

```json
{
  "data_directory": "path/to/data",
  "queryfile": "path/to/queryfile",
  "output_directory": false,
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

### Additional Config Parameters

For optional fields and advanced configuration, see the [Config File Advanced](#config-file-advanced) section.


# Query File

The **query file** (`queryfile`) defines the MassQL queries to be executed, along with a user-defined name for each query.

### Example Files

Example query files are provided in the root of the MassQLab repository:

- `MassQL_Queries.json`
- `MassQL_Queries.csv` (xlsx also permitted with same schema)

### Required Fields

Each query in the file must include the following:

- **`name`**:  
  A unique name for each query (per MS1/MS2 type). Avoid very long names to prevent output issues.

- **`query`**:  
  The [MassQL query](#massql-queries) string itself.

### Additional Query Parameters

For optional fields and advanced query configuration, see the [Query File Advanced](#query-file-advanced) section.


# MassQL Queries

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

# Launch

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

Results are saved in `MassQLab_Output` inside the defined `output_directory` (`data_directory` is used by default), organized by timestamp.

**MS1 Outputs**
- `ms1_raw_df.csv`: Raw query results + metadata  
- `ms1_analysis_df.csv`: Peak fitting analysis  
- `ms1_RT_analysis_df.csv`: RT stats vs. query RT
- `ms1_traces.pdf`: All Gaussian fits per query/file
- `ms1_summary_traces.pdf`: Summary of Gaussian fits per query
- `ms1_summary_traces_inverse.pdf`: Summary of Gaussian fits per file
- `ms1_summary_areas.pdf`: Peak areas per query  
- `ms1_summary_areas_inverse.pdf`: Peak areas per file
- `ms1_consolidated_traces.png`: Consolidated image of all ms1_traces with Gaussian fit

**MS2 Outputs**
- `ms2_raw_df.csv`: Raw query results + metadata  
- `ms2_analysis_df.csv`: Peak picking analysis  
- `ms2_plots.pdf`: Scan intensities (lines connect shared MS1 scans)  
- `ms2_summary_plots.pdf`: Top scan per query
- `ms2_cluster_plots.pdf`: Summary of MS2 analysis clustered by energy level if applicable
- `ms2_cluster_plots_group.pdf`: Summary of MS2 analysis by query group (if defined)



## Config File Advanced

- **`output_directory`**  
  Saves output to custom directory. Othewrise saved within data_directory.

- **`analysis`**  
  Whether to run downstream analysis on results returned from the MassQL queries (`true` or `false`).

- **`cache_setting`**  
  If `true`, saves an intermediate [Feather](https://arrow.apache.org/docs/python/feather.html) file to disk and uses it for faster re-querying.

- **`convert_raw`**  
  If `true`, attempts to convert raw files to `.mzML`.  
  *Note:* This feature is experimental. It is recommended to convert files using [ProteoWizard msconvert](https://proteowizard.sourceforge.io/) before using MassQLab.

- **`msconvertexe`**  
  Path to the `msconvert.exe` executable used for raw-to-mzML conversion (only used if `convert_raw` is `true`). [MSConvert](#msconvert)
  
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

## Peak Detection and Integration

Peaks are first detected from `y-data` vs `x-data` using `scipy.signal.find_peaks()` with the following thresholds:

- **Height**: greater than `mean(ydata) * 1.1`
- **Distance**: at least 5 observations apart
- **Prominence**: at least `max(ydata) / 10`

### Peak Filtering

Detected peaks are then filtered based on their x-axis value within a target retention time window:

```
(rtmin - range_rt) ≤ peak_x ≤ (rtmax + range_rt)
```

Where:

```
range_rt = rtmax - rtmin
```

Peaks are further refined by excluding those with excessive width:

```
FWHM < range_rt
```

### Peak Ranking and Analysis

The remaining peaks are ranked by their height, and the most prominent one is selected for detailed analysis.  
The Full Width at Half Maximum (FWHM) is calculated and used to approximate the peak as a Gaussian function.

```
σ = fwhm / (2 * sqrt(2 * log(2)))
peak_area = height × σ × sqrt(2π)
```

- **`peak_area_alt`**: Integrates the total signal using trapezoidal integration, without relying on peak picking.
- **`peak_intensity`**: Records the maximum y-value (intensity) from the data, also independent of peak picking.

## MSConvert
mzML files can be created from raw files on the fly if MSConvert (part of ProteoWizard software) is installed separately

https://proteowizard.sourceforge.io/download.html

You will need to identify the path to "msconvert.exe" and select it in the "msconvertexe" field in the MassQLab GUI window. You will also need to check the "convert_raw" button.

On my system, the path to mscovert.exe looks like this:
`"C:\Users\username\AppData\Local\Apps\ProteoWizard 3.0.23118.b2ed96f 64-bit\msconvert.exe"`

____
