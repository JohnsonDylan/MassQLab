# SpectraSpectre
An Implementation of MassQL for Streamlining Common Mass Spectrometry Analysis Workflows

SpectraSpectra applies a series of queries (written in the language of MassQL) to a directory containing mass spectrometry data in mzML format. Results are tabulated and saved as images and xlsx documents.

MassQL: https://github.com/mwang87/MassQueryLanguage

## Setup
Project was developed using Python 3.9 (recommended). A "requirements.txt" file is included for creating a virtualenv (virtual environment). See file "venv_setup.txt" for basic instructions for setting up a virtualenv. 

Note: MassQL package is currently not available via conda, therefore anaconda environments are not recommended.

## Usage
Default usage:
>  1. Put mzML files into "SpectraSpectre/data/", replacing "placeholder.txt"
>  2. Open "MassQL_Queries.xlsx" and modify queries (see Query File section below)
>  3. Run "SpectraSpectre.exe" OR "SpectraSpectre.ipynb" in Jupyter OR run "SpectraSpectre.py" from command prompt
>  4. Results will be output in "SpectraSpectre/data/SpectraSpectre_Output/"

Custom usage:
>  The file "spectre_config.json" can be modified to change the default configuration.  
>  	* data_directory: change this default path to the directory where your mzML format data is stored  
>  	* queryfile: change this filename (or include a full path) to an alternate MassQL query file  

## Query File
The query file is used to define the parameters that will compose a MassQL query.

The query file to be used is defined in spectre_config.json (MassQL_Queries.json is used by default).

Example queries are included by default but should be modified to acommodate your data.

Only "name" and "query" are required to compose a valid query. All other parameters are optional.

### Parameters:

"name":
>  	* required
>  	* a name that will be associated with a query and used for identification
>  	* use a unique name to prevent errors

"query":
>  	* required
>  	* The query in MassQL format
>  	* https://mwang87.github.io/MassQueryLanguage_Documentation/

"abundance":
>  	* optional
>  	* an anticipated abundance of the peak area (ms1) or intensity (MS2)
>  	* abundance validity will be measured with a 10% (ms1) or 20% (ms2) threshold of this value


## spectre_config
spectre_config.json is used to define configuration settings that will be loaded upon program initialization.

### Parameters:

"data_directory":
>  	* path to data directory
>  	* defaults to relative directory "data/" if not provided
>  	* specify relative directory if within same directory as working directory, otherwise, use full path name. Json requires forward slashes in path name.
 
"queryfile":
>  	* path to json format query file
>  	* defaults to relative filepath "MassQL_Queries.json" if not provided
>  	* specify relative filepath if within same directory as working directory, otherwise, use full path name
 
"convert_raw": false,
>  	* set to true if you have MSConvert installed and want to convert raw files to mzML
 
"msconvert_exe": false,
>  	* path to msconvert.exe
>  	* required if you want to convert raw files to mzML

"cache": true
>  	* use MassQL cache function (default)
>  	* if true, feather files will be saved in your data directory alongside mzML files
>  	* this makes subsequent analysis of the same data files much faster
