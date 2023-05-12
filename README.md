# SpectraSpectre
An Implementation of MassQL for Streamlining Common Mass Spectrometry Analysis Workflows

SpectraSpectra applies a series of queries (written in the language of MassQL) to a directory containing mass spectrometry data in mzML format. Results are tabulated and saved as images and xlsx documents.

MassQL: https://github.com/mwang87/MassQueryLanguage

## Setup
Project was developed using Python 3.9 (reccomended). A "requirements.txt" file is included for creating a virtualenv (virtual environment). See file "venv_setup.txt" for basic instructions for setting up a virtualenv. 

Note: MassQL package currently not available via conda and therefor anaconda environments are not reccomended.

## Usage
Default usage:
  1. Put mzML files into "SpectraSpectre/data/", replacing "placeholder.txt"
  2. Open "MassQL_Queries.xlsx" and modify queries (see Query File section below)
  3. Run "SpectraSpectre.ipynb" in Jupyter OR run "SpectraSpectre.py" from command prompt
  4. Results will be output in "SpectraSpectre/data/SpectraSpectre_Output/"

Custom usage:
  The file "spectre_config.json" can be modified to change the default configuration.
  data_directory: change this default path to the directory where your mzML format data is stored
	queryfile: change this filename (or include a full path) to an alternate MassQL query file

## Query File
The query file (MassQL_Queries.xlsx by default) is used to define the parameters that will compose a MassQL query.
The column headers should not be changed as they are required by SpectraSpectre.
Example queries are included by default but should be modified to accomodate your data.

Description of column headers:
Name: Name of query, likely the component to be identified
KEGG: KEGG identifier (currently unused)
Formula: Molecular formula (currently unused)
Monoisotopic: The monoisotopic mass of component to be identified. 
ion_mode: Use 1 for positive mode, 0 for negative mode
TOLERANCEPPM: The tolerance for the m/z in parts per million
RTMIN: The minimum value in the retention time range
RTMAX: The maximum value in the retention time range
QC_threshold: The tolerance for validating files versus QC standard(s). ie. 0.2 refers to a 20% threshold. (Not required for regular usaged, see "Sytem Suitability" section below)

## System Suitability
SpectraSpectre is useful for rapidly validating analytical instrumentation performance by comparing new data versus an established quality control.
SpectraSpectre will automatically compare non-QC samples versus QC samples to determine if the queried data meets a required threshold for system suitability.

To use the system suitability function, define your QC samples in one of the following two ways:
1) Prefix your QC mzML filename with "QC_"
or
2) Open "spectre_config.json" and add the name(s) of the QC file(s) to the field "QC_files"

If multiple QC files are present, the QC files will be averaged

An additional xlsx file will be generated and saved in "SpectraSpectre/data/SpectraSpectre_Output/"

The QC output file contains 4 sheets. The first sheet gives an overall determination of system suitability for each non-QC file (TRUE means the file passed the system suitability check). A passing system suitability check means each query applied to the file produced a resultant peak area that was withing the QC_threshold tolerance of the average of the QC_files. 

