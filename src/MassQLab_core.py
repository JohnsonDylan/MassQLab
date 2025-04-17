# Core MassQLab functions

import sys
import io, os, json, time, fnmatch, glob, warnings, subprocess, pandas as pd, numpy as np, regex as re, contextlib, textwrap
from pathlib import Path
from pyteomics import mzxml, mzml, auxiliary
from scipy.integrate import trapz
from scipy.signal import find_peaks, peak_widths
from scipy.interpolate import interp1d
from itertools import product
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from reportlab.lib.pagesizes import letter, A4, landscape
from reportlab.lib.units import inch
from reportlab.lib import colors, utils
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.platypus import Frame, PageTemplate, BaseDocTemplate, Image, Table, TableStyle, Paragraph, NextPageTemplate, PageBreak
from massql import msql_fileloading, msql_engine
from PIL import Image as RLImage


warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=DeprecationWarning)
pd.options.mode.chained_assignment = None  # default='warn'


def find_file(filename):
    filename = os.path.normpath(filename)
    base_filename = os.path.basename(filename)
    try:
        script_path = os.path.abspath(__file__)
    except NameError:
        script_path = os.getcwd()
    script_dir = os.path.dirname(script_path)
    script_dir_one_level_up = os.path.dirname(script_dir)
    filename_path_one_level_up = os.path.join(script_dir_one_level_up, base_filename)
    filename_path_current_wd = os.path.join(os.getcwd(), base_filename)
    filename_path_script_dir = os.path.join(script_dir, base_filename)
    if os.path.exists(filename):
        return filename
    elif os.path.exists(filename_path_one_level_up):
        return filename_path_one_level_up
    elif os.path.exists(filename_path_current_wd):
        return filename_path_current_wd
    elif os.path.exists(filename_path_script_dir):
        return filename_path_script_dir
    else: 
        return None

def find_dir(filedir):
    try:
        script_path = os.path.abspath(__file__)
    except NameError:
        script_path = os.getcwd()
    script_dir = os.path.dirname(script_path)
    script_dir_one_level_up = os.path.dirname(script_dir)
    filedir_one_level_up = os.path.join(script_dir_one_level_up, filedir)
    filedir_current_wd = os.path.join(os.getcwd(), filedir)
    filedir_script_dir = os.path.join(script_dir, filedir)
    if os.path.exists(filedir):
        return filedir
    elif os.path.exists(filedir_one_level_up):
        return filedir_one_level_up
    elif os.path.exists(filedir_current_wd):
        return filedir_current_wd
    elif os.path.exists(filedir_script_dir):
        return filedir_script_dir
    else: 
        return None





"""Configure"""

def configure_MassQLab(config_file_path="massqlab_config.json"):
    sys.stdout.write(f"\nInitialize MassQLab\n")
    sys.stdout.flush()
    config = {}

    config_file_path_loc = find_file(config_file_path)
    if config_file_path_loc:
        try:
            with open(config_file_path_loc, 'r') as config_file:
                config = json.load(config_file)
                sys.stdout.write(f"\nLoaded config file: {os.path.normpath(config_file_path_loc)}\n")
                sys.stdout.flush()
        except json.JSONDecodeError:
            sys.stdout.write(f"\nError: {config_file_path_loc} contains invalid JSON.\n")
            sys.stdout.flush()
    else:
        # sys.stdout.write(f"\nconfig_file: {config_file_path} could not be located.\nLoading defaults...\n")
        # sys.stdout.flush()
        pass

    # Extract configuration values with defaults if not found
    data_directory = config.get('data_directory', 'data/')
    queryfile = config.get('queryfile', 'MassQL_Queries.json')
    metadata_file = config.get('metadata_file', False)
    metadata_filename_column = config.get('metadata_filename_column', "CORE_Filename")
    metadata_group_columns = config.get('metadata_group_columns', "USER_Column1")
    kegg_path = config.get('kegg_path', False)
    convert_raw = config.get('convert_raw', False)
    msconvertexe = config.get('msconvert_exe', False)
    cache_setting = config.get('cache', True)
    datasaver = config.get('datasaver', False)
    analysis = config.get('analysis', True)

    # Check if data_directory and queryfile paths exist
    data_directory = os.path.normpath(data_directory)
    data_directory_loc = find_dir(data_directory)
    if data_directory_loc:
        sys.stdout.write(f"\ndata_directory: {os.path.normpath(data_directory_loc)}\n")
        sys.stdout.flush()
    else:
        # sys.stdout.write(f"\nWarning: data_directory '{os.path.normpath(data_directory)}' could not be located.\n")
        # sys.stdout.flush()
        pass

    queryfile = os.path.normpath(queryfile)
    queryfile_loc = find_file(queryfile)
    if queryfile_loc:
        sys.stdout.write(f"\nqueryfile: {os.path.normpath(queryfile_loc)}\n")
        sys.stdout.flush()
    else:
        # sys.stdout.write(f"\nWarning: queryfile '{queryfile}' could not be located.\n")
        # sys.stdout.flush()
        pass
    sys.stdout.write(f"\n")
    sys.stdout.flush()
    return data_directory_loc, queryfile_loc, metadata_file, metadata_filename_column, metadata_group_columns, kegg_path, convert_raw, msconvertexe, cache_setting, datasaver, analysis
    

"""Create Queries from Query File"""

# Function to extract RTMIN value from a given input string
def extract_rtmin_value(input_string):
    pattern = r'RTMIN=(\d+(?:\.\d+)?)'
    match = re.search(pattern, input_string)
    if match:
        return float(match.group(1))
    else:
        return 0

# Function to extract RTMAX value from a given input string
def extract_rtmax_value(input_string):
    pattern = r'RTMAX=(\d+(?:\.\d+)?)'
    match = re.search(pattern, input_string)
    if match:
        return float(match.group(1))
    else:
        return 99999
        
def create_queries(queryfile, queries=None, query_groups=None, name_kegg_dict=None):
    if name_kegg_dict is None:
        name_kegg_dict = {}
    if query_groups is None:
        query_groups = {}
    if queries is None:
        queries = []

    ms1_query_df = pd.DataFrame()
    ms2_query_df = pd.DataFrame()

    if queryfile:
        try:
            # Determine the file extension
            _, file_extension = os.path.splitext(queryfile)

            if file_extension.lower() == '.json':
                with open(queryfile) as queryfilej:
                    queryjson = json.load(queryfilej)
            elif file_extension.lower() == '.csv':
                df = pd.read_csv(queryfile)
                queryjson = df.to_dict('records')
            elif file_extension.lower() == '.xlsx':
                df = pd.read_excel(queryfile)
                queryjson = df.to_dict('records')
            else:
                raise ValueError(f"Unsupported file type: {file_extension}")

            # Process each entry in the data
            for entry in queryjson:
                if entry.get('query') and entry.get('name'):
                    # Determine MS level based on query content
                    mslevel = 2 if "scaninfo(MS2DATA)" in entry['query'] else 1
                    entry.update({"mslevel": mslevel})
                    # Extract RT min and max values
                    rtmin_val = extract_rtmin_value(entry['query'])
                    rtmax_val = extract_rtmax_value(entry['query'])
                    entry.update({'rtmin': rtmin_val, 'rtmax': rtmax_val})
                    
                    # Update dictionaries with KEGG and group info
                    if entry.get('KEGG'):
                        name_kegg_dict.update({entry['name']: entry.get('KEGG')})
                    if entry.get('group'):
                        query_groups.update({entry['name']: entry.get('group')})
                    else:
                        entry.update({"group": None})
                    queries.append(entry)
                else:
                    sys.stdout.write(f"\n\nInvalid query:\n {entry}\n")
                    sys.stdout.flush()

            # Convert queries to DataFrame for MS1 and MS2
            query_df = pd.DataFrame(data=queries)
            ms1_query_df = query_df[query_df['mslevel'] == 1]
            ms2_query_df = query_df[query_df['mslevel'] == 2]
            sys.stdout.write(f"\nCreated {str(len(queries))} MassQL Queries from {queryfile}")
            sys.stdout.flush()
        except FileNotFoundError:
            sys.stdout.write(f"\n\nError: The file '{queryfile}' was not found.")
            sys.stdout.flush()
        except (json.JSONDecodeError, ValueError) as e:
            sys.stdout.write(f"\n\nError: The file '{queryfile}' contains invalid JSON. {str(e)}")
            sys.stdout.flush()
    else:
        sys.stdout.write("\nNo queryfile specified")
        sys.stdout.flush()

    return queries, ms1_query_df, ms2_query_df, query_groups, name_kegg_dict


"""Convert raw files"""

def convert_raw_files(convert_raw, msconvertexe, data_directory, convert_count = 0):
    if convert_raw and msconvertexe:
        try:
            subprocess.run(msconvertexe, 
                           stdout=subprocess.DEVNULL, 
                           stderr=subprocess.STDOUT, 
                           creationflags=subprocess.CREATE_NO_WINDOW)
            for fn in os.listdir(data_directory):
                if ".raw" in fn:
                    if os.path.isfile(data_directory + '\\' + fn.replace('.raw','.mzML')):
                        pass
                    else:
                        sys.stdout.write(f"\n{data_directory} \\ {fn} converting")
                        sys.stdout.flush()
                        subprocess.run(msconvertexe + " " + fn +  " --zlib", 
                                      stdout=subprocess.DEVNULL,
                                       stderr=subprocess.STDOUT,
                                       creationflags=subprocess.CREATE_NO_WINDOW)
                        if os.path.isfile(data_directory + '\\' + fn.replace('.raw','.mzML')):
                            sys.stdout.write(f", conversion complete!")
                            sys.stdout.flush()
                            convert_count += 1
                        else:
                            sys.stdout.write(f", conversion FAILED!")
                            sys.stdout.flush()
            if convert_count == 0:
                sys.stdout.write(f"\nNo files converted")
                sys.stdout.flush()
            elif convert_count > 0:
                sys.stdout.write(f"\n{str(convert_count)} file(s) converted")
                sys.stdout.flush()
        except Exception:
            sys.stdout.write(f"\nError. No raw files will be converted. Check path to MSConvert executable.")
            sys.stdout.flush()
    else:
        sys.stdout.write(f"\nNot converting raw files (if any)")
        sys.stdout.flush()


"""Verify files are present"""

def mzml_file_count(data_directory, file_count=0):
    try:
        file_count = len(fnmatch.filter(os.listdir(data_directory), '*.mzml'))
        if file_count == 0:
            sys.stdout.write(f"\n\nWarning: No mzml files found in {data_directory}")
            sys.stdout.flush()
            # input("Press enter to exit...")
            # exit()
        else:
            sys.stdout.write(f"\n{file_count} mzML files found in {data_directory}\n")
            sys.stdout.flush()
    except FileNotFoundError as e:
        sys.stdout.write(f"\nFileNotFoundError\n{e}\nNot found in {data_directory}")
        sys.stdout.flush()
        # input("Press enter to exit...")
        # exit()
    return file_count


"""Query files"""

def query_files(data_directory, queries, datasaver, scan_attributes = True, cache_setting=True):
     #scan_attributes: True/False. Use True if multiple collision parameters per file. Kills performance.
    timestr = time.strftime("%Y_%m_%d_%H%M")

    #Function for getting collision energy information from scan
    def lookup_scan(row):
        try:
            collision_type = list(mzml_reader[row['scan']-1]['precursorList']['precursor'][0]['activation'].keys())[0]
            energy = str(float(list(mzml_reader[row['scan']-1]['precursorList']['precursor'][0]['activation'].values())[1]))
        except:
            collision_type = energy = "None"
        return collision_type, energy
    
    raw_df_ms1 = pd.DataFrame()
    raw_df_ms2 = pd.DataFrame()
    raw_ms1_df_list = []
    raw_ms2_df_list = []
    filename_groups = {}
    
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    counter = 0
    fail_log = []

    filename_list = sorted(glob.iglob(os.path.join(data_directory, '*.mzML')))
    filename_list = sorted([os.path.normpath(path) for path in filename_list])
    # filename_list = sorted([path.replace('\\', '/') for path in filename_list])

    if queries:
        # filename_list = filename_list[::-1]
        # filename_list = filename_list[0:4] #subset files
        sys.stdout.write(f"\nApplying {len(queries)} queries to {len(filename_list)} files")
        sys.stdout.flush()
        for filename in filename_list:
            ms1_df = pd.DataFrame()
            ms2_df = pd.DataFrame()
            filename_base = os.path.basename(filename)
            fail_count = 0
            counter += 1
            sys.stdout.write(f"\nProcessing file {counter}: {filename_base}")
            sys.stdout.flush()
            ms1_df, ms2_df = msql_fileloading.load_data(filename, cache=cache_setting)
            if scan_attributes: 
                mzml_reader = mzml.MzML(filename)
            for i, query in enumerate(queries):
                i = i + 1
                # sys.stdout.write(f"\rApplying {i}/{len(queries)} MassQL Queries to File {counter} of {len(filename_list)} ({fail_count} Failed)")
                # sys.stdout.flush()
        
                # scannum_query = query['query'].replace("scaninfo", "scannum")
                results_df = pd.DataFrame()
                try:
                    with contextlib.redirect_stdout(stdout_buffer), contextlib.redirect_stderr(stderr_buffer):
                        results_df = msql_engine.process_query(query['query'], filename, cache=cache_setting, ms1_df=ms1_df, ms2_df=ms2_df)
                    if not results_df.empty:
                        results_df['filename'] = filename_base
                        results_df['query_name'] = query['name']
                        results_df['query'] = query['query']
                        
                        if scan_attributes:
                            if query['mslevel'] == 2:
                                results_df["collision_type"], results_df["energy"] = zip(*results_df.apply(lambda row: lookup_scan(row), axis=1))
                            else:
                                results_df["collision_type"], results_df["energy"] = 'NULL', 'NULL'
                    
                except Exception:
                    fail_count += 1
                    fail_log.append(str(query['name'] + " for " +str(filename_base)))
                    results_df = pd.DataFrame()
                    pass 

                
                if query['mslevel'] == 1:
                    raw_ms1_df_list.append(results_df)
                if query['mslevel'] == 2:
                    raw_ms2_df_list.append(results_df)

            if fail_count:
                sys.stdout.write(f"\n {fail_count} queries failed for file {counter}")
                sys.stdout.flush()
        
        if fail_log:
            fail_log = '\n'.join(fail_log)
            sys.stdout.write(f"\n\nQueries Failed:\n{fail_log}\n") 
            sys.stdout.flush()   
        
        results_df = pd.DataFrame()
        if raw_ms1_df_list:
            if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/scans/"):
                os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/scans/")
            raw_df_ms1 = pd.concat(raw_ms1_df_list)
            
        if raw_ms2_df_list:
            if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/scans/"):
                os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/scans/")
            raw_df_ms2 = pd.concat(raw_ms2_df_list)
            if not raw_df_ms2.empty:
                raw_df_ms2['collision_type'] = raw_df_ms2['collision_type'].apply(lambda x: 'CID' if x == 'collision-induced dissociation' else x)
                raw_df_ms2['collision_type'] = raw_df_ms2['collision_type'].apply(lambda x: 'HCD' if x == 'beam-type collision-induced dissociation' else x)
                raw_df_ms2['collision_type_energy'] = raw_df_ms2['collision_type'].astype(str) + "__" + raw_df_ms2['energy'].astype(str)
            
        if not raw_df_ms1.empty:
            raw_df_ms1.to_csv(data_directory + "/MassQLab_Output/"+timestr+"/ms1_raw_df.csv")
            sys.stdout.write(f"\nCreated raw_df_ms1 and exported as csv.") 
            sys.stdout.flush()
        if not raw_df_ms2.empty:
            raw_df_ms2.to_csv(data_directory + "/MassQLab_Output/"+timestr+"/ms2_raw_df.csv")
            sys.stdout.write(f"\nCreated raw_df_ms2 and exported as csv.") 
            sys.stdout.flush()

    return raw_df_ms1, raw_df_ms2, filename_groups, timestr

