import sys
# sys.stdout.write(f"--------------------------------\n   Initialize MassQLab   \n--------------------------------\n")
# sys.stdout.flush()
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
    # Initialize config as an empty dictionary in case the file doesn't exist or is invalid
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
    metadata_file = config.get('metadata_file', None)
    metadata_filename_column = config.get('metadata_filename_column', "CORE_Filename")
    metadata_group_columns = config.get('metadata_group_columns', "USER_Column1")
    kegg_path = config.get('kegg_path', None)
    convert_raw = config.get('convert_raw', False)
    msconvertexe = config.get('msconvert_exe', None)
    cache_setting = config.get('cache', True)
    datasaver = config.get('datasaver', False)

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
    return data_directory_loc, queryfile_loc, metadata_file, metadata_filename_column, metadata_group_columns, kegg_path, convert_raw, msconvertexe, cache_setting, datasaver
    





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


import numpy as np
import pandas as pd
from scipy.interpolate import interp1d
from scipy.integrate import trapz
from scipy.signal import find_peaks, peak_widths
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import sys

def analysis_ms1(raw_df_ms1, data_directory, timestr):

    def index_to_xdata(xdata, indices):
        "interpolate the values from signal.peak_widths to xdata"
        ind = np.arange(len(xdata))
        f = interp1d(ind, xdata)
        return f(indices)

    def custom_score(peak):
        normalized_height = ydata_all[peak] / max(ydata_all)  # Normalize height by the maximum height in the signal
        height_score = normalized_height  
        closeness_score = 1 / ((abs(xdata_all[peak] - avg_rt))**(0.5))
        score = height_score  # Update this line if you want to use closeness_score as well
        return score

    def gaussian(x, A, mu, sigma):
        return A * np.exp(-((x - mu)**2) / (2 * sigma**2))

    def generate_linear_array(rtmin, rtmax):
        start = rtmin - (rtmax - rtmin)
        end = rtmax + (rtmax - rtmin)
        points = np.linspace(start, end, 1000)
        return points

    ms1_analysis_df_list = []
    plots = {}  # Dictionary to store individual plot file paths by (filename, query)

    if not raw_df_ms1.empty:
        unique_filenames = raw_df_ms1['filename'].unique()
        unique_queries = raw_df_ms1['query'].unique()
        num_files = len(unique_filenames)
        num_queries = len(unique_queries)

        output_dir = data_directory + "/MassQLab_Output/" + timestr + "/ms1_traces/"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        with PdfPages(data_directory + "/MassQLab_Output/" + timestr + "/ms1_traces.pdf") as pdf:
            for group_name, grouped_df in raw_df_ms1.groupby(['filename', 'query_name', 'query']):
                if len(grouped_df) < 5:
                    continue

                fig, ax = plt.subplots()  # Create a figure and axis for each plot

                analysis_df_dict = {'filename': {0: group_name[0]},
                                    'query_name': {0: group_name[1]},
                                    'query': {0: group_name[2]},
                                    'peak_area': {0: None},
                                    'peak_area_alt': {0: None},
                                    'fwhm': {0: None},
                                    'rtmin': {0: None},
                                    'rtmax': {0: None},
                                    'measured_RT': {0: None},
                                    'rt_error': {0: None},
                                    'datapoint_count': {0: len(grouped_df)}}

                rtmin = extract_rtmin_value(group_name[2])
                rtmax = extract_rtmax_value(group_name[2])
                avg_rt = (rtmin + rtmax) / 2
                range_rt = rtmax - rtmin
                xdata_all = np.array(grouped_df.rt)
                ydata_all = np.array(grouped_df.i)
                height_min = abs(ydata_all.mean() * 1.1)
                prominence_min = abs(np.max(ydata_all) / 10)

                peaks, properties = find_peaks(ydata_all, height=height_min, distance=5, prominence=prominence_min)
                filtered_peaks = [peak for peak in peaks if (rtmin - range_rt) < xdata_all[peak] < (rtmax + range_rt)]
                
                widths, width_heights, left_ips, right_ips = peak_widths(ydata_all, filtered_peaks, rel_height=0.5)
                xleft_ips = index_to_xdata(xdata_all, left_ips)
                xright_ips = index_to_xdata(xdata_all, right_ips)
                real_width = xright_ips - xleft_ips

                final_peaks = [peak for i, peak in enumerate(filtered_peaks) if real_width[i] < range_rt]
                top_n_peaks_prominence = sorted(final_peaks, key=lambda x: ydata_all[x], reverse=True)[:2]
                top_n_peaks = sorted(top_n_peaks_prominence, key=custom_score, reverse=True)
                
                peak_x_values = xdata_all[top_n_peaks]
                peak_y_values = ydata_all[top_n_peaks]

                if peak_x_values.any():
                    peak_x = peak_x_values[0]
                    peak_y = peak_y_values[0]
                    Awidths, Awidth_heights, Aleft_ips, Aright_ips = peak_widths(ydata_all, top_n_peaks, rel_height=0.5)
                    Aleft_ips = index_to_xdata(xdata_all, Aleft_ips)
                    Aright_ips = index_to_xdata(xdata_all, Aright_ips)
                    Areal_width = Aright_ips - Aleft_ips
                    Amidpoint = (Aright_ips + Aleft_ips) / 2

                    fwhm = Areal_width[0]
                    measured_RT = Amidpoint[0]
                    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
                    interpolated_array = generate_linear_array(rtmin, rtmax)
                    y_values = gaussian(interpolated_array, peak_y, Amidpoint[0], sigma)

                    ax.plot(interpolated_array, y_values, label='Gaussian Peak', color='r', lw=2)
                    peak_area = peak_y * sigma * np.sqrt(2 * np.pi)
                    peak_area_alt = trapz(grouped_df.i, x=grouped_df.rt)
                    rt_error = 0 if rtmin <= measured_RT <= rtmax else measured_RT - rtmin if measured_RT < rtmin else measured_RT - rtmax if measured_RT > rtmax else None
                    ax.hlines(Awidth_heights, Aleft_ips, Aright_ips, color='blue')
                    
                    analysis_df_dict.update({'peak_area': {0: peak_area},
                                             'peak_area_alt': {0: peak_area_alt},
                                             'fwhm': {0: fwhm},
                                             'rtmin': {0: rtmin},
                                             'rtmax': {0: rtmax},
                                             'measured_RT': {0: measured_RT},
                                             'rt_error': {0: rt_error}})

                ax.plot(xdata_all, ydata_all, '-o', label='Data', color='black', markerfacecolor='black', markeredgecolor='black')
                plot_symbols = ['x', '+', '+', '+', '+']
                for peak_index, peak_num in enumerate(top_n_peaks):
                    ax.plot(xdata_all[peak_num], ydata_all[peak_num], str(plot_symbols[peak_index]), color='r')
                ax.set_xlim(xmin=(rtmin - range_rt), xmax=(rtmax + range_rt))
                ax.axvline(x=rtmin, color='black', linestyle=':')
                ax.axvline(x=rtmax, color='black', linestyle=':')
                ax.set_xlabel('rt')
                ax.set_ylabel('intensity')
                plt_title = str(group_name[0]) + "_" + str(group_name[1])
                ax.set_title(str(group_name[0]) + '\n' + str(group_name[1]))
                ax.legend()
                pdf.savefig(fig)
                
                plot_path = output_dir + plt_title + '.png'
                fig.savefig(plot_path)
                plt.close(fig)

                if os.path.exists(plot_path):  # Check if the plot file was created
                    plots[(group_name[0], group_name[2])] = plot_path  # Save the plot file path

                analysis_df = pd.DataFrame(data=analysis_df_dict)
                ms1_analysis_df_list.append(analysis_df)

    if ms1_analysis_df_list:
        ms1_analysis_df = pd.concat(ms1_analysis_df_list)
        # ms1_analysis_df.to_csv("MassQLab_Output/"+timestr+"/"+timestr+"_ms1_analysis_df.csv")

    # Create consolidated grid of plots
    fig, axes = plt.subplots(num_files, num_queries, figsize=(6 * num_queries, 5 * num_files))
    
    for i, filename in enumerate(unique_filenames):
        for j, query in enumerate(unique_queries):
            if (filename, query) in plots:
                plot_path = plots[(filename, query)]
                img = plt.imread(plot_path)
                axes[i, j].imshow(img)
                axes[i, j].axis('off')  # Hide the axes
            else:
                axes[i, j].text(0.5, 0.5, 'No Data', ha='center', va='center', color='red')
                axes[i, j].axis('off')  # Hide the axes

    plt.tight_layout()
    consolidated_plot_path = data_directory + "/MassQLab_Output/" + timestr + "/consolidated_ms1_traces.png"
    plt.savefig(consolidated_plot_path)
    plt.close()

    sys.stdout.write(f"\nCreated ms1 analysis dataframe") 
    sys.stdout.flush()
    return ms1_analysis_df



"""merge ms1 query dataframe with ms1 analysis dataframe"""

def ms1_query_analysis_merge(ms1_analysis_df, ms1_query_df):
    if not ms1_analysis_df.empty and not ms1_query_df.empty:
        ms1_analysis_df = pd.merge(ms1_analysis_df, ms1_query_df.rename(columns={'name': 'query_name'}), on='query_name', how='inner', suffixes=('', '_duplicate'))
        ms1_analysis_df = ms1_analysis_df.drop(columns=[col for col in ms1_analysis_df.columns if "_duplicate" in col])
    return ms1_analysis_df




"""MS1 Peak Validity and QC Check"""

def ms1_validity_and_QC(ms1_analysis_df, fwhm_thresh=1, rt_thresh=0.1, datapoint_thresh=5, abundance_thresh_ms1=10, impute_ms1_data=True):
    
    if not ms1_analysis_df.empty:
        # Define peak validity conditions
        peak_condition = (
            (ms1_analysis_df['fwhm'] < fwhm_thresh) &
            (ms1_analysis_df.apply(
                lambda row: (np.abs(row['rt_error']) <= (row['rt_tolerance'] if 'rt_tolerance' in ms1_analysis_df.columns and row['rt_tolerance'] is not None else rt_thresh))
                if row['rt_error'] is not None else False,
                axis=1
            )) &
            (ms1_analysis_df['datapoint_count'] >= datapoint_thresh)
        )
        ms1_analysis_df['peak_valid'] = np.where(peak_condition, True, False)

        # Abundance check
        if 'abundance' in ms1_analysis_df.columns:
            ms1_analysis_df['abundance_error'] = np.where(
                (ms1_analysis_df['peak_area'].notnull()) &
                (ms1_analysis_df['abundance'].notnull()) &
                (ms1_analysis_df['abundance'] != 0),
                ((ms1_analysis_df['peak_area'] - ms1_analysis_df['abundance']) / ms1_analysis_df['abundance']) * 100,
                None
            )
            abundance_condition = ms1_analysis_df.apply(
                lambda row: np.abs(row['abundance_error']) <= (row['abundance_tolerance'] if 'abundance_tolerance' in ms1_analysis_df.columns and row['abundance_tolerance'] is not None else abundance_thresh_ms1)
                if row['abundance_error'] is not None else False,
                axis=1
            )
            ms1_analysis_df['abundance_valid'] = np.where(abundance_condition, True, False)
        else:
            ms1_analysis_df['abundance_error'] = ms1_analysis_df['abundance_valid'] = None

    return ms1_analysis_df





"""Impute MS1 function"""

def impute_ms1(ms1_analysis_df, c1='query_name', c3='filename', c4='peak_valid', c5='abundance_valid'):
    result_xdf = pd.DataFrame()
    if not ms1_analysis_df.empty:
        # Get unique values of c1 and c2
        unique_c1 = ms1_analysis_df[c1].unique()
        
        # Generate all combinations of c1, c2, and c3
        all_combinations = list(product(unique_c1, ms1_analysis_df[c3].unique()))
        
        # Create a DataFrame with all combinations
        result_xdf = pd.DataFrame(all_combinations, columns=[c1, c3])
        
        # Merge with the original DataFrame to get the 'Value' column
        result_xdf = pd.merge(result_xdf, ms1_analysis_df, on=[c1, c3], how='left', indicator=True)
        # Create 'imputed' column based on the '_merge' indicator
        result_xdf['imputed'] = result_xdf['_merge'].eq('left_only')
        
        # Drop the '_merge' column
        result_xdf.drop('_merge', axis=1, inplace=True)
        result_xdf[c4] = result_xdf[c4].fillna(False)
        result_xdf[c5] = result_xdf[c5].fillna(False)
    
    return result_xdf


"""Export MS1 analysis dataframe"""

def export_ms1_analysis_df(ms1_analysis_df, data_directory, timestr):
    if not ms1_analysis_df.empty:
        ms1_analysis_df.to_csv(data_directory + "/MassQLab_Output/"+timestr+"/ms1_analysis_df.csv")
        sys.stdout.write(f"\nCreated ms1_analysis_df and exported as csv.") 
        sys.stdout.flush()
    else:
        sys.stdout.write(f"\nms1_analysis_df empty. Not writing to csv.")
        sys.stdout.flush()
        

"""RT Analysis MS1"""

def rt_analysis_ms1(ms1_analysis_df, data_directory, timestr):
    if not ms1_analysis_df.empty:
        df = ms1_analysis_df.copy()
        df.loc[~df['peak_valid'], 'measured_RT'] = None
        pivot_df = df.pivot(index='query_name', columns='filename', values='measured_RT')
        original_columns = pivot_df.columns.tolist()
        
        # Calculate the average of each row based on original 'measured_RT' values
        pivot_df['average'] = pivot_df[original_columns].mean(axis=1)
        
        # Calculate the standard deviation of each row based on original 'measured_RT' values
        pivot_df['std_dev'] = pivot_df[original_columns].std(axis=1)
        
        # Add a new column with the RSD percentage of each row
        pivot_df['RSD_%'] = (pivot_df['std_dev'] / pivot_df['average']) * 100
        
        # Add a new column 'valid_%' with the percentage of valid (non-None) values
        pivot_df['valid_%'] = (pivot_df[original_columns].count(axis=1) / len(original_columns)) * 100
        
        if not pivot_df.empty:
            pivot_df.to_csv(data_directory + "/MassQLab_Output/"+timestr+"/ms1_RT_analysis_df.csv")
            sys.stdout.write(f"\nCreated ms1_RT_analysis_df and exported as csv.") 
            sys.stdout.flush()
            
    else:
        sys.stdout.write(f"\nms1_analysis_df empty. Not doing rt_analysis_ms1.")
        sys.stdout.flush()



"""Summary MS1 traces"""

def summary_ms1_traces(raw_df_ms1, data_directory, timestr):
    if not raw_df_ms1.empty:
        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr)
        pdf_path = os.path.join(output_dir, "ms1_summary_traces.pdf")
        all_data_dir = os.path.join(output_dir, "ms1_summary_traces", "all_data")
        rt_range_dir = os.path.join(output_dir, "ms1_summary_traces", "rt_range")
        
        if not os.path.exists(all_data_dir):
            os.makedirs(all_data_dir)
        if not os.path.exists(rt_range_dir):
            os.makedirs(rt_range_dir)

        with PdfPages(pdf_path) as pdf:
            # Iterate through unique 'query_name' groups
            for i, (frame_group, frame_data) in enumerate(raw_df_ms1.groupby(['query_name', 'query'])):
                rtmin = extract_rtmin_value(frame_group[1])
                rtmax = extract_rtmax_value(frame_group[1])
                
                plt_title = frame_group[0]
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                
                # Create the first figure
                fig1, ax1 = plt.subplots(figsize=(10, 5))
                ax1.set_title(f'{plt_title} - All Data')
                
                # Iterate through unique 'filename' groups within each 'query_name' group
                for group_name, group_data in frame_data.groupby('filename'):
                    c = next(colors_plot)["color"]
                    # Plot 'rt' versus 'i' for each 'filename' group within the 'query_name' group
                    group_data.plot(x='rt', y='i', kind='line', label=group_name, ax=ax1, color=c, style='.-')
                
                ax1.set_ylabel('intensity')
                # Place the legend outside the plot area
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig1.tight_layout(rect=[0, 0, 0.85, 1])
                pdf.savefig(fig1, bbox_inches='tight')
                fig1.savefig(os.path.join(all_data_dir, f"{plt_title}.png"), bbox_inches='tight')
                plt.close(fig1)

                # Create the second figure
                fig2, ax2 = plt.subplots(figsize=(10, 5))
                ax2.set_title(f'{plt_title} - RT Range')

                # Iterate through unique 'filename' groups within each 'query_name' group
                for group_name, group_data in frame_data.groupby('filename'):
                    c = next(colors_plot)["color"]
                    group_data.plot(x='rt', y='i', kind='line', label=group_name, ax=ax2, color=c, style='.-')
                
                range_rt = rtmax - rtmin
                ax2.set_xlim(rtmin - range_rt, rtmax + range_rt)
                ax2.axvline(x=rtmin, color='black', linestyle=':')
                ax2.axvline(x=rtmax, color='black', linestyle=':')

                ax2.set_ylabel('intensity')
                # Place the legend outside the plot area
                ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig2.tight_layout(rect=[0, 0, 0.85, 1])
                pdf.savefig(fig2, bbox_inches='tight')
                fig2.savefig(os.path.join(rt_range_dir, f"{plt_title}.png"), bbox_inches='tight')
                plt.close(fig2)

        sys.stdout.write(f"\nCreated ms1_summary_traces.")
        sys.stdout.flush()





"""Summary MS1 traces inverse"""

def summary_ms1_traces_inverse(raw_df_ms1, data_directory, timestr):
    if not raw_df_ms1.empty:
        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr)
        pdf_path = os.path.join(output_dir, "ms1_summary_traces_inverse.pdf")
        all_data_dir = os.path.join(output_dir, "ms1_summary_traces_inverse", "all_data")
        rt_range_dir = os.path.join(output_dir, "ms1_summary_traces_inverse", "rt_range")
        
        if not os.path.exists(all_data_dir):
            os.makedirs(all_data_dir)
        if not os.path.exists(rt_range_dir):
            os.makedirs(rt_range_dir)

        with PdfPages(pdf_path) as pdf:
            # Iterate through unique 'filename' groups
            for i, (frame_group, frame_data) in enumerate(raw_df_ms1.groupby('filename')):
                
                rtmins = []
                rtmaxs = []
                plt_title = frame_group
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                
                # Create the first figure
                fig1, ax1 = plt.subplots(figsize=(10, 5))
                ax1.set_title(f'{plt_title} - All Data')
                
                # Iterate through unique 'query_name' groups within each 'filename' group
                for group_name, group_data in frame_data.groupby(['query_name', 'query']):
                    rtmins.append(extract_rtmin_value(group_name[1]))
                    rtmaxs.append(extract_rtmax_value(group_name[1]))
                    c = next(colors_plot)["color"]
                    # Plot 'rt' versus 'i' for each 'query_name' group within the 'filename' group
                    group_data.plot(x='rt', y='i', kind='line', label=group_name[0], ax=ax1, color=c, style='.-')
                
                ax1.set_ylabel('intensity')
                # Place the legend outside the plot area
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig1.tight_layout(rect=[0, 0, 0.85, 1])
                pdf.savefig(fig1, bbox_inches='tight')
                fig1.savefig(os.path.join(all_data_dir, f"{plt_title}_all_data.png"), bbox_inches='tight')
                plt.close(fig1)

                # Create the second figure
                fig2, ax2 = plt.subplots(figsize=(10, 5))
                ax2.set_title(f'{plt_title} - RT Range')

                # Iterate through unique 'query_name' groups within each 'filename' group
                for group_name, group_data in frame_data.groupby(['query_name', 'query']):
                    c = next(colors_plot)["color"]
                    group_data.plot(x='rt', y='i', kind='line', label=group_name[0], ax=ax2, color=c, style='.-')
                
                if rtmins and rtmaxs:
                    ax2.set_xlim(min(rtmins)-1, max(rtmaxs)+1)
                
                ax2.set_ylabel('intensity')
                # Place the legend outside the plot area
                ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig2.tight_layout(rect=[0, 0, 0.85, 1])
                pdf.savefig(fig2, bbox_inches='tight')
                fig2.savefig(os.path.join(rt_range_dir, f"{plt_title}_rt_range.png"), bbox_inches='tight')
                plt.close(fig2)

        sys.stdout.write(f"\nCreated ms1_summary_traces_inverse.")
        sys.stdout.flush()


"""Summary MS1 Areas"""

def summary_ms1_areas(ms1_analysis_df, data_directory, timestr, abundance_thresh_ms1=10):
    if not ms1_analysis_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas.pdf") as pdf:
            for i, (frame_group, frame_data) in enumerate(ms1_analysis_df.groupby(['query_name'])):
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                
                frame_data.loc[frame_data['peak_valid'] == False, 'peak_area'] = np.nan
                
                plt.figure(figsize=(5, 5))
                plt_title = frame_group[0]
                plt.title(plt_title)
                frame_data['filename_short'] = '..' + frame_data['filename'].str[-15:-5]
                frame_data['filename_short_Index'] = range(1, len(frame_data) + 1)
                
                # Iterate through unique 'filename' groups within each 'query_name' group
                for group_name, group_data in frame_data.groupby('filename'):
                    c = next(colors_plot)["color"]
                    plt.bar(group_data['filename_short_Index'], group_data['peak_area'], color='gray')
        
                plt.xticks(frame_data['filename_short_Index'])
                plt.gca().set_xticklabels(frame_data['filename_short'])
                plt.xlim(frame_data['filename_short_Index'].min() - 1, frame_data['filename_short_Index'].max() + 1)
        
                if 'abundance' in frame_data.columns:
                    frame_abundance = frame_data['abundance'].mean()
                    
                    # Dynamically calculate tolerance
                    tolerance = frame_data['abundance_tolerance'].mean() if 'abundance_tolerance' in frame_data.columns and frame_data['abundance_tolerance'].notnull().any() else abundance_thresh_ms1
                    
                    plt.axhline(y=(frame_abundance + (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
                    plt.axhline(y=(frame_abundance - (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
                
                plt.tick_params(axis='x', rotation=90)
                plt.ylabel('peak area')
                plt.tight_layout(rect=[0, 0, 1, 0.96])
        
                pdf.savefig()
                if not os.path.exists(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas/"):
                    os.makedirs(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas/")
                plt.savefig(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas/" + plt_title + '.png')
                plt.close('all')
        
        sys.stdout.write(f"\nCreated ms1_summary_areas.") 
        sys.stdout.flush()


"""Summary MS1 Areas system test"""

def summary_ms1_areas_system_test(ms1_analysis_df, data_directory, timestr, abundance_thresh_ms1=10):
    system_test_df = ms1_analysis_df[ms1_analysis_df['filename'].str.startswith("system_test", na=False)]
    if not system_test_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas_system_test.pdf") as pdf:
            for i, (frame_group, frame_data) in enumerate(system_test_df.groupby(['query_name'])):
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                
                frame_data.loc[frame_data['peak_valid'] == False, 'peak_area'] = np.nan
                
                plt.figure(figsize=(5, 5))
                plt_title = frame_group[0]
                plt.title(plt_title)
                frame_data['filename_short'] = '..' + frame_data['filename'].str[-15:-5]
                frame_data['filename_short_Index'] = range(1, len(frame_data) + 1)
                
                # Iterate through unique 'filename' groups within each 'query_name' group
                for group_name, group_data in frame_data.groupby('filename'):
                    c = next(colors_plot)["color"]
                    plt.bar(group_data['filename_short_Index'], group_data['peak_area'], color='gray')
        
                plt.xticks(frame_data['filename_short_Index'])
                plt.gca().set_xticklabels(frame_data['filename_short'])
                plt.xlim(frame_data['filename_short_Index'].min() - 1, frame_data['filename_short_Index'].max() + 1)
        
                if 'abundance' in frame_data.columns:
                    frame_abundance = frame_data['abundance'].mean()
                    
                    # Dynamically calculate tolerance
                    tolerance = frame_data['abundance_tolerance'].mean() if 'abundance_tolerance' in frame_data.columns and frame_data['abundance_tolerance'].notnull().any() else abundance_thresh_ms1
                    
                    plt.axhline(y=(frame_abundance + (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
                    plt.axhline(y=(frame_abundance - (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
                
                plt.tick_params(axis='x', rotation=90)
                plt.ylabel('peak area')
                plt.tight_layout(rect=[0, 0, 1, 0.96])
        
                pdf.savefig()
                if not os.path.exists(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas_system_test/"):
                    os.makedirs(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas_system_test/")
                plt.savefig(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_areas_system_test/" + plt_title + '.png')
                plt.close('all')
        
        sys.stdout.write(f"\nCreated ms1_summary_areas_system_test.") 
        sys.stdout.flush()

"""Summary MS1 Areas Inverse"""

def summary_ms1_areas_inverse(ms1_analysis_df, data_directory, timestr, abundance_thresh_ms1 = 10):
    if not ms1_analysis_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_areas_inverse.pdf") as pdf:
            for i, (frame_group, frame_data) in enumerate(ms1_analysis_df.groupby(['filename'])):
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                
                frame_data.loc[frame_data['peak_valid'] == False, 'peak_area'] = np.nan
                
                plt.figure(figsize=(5, 5))
                plt_title = frame_group[0]
                plt.title(plt_title)
                frame_data['queryname_short'] = '..' + frame_data['query_name'].str[-15:-5]
                frame_data['queryname_short_Index'] = range(1, len(frame_data) + 1)
                # Iterate through unique 'filename' groups within each 'query_name' group
                for group_name, group_data in frame_data.groupby('query_name'):
                    c = next(colors_plot)["color"]
                    # Plot 'rt' versus 'i' for each 'filename' group within the 'query_name' group
                    
                    plt.bar(group_data['queryname_short_Index'], group_data['peak_area'], color='gray')
        
                plt.xticks(frame_data['queryname_short_Index'])
                plt.gca().set_xticklabels(frame_data['queryname_short'])
                plt.xlim(frame_data['queryname_short_Index'].min() - 1, frame_data['queryname_short_Index'].max() + 1)
        
                
                # if 'abundance' in frame_data.columns:
                #     frame_abundance = frame_data['abundance'].mean()
                #     plt.axhline(y=(frame_abundance+(frame_abundance*(abundance_thresh_ms1/100))), color='black', linestyle='--')
                #     plt.axhline(y=(frame_abundance-(frame_abundance*(abundance_thresh_ms1/100))), color='black', linestyle='--')
                plt.tick_params(axis='x', rotation=90)
                plt.ylabel('peak area')
                plt.tight_layout(rect=[0, 0, 1, 0.96])
        
                pdf.savefig()
                if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_areas_inverse/"):
                    os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_areas_inverse/")
                plt.savefig(data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_areas_inverse/"+plt_title+'.png')
                # plt.show()
                plt.close('all')
        sys.stdout.write(f"\nCreated ms1_summary_areas_inverse.") 
        sys.stdout.flush()


# """Summary MS1 RT system test"""

# def summary_ms1_rt_system_test(ms1_analysis_df, data_directory, timestr, abundance_thresh_ms1=10):
#     system_test_df = ms1_analysis_df[ms1_analysis_df['filename'].str.startswith("system_test", na=False)]
#     if not system_test_df.empty:
#         with PdfPages(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_rt_system_test.pdf") as pdf:
#             for i, (frame_group, frame_data) in enumerate(system_test_df.groupby(['query_name'])):
#                 colors_plot = plt.rcParams["axes.prop_cycle"]()
                
#                 frame_data.loc[frame_data['peak_valid'] == False, 'peak_area'] = np.nan
                
#                 plt.figure(figsize=(5, 5))
#                 plt_title = frame_group[0]
#                 plt.title(plt_title)
#                 frame_data['filename_short'] = '..' + frame_data['filename'].str[-15:-5]
#                 frame_data['filename_short_Index'] = range(1, len(frame_data) + 1)
                
#                 # Iterate through unique 'filename' groups within each 'query_name' group
#                 for group_name, group_data in frame_data.groupby('filename'):
#                     c = next(colors_plot)["color"]
#                     plt.bar(group_data['filename_short_Index'], group_data['measured_RT'], color='gray')
        
#                 plt.xticks(frame_data['filename_short_Index'])
#                 plt.gca().set_xticklabels(frame_data['filename_short'])
#                 plt.xlim(frame_data['filename_short_Index'].min() - 1, frame_data['filename_short_Index'].max() + 1)
        
#                 if 'abundance' in frame_data.columns:
#                     frame_abundance = frame_data['abundance'].mean()
                    
#                     # Dynamically calculate tolerance
#                     tolerance = frame_data['abundance_tolerance'].mean() if 'abundance_tolerance' in frame_data.columns and frame_data['abundance_tolerance'].notnull().any() else abundance_thresh_ms1
                    
#                     plt.axhline(y=(frame_abundance + (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
#                     plt.axhline(y=(frame_abundance - (frame_abundance * (tolerance / 100))), color='black', linestyle='--')
                
#                 plt.tick_params(axis='x', rotation=90)
#                 plt.ylabel('peak area')
#                 plt.tight_layout(rect=[0, 0, 1, 0.96])
        
#                 pdf.savefig()
#                 if not os.path.exists(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_rt_system_test/"):
#                     os.makedirs(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_rt_system_test/")
#                 plt.savefig(data_directory + "/MassQLab_Output/" + timestr + "/ms1_summary_rt_system_test/" + plt_title + '.png')
#                 plt.close('all')
        
#         sys.stdout.write(f"\nCreated ms1_summary_rt_system_test.") 
#         sys.stdout.flush()
        

def plot_ms2(raw_df_ms2, data_directory, timestr):
    if not raw_df_ms2.empty:
        with PdfPages(data_directory + "/MassQLab_Output/"+timestr+"/ms2_plots.pdf") as pdf:
            
            raw_df_ms2['collision_type'] = raw_df_ms2['collision_type'].apply(lambda x: 'CID' if x == 'collision-induced dissociation' else x)
            raw_df_ms2['collision_type'] = raw_df_ms2['collision_type'].apply(lambda x: 'HCD' if x == 'beam-type collision-induced dissociation' else x)
            raw_df_ms2['collision_type_energy'] = raw_df_ms2['collision_type'].astype(str) + "__" + raw_df_ms2['energy'].astype(str)
        
            for group_name_x, grouped_df_x in raw_df_ms2.groupby(['filename', 'query_name', 'query']):
                plt_title = str(group_name_x[0]) + "_" + str(group_name_x[1])
                colors = plt.cm.get_cmap('tab10')  # Example using a colormap; adjust accordingly
                color_dict = {}  # Dictionary to store assigned colors
                fig, ax = plt.subplots()
                ax.set_title(plt_title)

                for group_name, grouped_df in grouped_df_x.groupby('collision_type_energy'):
                    color = colors(len(color_dict))  # Assign color from colormap
                    color_dict[group_name] = color  # Store color in dictionary for consistency
                    
                    # Find the highest intensity point in the entire grouped_df
                    highest_point_idx = grouped_df['i'].idxmax()
                    highest_point = grouped_df.loc[highest_point_idx]
                    
                    label_added = False  # Track if label has been added to legend for this group
                    for group_name2, grouped_df2 in grouped_df.groupby('ms1scan'):
                        if not grouped_df2.empty:
                            # Plot each grouped_df2 with open circles. Only add label for the first subplot of each group_name to avoid duplicate legend entries
                            if not label_added:
                                ax.plot(grouped_df2.rt, grouped_df2.i, '-o', label=group_name, color=color_dict[group_name])
                                label_added = True
                            else:
                                ax.plot(grouped_df2.rt, grouped_df2.i, '-o', color=color_dict[group_name])
                            
                            # Check if the highest point is in this subgroup and mark it
                            if highest_point_idx in grouped_df2.index:
                                ax.plot(highest_point.rt, highest_point.i, 'o', markersize=12, color=color_dict[group_name])  # Mark the highest point distinctly as a solid circle
                                ax.text(highest_point.rt, highest_point.i, highest_point.scan, color='black', fontsize=10)  # Label the highest point with scan value
            
                ax.legend(title='CTE')
                pdf.savefig()
                if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms2_plots/"):
                    os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms2_plots/")
                plt.savefig(data_directory + "/MassQLab_Output/"+timestr+"/ms2_plots/"+plt_title+'.png')
                # plt.show()
                plt.close('all')
    sys.stdout.write(f"\nCreated ms2 plots.")





"""Analysis MS2"""

def analysis_ms2(raw_df_ms2, ms2_query_df, abundance_thresh_ms2 = 20):
    ms2_analysis_df = pd.DataFrame()
    ms2_analysis_df_list = []
    
    if not raw_df_ms2.empty:
        for group_name, grouped_df in raw_df_ms2.groupby(['filename', 'query_name', 'query', 'collision_type_energy']):
            grouped_df['i'].fillna(0, inplace=True)
    
            if not grouped_df.empty and len(grouped_df) > 0:
                
                grouped_df = grouped_df.loc[grouped_df['i'].idxmax()].to_frame().T
                grouped_df['i'] = pd.to_numeric(grouped_df['i'])
                grouped_df['measured_RT'] = grouped_df.loc[grouped_df['i'].idxmax()]['rt']
                grouped_df['scan'] = grouped_df.iloc[0]["scan"]
                grouped_df['collision_type_energy'] = grouped_df.iloc[0]["collision_type_energy"]
    
                ms2_analysis_df_list.append(grouped_df)
                    
            else:
                #TODO check what shows up here
                i = measured_RT = scan = collision_type_energy  = None
    
    if ms2_analysis_df_list:
        ms2_analysis_df = pd.concat(ms2_analysis_df_list).reset_index(drop=True)

    if not ms2_analysis_df.empty and not ms2_query_df.empty:
        ms2_analysis_df = pd.merge(ms2_analysis_df, ms2_query_df.rename(columns={'name': 'query_name'}), on='query_name', how='inner', suffixes=('', '_duplicate'))
        ms2_analysis_df = ms2_analysis_df.drop(columns=[col for col in ms2_analysis_df.columns if "_duplicate" in col])
    
    # if not ms2_analysis_df.empty and ms2_analysis_df.get('abundance'):
    if not ms2_analysis_df.empty and 'abundance' in ms2_analysis_df.columns:
        ms2_analysis_df['abundance_error'] = np.where(
            (ms2_analysis_df['i'].notnull()) & (ms2_analysis_df['abundance'].notnull()) & (ms2_analysis_df['abundance'] != 0),
            ((ms2_analysis_df['i'] - ms2_analysis_df['abundance']) / ms2_analysis_df['abundance']) * 100, None)
        abundance_condition = ms2_analysis_df['abundance_error'].apply(
            lambda x: np.abs(x) <= abundance_thresh_ms2 if x is not None else False)
        ms2_analysis_df['abundance_valid'] = np.where(abundance_condition, True, False)
    else:
         ms2_analysis_df['abundance_error'] = ms2_analysis_df['abundance_valid'] = None

    return ms2_analysis_df





"""Impute MS2 function"""

def impute_ms2(ms2_analysis_df, c1='query_name', c2='collision_type_energy', c3='filename', c4='query'):
    result_xdf = pd.DataFrame()
    if not ms2_analysis_df.empty:
        # Get unique values of c1 and c2
        unique_c1 = ms2_analysis_df[c1].unique()
        unique_c2 = ms2_analysis_df[c2].unique()
    
        # Map c1 to c4
        c1_to_c4 = dict(zip(ms2_analysis_df[c1], ms2_analysis_df[c4]))
        
        # Generate all combinations of c1, c2, and c3
        all_combinations = list(product(unique_c1, unique_c2, ms2_analysis_df[c3].unique()))
        
        # Create a DataFrame with all combinations
        result_xdf = pd.DataFrame(all_combinations, columns=[c1, c2, c3])
    
        # Map c4 values based on c1
        result_xdf[c4] = result_xdf[c1].map(c1_to_c4)
        
        # Merge with the original DataFrame to get the 'Value' column
        result_xdf = pd.merge(result_xdf, ms2_analysis_df, on=[c1, c2, c3, c4], how='left', indicator=True)
        result_xdf['imputed'] = result_xdf['_merge'].eq('left_only')
        result_xdf.drop('_merge', axis=1, inplace=True)
        
    return result_xdf





"""export_ms2_analysis_df"""

def export_ms2_analysis_df(ms2_analysis_df, data_directory, timestr):
    if not ms2_analysis_df.empty:
        ms2_analysis_df.to_csv(data_directory + "/MassQLab_Output/"+timestr+"/ms2_analysis_df.csv")
        sys.stdout.write(f"\nCreated ms2_analysis_df and exported as csv.")
        sys.stdout.flush()
    else:
        sys.stdout.write(f"\nms2_analysis_df empty. Not writing to csv.")
        sys.stdout.flush()





"""Cluster Plot MS2"""

def cluster_plot_ms2(ms2_analysis_df, data_directory, timestr):
    if not ms2_analysis_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots.pdf") as pdf:
            for group_name, grouped_df in ms2_analysis_df.groupby('collision_type'):
                for g_group_name, g_grouped_df in grouped_df.groupby(['query_name', 'query']):
                    if not grouped_df.empty and len(grouped_df) > 0:
                        pivot_df = g_grouped_df.pivot(index='energy', columns='filename', values='i')
                        pivot_df.plot(kind='bar', figsize=(10, 6))
                        plt.xlabel('Energy')
                        plt.ylabel('Intensity')
                        plt_title = f'{g_group_name[0]}_{group_name}'
                        plt_title_w_query = f'{g_group_name[0]}, {group_name}, {g_group_name[1]}'
                        wrapped_title = textwrap.fill(plt_title_w_query, width=100)
                        plt.title(wrapped_title, fontsize=8)
                        # plt.title(f'{g_group_name[1]}', fontsize=8)
                        pdf.savefig()
                        
                        if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots/"):
                            os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots/")
                        plt.savefig(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots/"+plt_title+'.png')
                        plt.close('all')
    sys.stdout.write(f"\nCreated ms2 cluster plot.")
    sys.stdout.flush()





"""Cluster Plot MS2 alt"""

def cluster_plot_ms2_alt(ms2_analysis_df, data_directory, timestr):
    if not ms2_analysis_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots_alt.pdf") as pdf:
            for group_name, grouped_df in ms2_analysis_df.groupby('collision_type'):
                for g_group_name, g_grouped_df in grouped_df.groupby(['query_name', 'query']):
                    if not grouped_df.empty and len(grouped_df) > 0:
                        pivot_df = g_grouped_df.pivot(index='filename', columns='energy', values='i')
                        pivot_df.plot(kind='bar', figsize=(10, 6))
                        plt.xlabel('Filename')
                        plt.ylabel('Intensity')
                        plt_title = f'{g_group_name[0]}_{group_name}'
                        plt_title_w_query = f'{g_group_name[0]}, {group_name}, {g_group_name[1]}'
                        wrapped_title = textwrap.fill(plt_title_w_query, width=100)
                        plt.title(wrapped_title, fontsize=8)
                        # plt.title(f'{g_group_name[1]}', fontsize=8)
                        pdf.savefig()
                        
                        if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots_alt/"):
                            os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots_alt/")
                        plt.savefig(data_directory + "/MassQLab_Output/"+timestr+"/ms2_cluster_plots_alt/"+plt_title+'.png')
                        plt.close('all')
    sys.stdout.write(f"\nCreated ms2 cluster plot alt.")
    sys.stdout.flush()





"""Cluster Plot MS2 group"""

def cluster_plot_ms2_group(ms2_analysis_df, data_directory, timestr):

    def sort_dataframe_column_unique(df, column_name):
        def find_lowest_number(s):
            numbers = re.findall(r'\d{2,}', s)
            numbers = [int(n) for n in numbers]
            return min(numbers) if numbers else float('inf')
        unique_values = df[column_name].unique()
        sorted_unique_values = sorted(unique_values, key=find_lowest_number)
        return sorted_unique_values

    if not ms2_analysis_df.empty:
        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr, "ms2_cluster_plots_group")
        output_dir_2 = os.path.join(data_directory, "MassQLab_Output", timestr)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        # Prepare a color map for filenames
        filenames = ms2_analysis_df['filename'].unique()
        colors = plt.cm.tab10(np.linspace(0, 1, len(filenames)))
        filename_to_color = {filename: color for filename, color in zip(filenames, colors)}
        
        # Define a list of marker styles
        markers = ['o', 's', '^', '>', '<', 'v', 'd', 'p', 'h', 'x']
        
        # Filter out NaN values from energies and ensure there are enough markers
        energies = ms2_analysis_df['energy'].dropna().unique()
        if len(energies) > len(markers):
            markers *= (len(energies) // len(markers)) + 1
        
        # Map non-NaN energies to markers
        energy_to_marker = {energy: markers[i] for i, energy in enumerate(sorted(energies))}

        with PdfPages(os.path.join(output_dir_2, "ms2_cluster_plots_group.pdf")) as pdf:
            for (group_name, collision_type), grouped_df in ms2_analysis_df.groupby(['group', 'collision_type']):
                fig, ax = plt.subplots(figsize=(12, 6))
                
                # Sort query names by length and alphabetically
                query_names = sort_dataframe_column_unique(grouped_df, 'query_name')
                
                energies = grouped_df['energy'].dropna().unique()
                
                # Closer spacing of energy levels within 'group' clusters
                energy_offset = np.linspace(-0.1, 0.1, len(energies))
                spacing = 0.25  # Adjusted spacing between groups for compactness
                query_positions = {name: i * spacing for i, name in enumerate(query_names)}

                for query_name in query_names:
                    query_group = grouped_df[grouped_df['query_name'] == query_name]
                    for filename, file_group in query_group.groupby('filename'):
                        color = filename_to_color[filename]
                        x_values = []
                        y_values = []
                        for energy_index, energy in enumerate(sorted(energies)):
                            marker = energy_to_marker.get(energy, None)  # Use get to handle potential NaNs safely
                            if marker:  # Only plot if energy is not NaN
                                energy_group = file_group[file_group['energy'] == energy]
                                x_vals = [query_positions[query_name] + energy_offset[energy_index]] * len(energy_group)
                                y_vals = energy_group['i'].values
                                ax.scatter(x_vals, y_vals, color=color, marker=marker)
                                x_values.extend(x_vals)
                                y_values.extend(y_vals)

                        # Connect points with lines
                        if len(x_values) > 0 and len(y_values) > 0:
                            ax.plot(x_values, y_values, color=color)

                ax.set_yscale('log')
                ax.set_xticks(list(query_positions.values()))
                ax.set_xticklabels(query_names, fontsize=6, rotation=45)
                plt.xlabel('Query Name', fontsize=8)
                plt.ylabel('Intensity', fontsize=8)
                plt.title(f'Group: {group_name}, Collision Type: {collision_type}', fontsize=10)

                # Adjust the legend for filenames
                custom_legend1 = [plt.Line2D([0], [0], color=color, marker='o', linestyle='', label=filename) for filename, color in filename_to_color.items()]
                legend1 = ax.legend(handles=custom_legend1, title='Filename', title_fontsize='8', loc='upper left', bbox_to_anchor=(1, 1), fontsize=6)

                # Create a second legend for marker shapes indicating energy levels, excluding NaN
                custom_legend2 = [plt.Line2D([0], [0], color='black', marker=marker, linestyle='', label=f'Energy {energy}') for energy, marker in energy_to_marker.items() if pd.notna(energy)]
                if custom_legend2:  # Only add legend if there are non-NaN energies
                    legend2 = ax.legend(handles=custom_legend2, title='Energy', title_fontsize='8', loc='upper left', bbox_to_anchor=(1, 0.5), fontsize=6)
                    ax.add_artist(legend1)  # Ensure the first legend is still displayed

                fig.tight_layout(rect=[0, 0, 0.75, 1])  # Adjust to make room for the legends
                pdf.savefig(fig, bbox_extra_artists=(legend1, legend2,) if custom_legend2 else (legend1,), bbox_inches='tight')
                plt.savefig(os.path.join(output_dir, f"{group_name}_{collision_type}.png"), bbox_extra_artists=(legend1, legend2,) if custom_legend2 else (legend1,), bbox_inches='tight')
                plt.close(fig)
    sys.stdout.write(f"\nCreated ms2 cluster grouped.")
    sys.stdout.flush()





"""Data Plot MS2"""

def summary_ms2(ms2_analysis_df, data_directory, timestr, abundance_thresh_ms2=20):
    if not ms2_analysis_df.empty:
        with PdfPages(data_directory + "/MassQLab_Output/"+timestr+"/ms2_summary_plots.pdf") as pdf:
            for i, (frame_group, frame_data) in enumerate(ms2_analysis_df.groupby(['query_name', 'collision_type_energy'])):
                colors_plot = plt.rcParams["axes.prop_cycle"]()
                plt.figure(figsize=(10, 6))
                if frame_group[1] == 'NA':
                    plt_title = frame_group[0]
                else:
                    plt_title = frame_group[0] + '_' + frame_group[1]
                plt.suptitle(plt_title, y=0.92)
                left_ax = plt.subplot(121)
                right_ax = plt.subplot(122)
                
                frame_data['filename_short'] = '..' + frame_data['filename'].str[-15:-5]
                frame_data['filename_short_Index'] = range(1, len(frame_data) + 1)
                for group_name, group_data in frame_data.groupby('filename'):
                    c = next(colors_plot)["color"]
                    group_data.plot(x='rt', y='i', kind='scatter', label=group_name, ax=left_ax, color=c, style='.-')
                    left_ax.set_ylabel('intensity')  # Add y-axis label
                    right_ax.bar(group_data['filename_short_Index'], group_data['i'])
        
                right_ax.set_xticks(frame_data['filename_short_Index'])
                right_ax.set_xticklabels(frame_data['filename_short'])
                right_ax.set_xlim(frame_data['filename_short_Index'].min() - 1, frame_data['filename_short_Index'].max() + 1)

                if 'abundance' in frame_data.columns:
                    frame_abundance = frame_data['abundance'].mean()
                    plt.tick_params(axis='x', rotation=90)
                    plt.axhline(y=(frame_abundance+(frame_abundance*(abundance_thresh_ms2/100))), color='black', linestyle='--')
                    plt.axhline(y=(frame_abundance-(frame_abundance*(abundance_thresh_ms2/100))), color='black', linestyle='--')
                
                if len(left_ax.get_lines()) >4:
                    left_ax.get_legend().remove()
        
                right_ax.tick_params(axis='x', rotation=90)
                left_ax.axhline(y=0, color='black', linestyle='--')
                right_ax.set_xlabel('')
                right_ax.set_ylabel('')
            
                plt.tight_layout(rect=[0, 0, 1, 0.96])
        
                pdf.savefig()
                if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms2_summary_plots/"):
                    os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms2_summary_plots/")
                plt.savefig(data_directory + "/MassQLab_Output/"+timestr+"/ms2_summary_plots/"+plt_title+'.png')
                # plt.show()
                plt.close('all')
        sys.stdout.write(f"\nCreated ms2_summary_plots.") 
        sys.stdout.flush()





"Save MS2 Scan"""

def save_ms2_scans_old(ms2_analysis_df, data_directory, timestr):
    if not ms2_analysis_df.empty:
        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr, "scans")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for filename_group_name, filename_group_df in ms2_analysis_df.groupby(["filename"]):
            reader = mzml.MzML(os.path.join(data_directory, filename_group_name[0]))
            filename_group_df = filename_group_df.dropna(subset=['scan'])
            pdf_filename = os.path.join(output_dir, f"{filename_group_name[0]}_scans.pdf")

            with PdfPages(pdf_filename) as pdf:
                for group_name, group_df in filename_group_df.groupby(["query_name", "collision_type_energy"]):
                    for index, row in group_df.iterrows():
                        scan_num = int(row['scan']) - 1
                        spectrum = reader[scan_num]
                        mz, intensity = spectrum['m/z array'], spectrum['intensity array']
                        height_min = np.max(intensity) / 20
                        peaks, properties = find_peaks(intensity, height=height_min, prominence=height_min, distance=20)

                        top_peaks_indices = np.argsort(properties["prominences"])[-5:]

                        plt.figure(figsize=(10, 6))
                        plt.plot(mz, intensity, label='Spectrum')
                        plt.scatter(mz[peaks[top_peaks_indices]], intensity[peaks[top_peaks_indices]], color='red', label='Top Peaks')

                        for peak_idx in top_peaks_indices:
                            plt.annotate(f'{mz[peaks[peak_idx]]:.2f}', 
                                         xy=(mz[peaks[peak_idx]], intensity[peaks[peak_idx]]), 
                                         xytext=(0, 5), 
                                         textcoords='offset points', 
                                         ha='center', fontsize=8)

                        plt.xlabel('m/z')
                        plt.ylabel('Intensity')
                        plot_title = f'{group_name[0]}, {group_name[1]}, scan:{scan_num}'
                        plt.title(plot_title, fontsize=8)
                        pdf.savefig()
                        plt.close('all')





"Save MS2 Scan"""

def save_ms2_scans(ms2_analysis_df, data_directory, timestr):
    if not ms2_analysis_df.empty:
        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr, "scans")
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        for filename_group_name, filename_group_df in ms2_analysis_df.groupby(["filename", "collision_type"]):
            reader = mzml.MzML(os.path.join(data_directory, filename_group_name[0]))
            filename_group_df = filename_group_df.dropna(subset=['scan'])
            pdf_filename = os.path.join(output_dir, f"{filename_group_name[0]}_{filename_group_name[1]}_scans.pdf")

            with PdfPages(pdf_filename) as pdf:
                for group_name, group_df in filename_group_df.groupby(["query_name", "energy"]):
                    for index, row in group_df.iterrows():
                        scan_num = int(row['scan']) - 1
                        spectrum = reader[scan_num]
                        mz, intensity = spectrum['m/z array'], spectrum['intensity array']
                        height_min = np.max(intensity) / 20
                        peaks, properties = find_peaks(intensity, height=height_min, prominence=height_min, distance=20)

                        top_peaks_indices = np.argsort(properties["prominences"])[-5:]

                        plt.figure(figsize=(10, 6))
                        plt.plot(mz, intensity, label='Spectrum')
                        plt.scatter(mz[peaks[top_peaks_indices]], intensity[peaks[top_peaks_indices]], color='red', label='Top Peaks')

                        for peak_idx in top_peaks_indices:
                            plt.annotate(f'{mz[peaks[peak_idx]]:.2f}', 
                                         xy=(mz[peaks[peak_idx]], intensity[peaks[peak_idx]]), 
                                         xytext=(0, 5), 
                                         textcoords='offset points', 
                                         ha='center', fontsize=8)

                        plt.xlabel('m/z')
                        plt.ylabel('Intensity')
                        plot_title = f'{group_name[0]}, Energy: {group_name[1]}, scan:{scan_num}'
                        plt.title(plot_title, fontsize=8)
                        pdf.savefig()
                        plt.close('all')





"""Initialize ReportLab"""

def on_page(canvas, doc):
    page_num = canvas.getPageNumber()
    canvas.drawCentredString(A4[0]/2, 50, str(page_num))

def on_page_landscape(canvas, doc):
  return on_page(canvas, doc)

def fig2image(f):
    buf = io.BytesIO()
    f.savefig(buf, format='png', dpi=300)
    buf.seek(0)
    x, y = f.get_size_inches()
    return Image(buf, x * inch, y * inch)

def df2table(df):
    return Table(
      [[Paragraph(col) for col in df.columns]] + df.values.tolist(), 
      style=[
        ('FONTNAME', (0,0), (-1,0), 'Helvetica-Bold'),
        ('FONTSIZE', (0, 0), (-1, -1), 8),
        ('LINEBELOW',(0,0), (-1,0), 1, colors.black),
        ('INNERGRID', (0,0), (-1,-1), 0.25, colors.black),
        ('BOX', (0,0), (-1,-1), 1, colors.black),
        ('ROWBACKGROUNDS', (0,0), (-1,-1), [colors.lightgrey, colors.white])],
      hAlign = 'LEFT')

def get_image(path, width):
    img = utils.ImageReader(path)
    iw, ih = img.getSize()
    aspect = ih / float(iw)
    return Image(path, width=width, height=(width * aspect))





"""MS1 build reportlab doc"""

def reportlab_ms1(ms1_analysis_df, data_directory, timestr):
    styles = getSampleStyleSheet()
    
    padding = dict(
      leftPadding=24, 
      rightPadding=24,
      topPadding=24,
      bottomPadding=24)
    
    portrait_frame = Frame(0, 0, *A4, **padding)
    landscape_frame = Frame(0, 0, *landscape(A4), **padding)
    
    portrait_template = PageTemplate(
      id='portrait', 
      frames=portrait_frame,
      onPage=on_page)
    
    landscape_template = PageTemplate(
      id='landscape', 
      frames=landscape_frame, 
      onPage=on_page_landscape)

    if not ms1_analysis_df.empty:
        for group_name, grouped_df in ms1_analysis_df.groupby(['query_name', 'query']):
        
            pdfname = data_directory + "/MassQLab_Output/"+timestr+"/ms1_reports/"+group_name[0]+"_ms1_report.pdf"
            doc = BaseDocTemplate(
              pdfname,
              pageTemplates=[
                portrait_template,
                landscape_template])
        
            image_path1 = data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_traces/rt_range/"+group_name[0]+".png"
            image_path2 = data_directory + "/MassQLab_Output/"+timestr+"/ms1_summary_areas/"+group_name[0]+".png"
        
            grouped_df1 = grouped_df[['filename', 'peak_area', 'fwhm', 'abundance_valid']]
            grouped_df2 = grouped_df[['filename', 'measured_RT', 'rt_error', 'peak_valid']]
        
            story = [
                Paragraph('MS1 Report: '+str(group_name[0]), styles['Heading1']),
                Paragraph(str(group_name[1]), styles['Heading2']),
                get_image(image_path1, width=7*inch),
                get_image(image_path2, width=5*inch),
                PageBreak(),
                df2table(grouped_df1),
                PageBreak(),
                df2table(grouped_df2)
            ]
        
            if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms1_reports/"):
                    os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms1_reports/")
                
            doc.build(story)
        sys.stdout.write(f"\nCreated ms1_reports.") 
        sys.stdout.flush()





"""MS2 build reportlab doc. Create"""

def reportlab_ms2(ms2_analysis_df, data_directory, timestr):
    
    styles = getSampleStyleSheet()
    
    padding = dict(
        leftPadding=24,
        rightPadding=24,
        topPadding=24,
        bottomPadding=24)
        
    portrait_frame = Frame(0, 0, *A4, **padding)
    landscape_frame = Frame(0, 0, *landscape(A4), **padding)
    
    portrait_template = PageTemplate(
      id='portrait', 
      frames=portrait_frame,
      onPage=on_page)
    
    landscape_template = PageTemplate(
      id='landscape', 
      frames=landscape_frame, 
      onPage=on_page_landscape)

    if not ms2_analysis_df.empty:
        for group_name, grouped_df in ms2_analysis_df.groupby(['query_name', 'query']):
        
            story = [Paragraph('MS2 Report: '+str(group_name[0]), styles['Heading1']),
                     Paragraph(str(group_name[1]), styles['Heading2'])
            ]
            
            pdfname = data_directory + "/MassQLab_Output/"+timestr+"/ms2_reports/"+group_name[0]+"_ms2_report.pdf"
            doc = BaseDocTemplate(
              pdfname,
              pageTemplates=[
                portrait_template,
                landscape_template])
        
            for group_name2, grouped_df2 in grouped_df.groupby('collision_type_energy'):    
                
                image_path1 = data_directory + "/MassQLab_Output/"+timestr+"/ms2_summary_plots/"+group_name[0]+"_"+group_name2+".png"
        
                grouped_data = grouped_df2[['filename', 'i', 'rt', 'abundance_error', 'abundance_valid']]
        
                story = story + (
                    [
                    # Paragraph(str(group_name[0]), styles['Heading2']),
                    Paragraph(str(group_name2), styles['Heading2']),
                    # Paragraph(str(group_name[1]), styles['Heading2']),
                    get_image(image_path1, width=7*inch),
                    df2table(grouped_data),
                    PageBreak(),])
                
            if not os.path.exists(data_directory + "/MassQLab_Output/"+timestr+"/ms2_reports/"):
                    os.makedirs(data_directory + "/MassQLab_Output/"+timestr+"/ms2_reports/")
                
            doc.build(story)
        sys.stdout.write(f"\nCreated ms2_reports.") 
        sys.stdout.flush()



""" system_test_summary """

def system_test_summary(data_directory, timestr, 
                        summary_traces="ms1_traces", 
                        summary_areas="ms1_summary_areas_system_test"):
    """
    Creates a consolidated grid image with rows aligned between subdir1 and subdir2.

    Parameters:
        data_directory (str): Base directory containing the output folders.
        timestr (str): Specific timestamp or folder name.
        summary_traces (str): Relative path to the summary traces directory.
        summary_areas (str): Relative path to the summary areas directory.

    Returns:
        str: Path to the saved PNG file.
    """
    def collect_images(subdir, filter_prefix=None):
        """
        Collect images from a directory.

        Parameters:
            subdir (str): Directory containing the images.

        Returns:
            list: List of RLImage objects.
        """
        images = []
        if os.path.exists(subdir):
            for filename in os.listdir(subdir):
                if filename.lower().endswith(('.png', '.jpg', '.jpeg')) and \
                   (filter_prefix is None or filename.startswith(filter_prefix)):
                    img_path = os.path.join(subdir, filename)
                    try:
                        img = RLImage.open(img_path)
                        images.append(img)
                    except Exception as e:
                        print(f"Error processing image {img_path}: {e}")
        return images

    # Collect images
    subdir1 = os.path.join(data_directory, "MassQLab_Output", timestr, summary_traces)
    subdir2 = os.path.join(data_directory, "MassQLab_Output", timestr, summary_areas)

    images1 = collect_images(subdir1, filter_prefix="system_test")
    images2 = collect_images(subdir2)

    # Ensure the same number of images in both directories
    max_images = max(len(images1), len(images2))
    placeholder_size = (100, 100) if not images1 else images1[0].size
    placeholder = RLImage.new("RGB", placeholder_size, (255, 255, 255))

    images1.extend([placeholder.copy() for _ in range(max_images - len(images1))])
    images2.extend([placeholder.copy() for _ in range(max_images - len(images2))])

    # Determine maximum dimensions for grid cells
    max_width = max(img.size[0] for img in images1 + images2)
    max_height = max(img.size[1] for img in images1 + images2)

    # Resize images while maintaining aspect ratio
    def resize_with_aspect(img, max_width, max_height):
        img_width, img_height = img.size
        scale = min(max_width / img_width, max_height / img_height)
        new_width = int(img_width * scale)
        new_height = int(img_height * scale)
        return img.resize((new_width, new_height))

    resized_images1 = [resize_with_aspect(img, max_width, max_height) for img in images1]
    resized_images2 = [resize_with_aspect(img, max_width, max_height) for img in images2]

    # Create a grid canvas
    grid_width = max_width * max_images
    grid_height = max_height * 2  # Two rows
    consolidated_img = RLImage.new("RGB", (grid_width, grid_height), (255, 255, 255))

    # Paste images into the grid with centering
    def paste_centered(base_img, img, x_offset, y_offset, cell_width, cell_height):
        x_center = x_offset + (cell_width - img.size[0]) // 2
        y_center = y_offset + (cell_height - img.size[1]) // 2
        base_img.paste(img, (x_center, y_center))

    # Top row
    for idx, img in enumerate(resized_images1):
        x_offset = idx * max_width
        paste_centered(consolidated_img, img, x_offset, 0, max_width, max_height)

    # Bottom row
    for idx, img in enumerate(resized_images2):
        x_offset = idx * max_width
        paste_centered(consolidated_img, img, x_offset, max_height, max_width, max_height)

    # Save the consolidated image
    output_png_path = os.path.join(data_directory, "MassQLab_Output", timestr, "aligned_grid.png")
    consolidated_img.save(output_png_path)
    sys.stdout.write(f"Consolidated image saved to {output_png_path}")
    sys.stdout.flush()



def main(data_directory=None, queryfile=None, metadata_file=None, metadata_filename_column=None, metadata_group_columns=None, kegg_path=None, convert_raw=None, msconvertexe=None, cache_setting=None, datasaver=None):
    if not data_directory and not queryfile:
        data_directory, queryfile, metadata_file, metadata_filename_column, metadata_group_columns, kegg_path, convert_raw, msconvertexe, cache_setting, datasaver = configure_MassQLab()
    if data_directory and queryfile:
        queries, ms1_query_df, ms2_query_df, query_groups, name_kegg_dict = create_queries(queryfile)
        if queries:
            convert_raw_files(convert_raw, msconvertexe, data_directory)
            file_count = mzml_file_count(data_directory)
            try:
                raw_df_ms1, raw_df_ms2, filename_groups, timestr = query_files(data_directory, queries, datasaver)
            except Exception as e:
                print(f"Exception caught: {e}")
                return
            if not raw_df_ms1.empty:
                ms1_analysis_df = analysis_ms1(raw_df_ms1, data_directory, timestr)
                ms1_analysis_df = ms1_query_analysis_merge(ms1_analysis_df, ms1_query_df)
                ms1_analysis_df = ms1_validity_and_QC(ms1_analysis_df)
                ms1_analysis_df = impute_ms1(ms1_analysis_df)
                rt_analysis_ms1(ms1_analysis_df, data_directory, timestr)
                summary_ms1_traces(raw_df_ms1, data_directory, timestr)
                summary_ms1_traces_inverse(raw_df_ms1, data_directory, timestr)
                if not ms1_analysis_df.empty:
                    export_ms1_analysis_df(ms1_analysis_df, data_directory, timestr)
                    summary_ms1_areas(ms1_analysis_df, data_directory, timestr)
                    summary_ms1_areas_inverse(ms1_analysis_df, data_directory, timestr)
                    reportlab_ms1(ms1_analysis_df, data_directory, timestr)
            if not raw_df_ms2.empty:
                plot_ms2(raw_df_ms2, data_directory, timestr)
                ms2_analysis_df = analysis_ms2(raw_df_ms2, ms2_query_df)
                ms2_analysis_df = impute_ms2(ms2_analysis_df)
                export_ms2_analysis_df(ms2_analysis_df, data_directory, timestr)
                cluster_plot_ms2(ms2_analysis_df, data_directory, timestr)      
                # cluster_plot_ms2_alt(ms2_analysis_df, data_directory, timestr)
                cluster_plot_ms2_group(ms2_analysis_df, data_directory, timestr)
                summary_ms2(ms2_analysis_df, data_directory, timestr)
                save_ms2_scans(ms2_analysis_df, data_directory, timestr)
                reportlab_ms2(ms2_analysis_df, data_directory, timestr)
        sys.stdout.write(f"\n--------------\n   Complete   \n--------------\n") 
        sys.stdout.flush()
        time.sleep(1)
        run()
    else:
        sys.stdout.write(f"\nNo data_directory and/or queryfile defined\n") 
        sys.stdout.flush()
        time.sleep(1)
        run()

def run():
    data_directory, queryfile, metadata_file, metadata_filename_column, metadata_group_columns, kegg_path, convert_raw, msconvertexe, cache_setting, datasaver = configure_MassQLab()
    if data_directory and queryfile:
        user_input = input("\n\n1) Enter 1 to run analysis.\n2) Enter 2 to reinitialize.\nAny other input closes.\nInput: ")
        if user_input == '1':
            main(data_directory, queryfile, metadata_file, metadata_filename_column, metadata_group_columns, kegg_path, convert_raw, msconvertexe, cache_setting, datasaver)
        elif user_input == '2':
            run()
        else:
            sys.stdout.write(f"\nExiting...\n") 
            sys.stdout.flush()
            time.sleep(3)
            pass
    else:
        user_input = input("\n\n2) Enter 2 to reinitialize.\nAny other input exits.\nInput: ")
        if user_input == '2':
            run()
        else:
            sys.stdout.write(f"\nExiting...\n") 
            sys.stdout.flush()
            time.sleep(3)
            pass