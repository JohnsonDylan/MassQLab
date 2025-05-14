# high-level MassQLab workflow operations

import sys
import io, os, json, time, fnmatch, glob, warnings, subprocess, pandas as pd, numpy as np, regex as re, contextlib, textwrap
from pathlib import Path
from pyteomics import mzxml, mzml, auxiliary
from scipy.integrate import trapezoid
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

from MassQLab_core import *
from MassQLab_MS1 import *
from MassQLab_MS2 import *

"""high-level initialization workflow sequence"""
def initialize_config(data_directory, queryfile, output_directory, convert_raw, msconvertexe,
                      cache_setting, analysis):
    if not data_directory and not queryfile:
        return configure_MassQLab()
    return (data_directory, queryfile, output_directory, convert_raw, msconvertexe,
            cache_setting, analysis)

"""high-level core workflow sequence"""
def core_workflow(data_directory, queryfile, convert_raw, msconvertexe, output_directory):
    raw_df_ms1 = pd.DataFrame()
    raw_df_ms2 = pd.DataFrame()
    ms1_query_df = pd.DataFrame()
    ms2_query_df = pd.DataFrame()
    timestr = ""

    if data_directory and queryfile:
        sys.stdout.write(f"\nRun Start\n")
        queries, ms1_query_df, ms2_query_df, query_groups, name_kegg_dict = create_queries(queryfile)

        if queries:
            convert_raw_files(convert_raw, msconvertexe, data_directory)
            file_count = mzml_file_count(data_directory)

            try:
                raw_df_ms1, raw_df_ms2, filename_groups, timestr = query_files(data_directory, queries, output_directory)
            except Exception as e:
                print(f"Exception caught: {e}")
                return pd.DataFrame(), pd.DataFrame(), "", pd.DataFrame(), pd.DataFrame()
    else:
        sys.stdout.write(f"\nNo data_directory and/or queryfile defined\n")

    return raw_df_ms1, raw_df_ms2, timestr, ms1_query_df, ms2_query_df

"""high-level ms1 processing workflow sequence"""
def process_ms1(raw_df_ms1, ms1_query_df, data_directory, timestr, output_directory):
    ms1_analysis_df = process_ms1_data(raw_df_ms1, ms1_query_df, data_directory, timestr, output_directory)
    rt_analysis_ms1(ms1_analysis_df, data_directory, timestr, output_directory)
    summary_ms1_traces(raw_df_ms1, data_directory, timestr, output_directory)
    summary_ms1_traces_inverse(raw_df_ms1, data_directory, timestr, output_directory)

    if not ms1_analysis_df.empty:
        summary_ms1_areas(ms1_analysis_df, data_directory, timestr, output_directory)
        summary_ms1_areas_inverse(ms1_analysis_df, data_directory, timestr, output_directory)
        reportlab_ms1(ms1_analysis_df, data_directory, timestr, output_directory)

"""high-level ms2 processing workflow sequence"""
def process_ms2(raw_df_ms2, ms2_query_df, data_directory, timestr, output_directory):
    ms2_analysis_df = process_ms2_data(raw_df_ms2, ms2_query_df, data_directory, timestr, output_directory)
    plot_ms2(raw_df_ms2, data_directory, timestr, output_directory)

    if not ms2_analysis_df.empty:
        cluster_plot_ms2(ms2_analysis_df, data_directory, timestr, output_directory)
        cluster_plot_ms2_group(ms2_analysis_df, data_directory, timestr, output_directory)
        summary_ms2(ms2_analysis_df, data_directory, timestr, output_directory)
        save_ms2_scans(ms2_analysis_df, data_directory, timestr, output_directory)
        reportlab_ms2(ms2_analysis_df, data_directory, timestr, output_directory)

"""complete main workflow sequence"""
def main(data_directory=None, queryfile=None, output_directory=None, convert_raw=None, msconvertexe=None, cache_setting=None, analysis=None):
    data_directory, queryfile, output_directory, convert_raw, msconvertexe, cache_setting, analysis = initialize_config(
            data_directory, queryfile, output_directory, convert_raw, msconvertexe, 
            cache_setting, analysis)

    raw_df_ms1, raw_df_ms2, timestr, ms1_query_df, ms2_query_df = core_workflow(
        data_directory, queryfile, convert_raw, msconvertexe, output_directory)
    if analysis:
        if not raw_df_ms1.empty:
            process_ms1(raw_df_ms1, ms1_query_df, data_directory, timestr, output_directory)
    
        if not raw_df_ms2.empty:
            process_ms2(raw_df_ms2, ms2_query_df, data_directory, timestr, output_directory)

    sys.stdout.write(f"\nRun Complete\n")
    sys.stdout.flush()

def run():
    data_directory, queryfile, output_directory, convert_raw, msconvertexe, cache_setting, analysis = configure_MassQLab()
    if data_directory and queryfile:
        user_input = input("\n\n1) Enter 1 to run analysis.\n2) Enter 2 to reinitialize.\nAny other input closes.\nInput: ")
        if user_input == '1':
            main(data_directory, queryfile, output_directory, convert_raw, msconvertexe, cache_setting, analysis)
        elif user_input == '2':
            run()
        else:
            sys.stdout.write(f"\nExiting...\n") 
            sys.stdout.flush()
            time.sleep(1)
            pass
    else:
        sys.stdout.write(f"!!! data_directory and/or queryfile not found\n") 
        sys.stdout.flush()
        user_input = input("\n2) Enter 2 to reinitialize.\nAny other input exits.\nInput: ")
        if user_input == '2':
            run()
        else:
            sys.stdout.write(f"\nExiting...\n") 
            sys.stdout.flush()
            time.sleep(1)
            pass