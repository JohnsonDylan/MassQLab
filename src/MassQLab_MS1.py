# MassQLab functions involved with downstram analysis of MS1 data in ms1_raw_df

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

def process_ms1_data(raw_df_ms1, ms1_query_df, data_directory, timestr,
                     fwhm_thresh=1, rt_thresh=0.1, datapoint_thresh=5,
                     abundance_thresh_ms1=10, impute=True, export=True):

    def index_to_xdata(xdata, indices):
        return interp1d(np.arange(len(xdata)), xdata)(indices)

    def gaussian(x, A, mu, sigma):
        return A * np.exp(-((x - mu) ** 2) / (2 * sigma ** 2))

    def generate_linear_array(rtmin, rtmax):
        return np.linspace(rtmin - (rtmax - rtmin), rtmax + (rtmax - rtmin), 1000)

    ms1_analysis_df_list = []
    plots = {}

    if not raw_df_ms1.empty:
        unique_filenames = raw_df_ms1['filename'].unique()
        unique_queries = raw_df_ms1['query'].unique()

        output_dir = os.path.join(data_directory, "MassQLab_Output", timestr, "ms1_traces")
        os.makedirs(output_dir, exist_ok=True)

        pdf_path = os.path.join(data_directory, "MassQLab_Output", timestr, "ms1_traces.pdf")
        with PdfPages(pdf_path) as pdf:
            for group_name, grouped_df in raw_df_ms1.groupby(['filename', 'query_name', 'query']):
                if len(grouped_df) < 5:
                    continue

                fig, ax = plt.subplots()
                analysis_df_dict = {
                    'filename': {0: group_name[0]},
                    'query_name': {0: group_name[1]},
                    'query': {0: group_name[2]},
                    'peak_area': {0: None},
                    'peak_area_alt': {0: None},
                    'fwhm': {0: None},
                    'rtmin': {0: None},
                    'rtmax': {0: None},
                    'measured_RT': {0: None},
                    'rt_error': {0: None},
                    'datapoint_count': {0: len(grouped_df)}
                }

                rtmin = extract_rtmin_value(group_name[2])
                rtmax = extract_rtmax_value(group_name[2])
                avg_rt = (rtmin + rtmax) / 2
                range_rt = rtmax - rtmin
                xdata_all = grouped_df.rt.values
                ydata_all = grouped_df.i.values

                def custom_score(peak):
                    normalized_height = ydata_all[peak] / max(ydata_all)
                    return normalized_height

                peaks, properties = find_peaks(
                    ydata_all,
                    height=abs(ydata_all.mean() * 1.1),
                    distance=5,
                    prominence=abs(np.max(ydata_all) / 10)
                )
                filtered_peaks = [p for p in peaks if (rtmin - range_rt) < xdata_all[p] < (rtmax + range_rt)]
                widths, height, l_ips, r_ips = peak_widths(ydata_all, filtered_peaks, rel_height=0.5)
                real_width = index_to_xdata(xdata_all, r_ips) - index_to_xdata(xdata_all, l_ips)
                final_peaks = [p for i, p in enumerate(filtered_peaks) if real_width[i] < range_rt]
                top_peaks = sorted(final_peaks, key=custom_score, reverse=True)[:2]

                if top_peaks:
                    peak = top_peaks[0]
                    peak_x, peak_y = xdata_all[peak], ydata_all[peak]
                    w, h, l, r = peak_widths(ydata_all, [peak], rel_height=0.5)
                    left, right = index_to_xdata(xdata_all, l)[0], index_to_xdata(xdata_all, r)[0]
                    fwhm = right - left
                    measured_rt = (left + right) / 2
                    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
                    interpolated_rt = generate_linear_array(rtmin, rtmax)
                    gaussian_y = gaussian(interpolated_rt, peak_y, measured_rt, sigma)
                    peak_area = peak_y * sigma * np.sqrt(2 * np.pi)
                    peak_area_alt = trapezoid(ydata_all, x=xdata_all)

                    ax.plot(interpolated_rt, gaussian_y, label='Gaussian Peak', color='r', lw=2)
                    ax.hlines(h, left, right, color='blue')

                    rt_error = 0 if rtmin <= measured_rt <= rtmax else measured_rt - (rtmin if measured_rt < rtmin else rtmax)
                    analysis_df_dict.update({
                        'peak_area': {0: peak_area},
                        'peak_area_alt': {0: peak_area_alt},
                        'fwhm': {0: fwhm},
                        'rtmin': {0: rtmin},
                        'rtmax': {0: rtmax},
                        'measured_RT': {0: measured_rt},
                        'rt_error': {0: rt_error}
                    })

                ax.plot(xdata_all, ydata_all, '-o', label='Data', color='black')
                for i, p in enumerate(top_peaks):
                    ax.plot(xdata_all[p], ydata_all[p], 'x' if i == 0 else '+', color='r')
                ax.axvline(rtmin, color='black', linestyle=':')
                ax.axvline(rtmax, color='black', linestyle=':')
                ax.set_xlim(rtmin - range_rt, rtmax + range_rt)
                ax.set_title(f"{group_name[0]}\n{group_name[1]}")
                ax.set_xlabel("rt")
                ax.set_ylabel("intensity")
                ax.legend()
                pdf.savefig(fig)

                plot_path = os.path.join(output_dir, f"{group_name[0]}_{group_name[1]}.png")
                fig.savefig(plot_path)
                plt.close(fig)

                plots[(group_name[0], group_name[2])] = plot_path
                ms1_analysis_df_list.append(pd.DataFrame(analysis_df_dict))

    ms1_analysis_df = pd.concat(ms1_analysis_df_list, ignore_index=True) if ms1_analysis_df_list else pd.DataFrame()

    # Merge with query dataframe
    if not ms1_analysis_df.empty and not ms1_query_df.empty:
        ms1_analysis_df = pd.merge(
            ms1_analysis_df,
            ms1_query_df.rename(columns={'name': 'query_name'}),
            on='query_name',
            how='inner',
            suffixes=('', '_duplicate')
        ).drop(columns=[col for col in ms1_analysis_df.columns if '_duplicate' in col])

    # QC Check
    if not ms1_analysis_df.empty:
        peak_condition = (
            (ms1_analysis_df['fwhm'] < fwhm_thresh) &
            (ms1_analysis_df['datapoint_count'] >= datapoint_thresh) &
            ms1_analysis_df.apply(lambda row: abs(row['rt_error']) <= row.get('rt_tolerance', rt_thresh)
                                  if row['rt_error'] is not None else False, axis=1)
        )
        ms1_analysis_df['peak_valid'] = peak_condition

        if 'abundance' in ms1_analysis_df.columns:
            ms1_analysis_df['abundance_error'] = np.where(
                (ms1_analysis_df['peak_area'].notnull()) &
                (ms1_analysis_df['abundance'].notnull()) &
                (ms1_analysis_df['abundance'] != 0),
                ((ms1_analysis_df['peak_area'] - ms1_analysis_df['abundance']) / ms1_analysis_df['abundance']) * 100,
                None
            )
            ms1_analysis_df['abundance_valid'] = ms1_analysis_df.apply(
                lambda row: abs(row['abundance_error']) <= row.get('abundance_tolerance', abundance_thresh_ms1)
                if row['abundance_error'] is not None else False,
                axis=1
            )
        else:
            ms1_analysis_df['abundance_error'] = ms1_analysis_df['abundance_valid'] = None

    # Impute
    if impute and not ms1_analysis_df.empty:
        c1, c3 = 'query_name', 'filename'
        all_combinations = list(product(ms1_analysis_df[c1].unique(), ms1_analysis_df[c3].unique()))
        result_xdf = pd.DataFrame(all_combinations, columns=[c1, c3])
        result_xdf = pd.merge(result_xdf, ms1_analysis_df, on=[c1, c3], how='left', indicator=True)
        result_xdf['imputed'] = result_xdf['_merge'] == 'left_only'
        result_xdf.drop(columns=['_merge'], inplace=True)
        result_xdf['peak_valid'] = result_xdf['peak_valid'].fillna(False)
        result_xdf['abundance_valid'] = result_xdf['abundance_valid'].fillna(False)
        ms1_analysis_df = result_xdf

    # Export CSV
    if export and not ms1_analysis_df.empty:
        output_csv = os.path.join(data_directory, "MassQLab_Output", timestr, "ms1_analysis_df.csv")
        ms1_analysis_df.to_csv(output_csv, index=False)
        sys.stdout.write("\nCreated ms1_analysis_df and exported as CSV.")
    else:
        sys.stdout.write("\nms1_analysis_df is empty or export skipped.")

    sys.stdout.flush()

    # === Create consolidated grid of plots ===
    if plots:
        fig, axes = plt.subplots(len(unique_filenames), len(unique_queries), figsize=(6 * len(unique_queries), 5 * len(unique_filenames)))
        
        if len(unique_filenames) == 1 and len(unique_queries) == 1:
            axes = np.array([[axes]])
        elif len(unique_filenames) == 1 or len(unique_queries) == 1:
            axes = axes.reshape((len(unique_filenames), len(unique_queries)))

        for i, filename in enumerate(unique_filenames):
            for j, query in enumerate(unique_queries):
                ax = axes[i][j]
                if (filename, query) in plots:
                    img = plt.imread(plots[(filename, query)])
                    ax.imshow(img)
                    ax.axis('off')
                else:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center', color='red')
                    ax.axis('off')

        plt.tight_layout()
        consolidated_path = os.path.join(data_directory, "MassQLab_Output", timestr, "consolidated_ms1_traces.png")
        plt.savefig(consolidated_path)
        plt.close()
        
    return ms1_analysis_df



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
                ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                fig1.tight_layout(rect=[0, 0, 0.85, 1])
                pdf.savefig(fig1, bbox_inches='tight')
                fig1.savefig(os.path.join(all_data_dir, f"{plt_title}_all_data.png"), bbox_inches='tight')
                plt.close(fig1)

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


