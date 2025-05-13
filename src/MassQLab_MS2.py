# MassQLab functions involved with downstram analysis of MS2 data in ms2_raw_df

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

"""Ingest raw ms2 dataframe and generate visualization"""
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



"""Ingest raw ms2 dataframe and identify top peak (by collision type and energy if applicable)"""
def process_ms2_data(raw_df_ms2, ms2_query_df, data_directory, timestr, abundance_thresh_ms2=20):
    ms2_analysis_df = pd.DataFrame()
    ms2_analysis_df_list = []

    # === MS2 Analysis ===
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

    if ms2_analysis_df_list:
        ms2_analysis_df = pd.concat(ms2_analysis_df_list).reset_index(drop=True)

    if not ms2_analysis_df.empty and not ms2_query_df.empty:
        ms2_analysis_df = pd.merge(
            ms2_analysis_df,
            ms2_query_df.rename(columns={'name': 'query_name'}),
            on='query_name',
            how='inner',
            suffixes=('', '_duplicate')
        )
        ms2_analysis_df.drop(columns=[col for col in ms2_analysis_df.columns if "_duplicate" in col], inplace=True)

    if not ms2_analysis_df.empty and 'abundance' in ms2_analysis_df.columns:
        ms2_analysis_df['abundance_error'] = np.where(
            (ms2_analysis_df['i'].notnull()) & (ms2_analysis_df['abundance'].notnull()) & (ms2_analysis_df['abundance'] != 0),
            ((ms2_analysis_df['i'] - ms2_analysis_df['abundance']) / ms2_analysis_df['abundance']) * 100,
            None
        )
        abundance_condition = ms2_analysis_df['abundance_error'].apply(
            lambda x: np.abs(x) <= abundance_thresh_ms2 if x is not None else False
        )
        ms2_analysis_df['abundance_valid'] = np.where(abundance_condition, True, False)
    else:
        ms2_analysis_df['abundance_error'] = ms2_analysis_df['abundance_valid'] = None

    # === Impute missing combinations ===
    if not ms2_analysis_df.empty:
        c1, c2, c3, c4 = 'query_name', 'collision_type_energy', 'filename', 'query'
        unique_c1 = ms2_analysis_df[c1].unique()
        unique_c2 = ms2_analysis_df[c2].unique()
        c1_to_c4 = dict(zip(ms2_analysis_df[c1], ms2_analysis_df[c4]))
        all_combinations = list(product(unique_c1, unique_c2, ms2_analysis_df[c3].unique()))
        result_xdf = pd.DataFrame(all_combinations, columns=[c1, c2, c3])
        result_xdf[c4] = result_xdf[c1].map(c1_to_c4)
        result_xdf = pd.merge(result_xdf, ms2_analysis_df, on=[c1, c2, c3, c4], how='left', indicator=True)
        result_xdf['imputed'] = result_xdf['_merge'].eq('left_only')
        result_xdf.drop('_merge', axis=1, inplace=True)
        ms2_analysis_df = result_xdf

    # === Export to CSV ===
    output_dir = os.path.join(data_directory, "MassQLab_Output", timestr)
    os.makedirs(output_dir, exist_ok=True)
    if not ms2_analysis_df.empty:
        ms2_analysis_df.to_csv(os.path.join(output_dir, "ms2_analysis_df.csv"), index=False)
        sys.stdout.write("\nCreated ms2_analysis_df and exported as csv.")
    else:
        sys.stdout.write("\nms2_analysis_df empty. Not writing to csv.")
    sys.stdout.flush()

    return ms2_analysis_df


"""visualize top ms2 peaks for each query clustered by filename, subclustered by energy"""
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


"""visualize top ms2 peaks for each query clustered by energy, subclustered by filename"""
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



"""visualize top ms2 peaks for all queries in query group (if present)"""
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



"""Create consolidated ms2 summary figures with traces and area barplot"""
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





"Load top ms2 scans and save as image"""
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
                        scan_num = int(row['scan'])
                        spectrum = reader[scan_num - 1]
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


"""Initialize ReportLab functions"""
def on_page(canvas, doc):
    page_num = canvas.getPageNumber()
    canvas.drawCentredString(A4[0]/2, 50, str(page_num))

def on_page_landscape(canvas, doc):
  return on_page(canvas, doc)

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



"""MS2 build reportlab doc"""
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

