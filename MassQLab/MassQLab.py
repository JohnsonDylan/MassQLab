from MassQLab_functions import *
import tkinter as tk
from tkinter import filedialog, scrolledtext
import threading

# Redirects sys.stdout to a tkinter Text or scrolledtext widget
class RedirectText:
    def __init__(self, text_widget):
        self.text_widget = text_widget

    def write(self, msg):
        # Schedule the `insert_text` operation to run in the Tkinter main loop
        self.text_widget.after(0, self.insert_text, msg)

    def insert_text(self, msg):
        self.text_widget.configure(state='normal')  # Temporarily make the widget writable
        self.text_widget.insert(tk.END, msg)
        self.text_widget.see(tk.END)
        self.text_widget.configure(state='disabled')  # Disable editing again

    def flush(self):
        pass

# Main application class
class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("MassQLab")
        self.console = scrolledtext.ScrolledText(self, state='disabled', height=14)
        sys.stdout = RedirectText(self.console)
        sys.stderr = RedirectText(self.console)
        
        self.init_ui()

    def init_ui(self):
        # Your existing UI setup code...

        # Adjust the layout configuration
        self.grid_columnconfigure(1, weight=1)  # Makes middle column expandable

        # Retrieve configuration from configure_MassQLab
        configurations = configure_MassQLab()
        config_keys = ['data_directory', 'queryfile', 'metadata_file', 'metadata_filename_column', 'metadata_group_columns', 'kegg_path', 'convert_raw', 'msconvertexe', 'use_cache']
        self.grid_rowconfigure(len(config_keys)+1, weight=1)

        self.entries = {}
        self.check_vars = {}  # Dictionary to store BooleanVars for checkboxes
        checkbuttons = []  # List to store checkbuttons for later placement

        def log_change(var_name, index, mode):
            # Log changes for BooleanVars
            print(f"Configuration changed: {var_name} set to {self.check_vars[var_name].get()}")
        
        def on_entry_focus_out(event, key):
            # Log changes for Entries
            print(f"Configuration changed: {key} set to {event.widget.get()}")
    
        row_offset = 0  # Offset to track where to start placing checkbuttons

        for i, key in enumerate(config_keys):
            if key in ['convert_raw', 'use_cache']:
                self.check_vars[key] = tk.BooleanVar(value=configurations[i])
                self.check_vars[key].trace_add('write', lambda name, index, mode, var_name=key: log_change(var_name, index, mode))

                checkbutton = tk.Checkbutton(self, text=key, variable=self.check_vars[key], onvalue=True, offvalue=False)
                # Instead of placing now, store it with the intended row index
                checkbuttons.append((i, checkbutton))
                row_offset += 1  # Increment offset for each checkbutton found
            elif key in ['metadata_file', 'metadata_filename_column', 'metadata_group_columns', 'kegg_path']:
                pass
            else:
                tk.Label(self, text=key).grid(row=i - row_offset, column=0, sticky='w')
                entry = tk.Entry(self)
                entry.grid(row=i - row_offset, column=1, sticky='ew')

                entry.bind("<FocusOut>", lambda event, k=key: on_entry_focus_out(event, k))
                
                # Convert non-string values to strings before inserting
                value = configurations[i]
                if isinstance(value, list):
                    value = ', '.join(value)  # Convert list to comma-separated string
                elif isinstance(value, bool):
                    value = str(value)  # Convert boolean to string
                elif value is None:
                    value = "None"  # Convert None to a placeholder or empty string
                    
                entry.insert(0, value)
                self.entries[key] = entry
    
                if key in ['data_directory', 'queryfile', 'msconvertexe']:
                    btn_text = 'Browse Directory' if key == 'data_directory' else 'Browse File'
                    tk.Button(self, text=btn_text, command=lambda k=key: self.browse(k)).grid(row=i - row_offset, column=2)
        
        # Place checkbuttons at the bottom
        for i, (original_row, checkbutton) in enumerate(checkbuttons):
            checkbutton.grid(row=len(config_keys) - len(checkbuttons) + i, column=1, sticky='w')
            
        self.run_button = tk.Button(self, text="Run", command=self.run_main_tk, bg="lightgreen", fg="black", font=("Arial", 10, "bold"))
        # Move the "Run" button to the bottom right. Assuming 12 is a new row index after your last widget.

        # Adjust the console's placement to allow it to expand properly
        self.console.grid(row=len(config_keys)+1, column=0, columnspan=3, sticky='nsew', padx=5, pady=5)
        
        self.run_button.grid(row=len(config_keys), column=2, sticky='ew', padx=5, pady=5)

        # Make the window and its widgets resizable
        self.resizable(True, True)

        
    def browse(self, key):
        # Retrieve current value from the entry widget
        current_value = self.entries[key].get()
    
        # Use the user's home directory as a fallback initial directory
        home_dir = os.path.expanduser('~')
    
        # Determine the initial directory based on the current value
        if os.path.isdir(current_value):
            # If current_value is a directory, use it as initial_dir
            initial_dir = current_value
        elif os.path.isfile(current_value):
            # If current_value is a file, use its directory as initial_dir
            initial_dir = os.path.dirname(current_value)
        else:
            # Fallback to the home directory if current_value is not a valid path
            initial_dir = home_dir
    
        if key in ['data_directory']:
            selected_directory = os.path.normpath(filedialog.askdirectory(initialdir=initial_dir))
            if selected_directory and selected_directory != '.':
                self.entries[key].delete(0, tk.END)
                self.entries[key].insert(0, selected_directory)
                print(f"Configuration changed: {key} set to {selected_directory}")
        elif key in ['queryfile', 'msconvertexe']:  # For queryfile and msconvertexe, which are file paths
            selected_file = os.path.normpath(filedialog.askopenfilename(initialdir=initial_dir))
            if selected_file and selected_file != '.':
                self.entries[key].delete(0, tk.END)
                self.entries[key].insert(0, selected_file)
                print(f"Configuration changed: {key} set to {selected_file}")

    def run_main_tk(self):
        self.run_button.config(state='disabled', bg='red', fg='white')
        args = [self.entries[key].get() for key in self.entries]
        threading.Thread(target=lambda: self.execute_main_tk(*args)).start()

    def execute_main_tk(self, *args):
        main_tk(*args)
        self.run_button.config(state='normal', bg='lightgreen', fg='black')

    def on_closing(self):
        sys.stdout = sys.__stdout__
        self.destroy()

# Ensure you replace this with the actual main_tk function code you provided earlier, with appropriate modifications if necessary.
def main_tk(data_directory=None, queryfile=None, metadata_file=None, metadata_filename_column=None, metadata_group_columns=None, kegg_path=None, convert_raw=None, msconvertexe=None, cache_setting=None, datasaver=None):
    if data_directory and queryfile:
        sys.stdout.write(f"\nRun Start\n") 
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
                    summary_ms1_areas_system_test(ms1_analysis_df, data_directory, timestr)
                    summary_ms1_areas_inverse(ms1_analysis_df, data_directory, timestr)
                    reportlab_ms1(ms1_analysis_df, data_directory, timestr)
            if not raw_df_ms2.empty:
                plot_ms2(raw_df_ms2, data_directory, timestr)
                ms2_analysis_df = analysis_ms2(raw_df_ms2, ms2_query_df)
                ms2_analysis_df = impute_ms2(ms2_analysis_df)
                export_ms2_analysis_df(ms2_analysis_df, data_directory, timestr)
                cluster_plot_ms2(ms2_analysis_df, data_directory, timestr)      
                cluster_plot_ms2_alt(ms2_analysis_df, data_directory, timestr)
                cluster_plot_ms2_group(ms2_analysis_df, data_directory, timestr)
                summary_ms2(ms2_analysis_df, data_directory, timestr)
                save_ms2_scans(ms2_analysis_df, data_directory, timestr)
                reportlab_ms2(ms2_analysis_df, data_directory, timestr)
        sys.stdout.write(f"\nRun Complete\n") 
    else:
        sys.stdout.write(f"\nNo data_directory and/or queryfile defined\n") 

if __name__ == "__main__":
    app = App()
    app.protocol("WM_DELETE_WINDOW", app.on_closing)
    app.mainloop()
