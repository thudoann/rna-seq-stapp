import streamlit as st
import subprocess
import streamlit.components.v1 as components
import os
import pandas as pd

# Function to run a command and capture output
def run_command(command):
    result = subprocess.run(command, shell=True, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

# Function to display MultiQC report
def display_multiqc_report(multiqc_report_path):
    if os.path.exists(multiqc_report_path):
        with open(multiqc_report_path, 'r') as f:
            multiqc_report = f.read()
        components.html(multiqc_report, height=800, scrolling=True)
    else:
        st.error("MultiQC report not found.")

# Set up the Streamlit app
st.set_page_config(page_title="RNA-seq Pipeline Interface", layout="wide")

st.title("RNA-seq Pipeline Interface")

# Collect inputs in a sidebar
st.sidebar.header("Pipeline Configuration")

input_method = st.sidebar.radio(
    "Select input method:",
    ('Manual Input', 'Upload CSV')
)

if input_method == 'Manual Input':
    file_type = st.sidebar.radio(
        "Select your file type:",
        ('Single-end', 'Paired-end')
    )
    SRA_study = st.sidebar.text_input("Enter the Study number:")
    run_numbers = st.sidebar.text_area("Enter the run numbers (separated by spaces):")
    if SRA_study and run_numbers:
        run_numbers_list = run_numbers.split()
    else:
        run_numbers_list = []
else:
    uploaded_file = st.sidebar.file_uploader("Upload CSV file", type=["csv"])
    if uploaded_file:
        df = pd.read_csv(uploaded_file, delimiter =',')
        st.sidebar.write("File uploaded successfully:")
        st.sidebar.write(df.head())
        if not all(column in df.columns for column in ['SRA_study', 'Run', 'LibraryLayout']):
            st.sidebar.error("The CSV file must contain 'SRA_study', 'Run', and 'LibraryLayout' columns.")
            df = None
    else:
        df = None

kmer_size = st.sidebar.number_input("Enter the k-mer size:", min_value=1, value=15)

# Provide a button to run the pipeline
if st.sidebar.button("Run Pipeline"):
    if input_method == 'Manual Input':
        if not SRA_study or not run_numbers_list:
            st.sidebar.error("Please enter both the study number and run numbers.")
        elif kmer_size % 2 == 0:
            st.sidebar.error("The k-mer size must be an odd number.")
        else:
            st.write("## Pipeline Execution")
            with st.spinner("Running Python script..."):
                run_numbers_str = ' '.join(run_numbers_list)
                command = f"python fastqdl.py manual {file_type} {SRA_study} {run_numbers_str}"
                stdout, stderr, returncode = run_command(command)
                if returncode != 0:
                    st.error("Error running Python script.")
                    st.write("### Python Script Error Details")
                    st.code(stderr)
                else:
                    nf_script = f"{file_type.lower().replace('-end', '_end')}_rnaseq.nf"
                    output_dir = os.path.join("Output", SRA_study)
                    os.makedirs(output_dir, exist_ok=True)
                    with st.spinner(f"Running Nextflow script: {nf_script}..."):
                        command = f"nextflow run {nf_script} --kmer_size {kmer_size} --study_id {SRA_study} --output_dir {output_dir}"
                        st.write(f"Executing: {command}")  # Print the exact command being run for debugging
                        stdout, stderr, returncode = run_command(command)
                        if returncode != 0:
                            st.error("Error running Nextflow script.")
                            st.write("### Nextflow Script Error Details")
                            st.code(stderr)
                        else:
                            st.success("Nextflow script ran successfully!")
                            st.write("## MultiQC Report")
                            multiqc_report_path = os.path.join(output_dir, "multiqc_report.html")
                            display_multiqc_report(multiqc_report_path)
    elif input_method == 'Upload CSV':
        if df is None:
            st.sidebar.error("Please upload a CSV file.")
        elif kmer_size % 2 == 0:
            st.sidebar.error("The k-mer size must be an odd number.")
        else:
            st.write("## Pipeline Execution")
            with st.spinner("Processing CSV file..."):
                # Break down the CSV file by study number
                study_groups = df.groupby('SRA_study')
                for SRA_study, group in study_groups:
                    study_csv_path = os.path.join("uploaded_files", f"{SRA_study}.csv")
                    os.makedirs("uploaded_files", exist_ok=True)
                    group.to_csv(study_csv_path, index=False)

                    with st.spinner(f"Running Python script for study {SRA_study}..."):
                        command = f"python fastqdl.py csv {study_csv_path}"
                        stdout, stderr, returncode = run_command(command)
                        if returncode != 0:
                            st.error(f"Error running Python script for study {SRA_study}.")
                            st.write(f"### Python Script Error Details for study {SRA_study}")
                            st.code(stderr)
                        else:
                            st.success(f"Data download and preprocessing complete for study {SRA_study}. Running the Nextflow pipeline...")
                            file_type = group['LibraryLayout'].iloc[0].upper()
                            nf_script = f"{file_type.lower()}_end_rnaseq.nf"
                            output_dir = os.path.join("Output", SRA_study)
                            os.makedirs(output_dir, exist_ok=True)
                            with st.spinner(f"Running Nextflow script for study {SRA_study}: {nf_script}..."):
                                command = f"nextflow run {nf_script} --kmer_size {kmer_size} --study_id {SRA_study} --output_dir {output_dir}"
                                st.write(f"Executing: {command}")  # Print the exact command being run for debugging
                                stdout, stderr, returncode = run_command(command)
                                if returncode != 0:
                                    st.error(f"Error running Nextflow script for study {SRA_study}.")
                                    st.write(f"### Nextflow Script Error Details for study {SRA_study}")
                                    st.code(stderr)
                                else:
                                    st.success(f"Nextflow pipeline ran successfully for study {SRA_study}!")
                            st.write("## MultiQC Report")
                            multiqc_report_path = os.path.join(output_dir, "multiqc_report.html")
                            display_multiqc_report(multiqc_report_path)

# Instructions and Information
st.sidebar.header("Instructions")
st.sidebar.info("""
1. Select the input method: manual input or upload CSV.
2. If manual input, select the file type and enter the study number and run numbers.
3. If uploading a CSV, ensure it contains columns for study number, SRA numbers, and file type (single or paired).
4. Enter the k-mer size for indexing (must be odd).
5. Click on "Run Pipeline" to start the process.
""")

st.sidebar.header("About")
st.sidebar.info("""
This interface allows you to run an RNA-seq pipeline and visualize the results directly within the application.
""")

