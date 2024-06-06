import subprocess
import os
import shutil
import gzip
import sys
import pandas as pd

# Define number of CPUs to use
num_cpus = 4

def check_tool_installed(tool_name):
    """Check if a tool is installed and available in the PATH."""
    result = subprocess.run(['which', tool_name], capture_output=True, text=True)
    return result.returncode == 0

def download_sra(sra_id, sra_directory, transport_method="https"):
    print(f"Currently downloading: {sra_id} using {transport_method}")
    prefetch_command = ["prefetch", "--transport", transport_method, sra_id]
    print("The command used was:", prefetch_command)
    result = subprocess.run(prefetch_command, cwd=sra_directory, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"Error downloading {sra_id} with {transport_method}: {result.stderr}")
        if transport_method == "https":
            print("Trying with FTP transport method...")
            return download_sra(sra_id, sra_directory, transport_method="ftp")
        else:
            return False  # Skip to the next SRA ID
    return True

def extract_fastq(sra_id, sra_directory, fastq_directory):
    sra_id_folder = os.path.join(sra_directory, sra_id)
    sra_file_path = os.path.join(sra_id_folder, f"{sra_id}.sra")
    if os.path.exists(sra_file_path):
        print("Generating fastq for:", sra_id)
        fastq_dump_command = [
            "parallel-fastq-dump",
            "--outdir", fastq_directory,
            "--gzip",
            "--skip-technical",
            "--readids",
            "--read-filter", "pass",
            "--dumpbase",
            "--split-3",
            "--clip",
            "--threads", str(num_cpus),
            "--sra-id", sra_file_path
        ]
        print("The command used was:", fastq_dump_command)
        subprocess.run(fastq_dump_command)
    else:
        print("SRA file not found:", sra_file_path)

# Ensure required tools are installed
if not check_tool_installed("prefetch"):
    raise EnvironmentError("prefetch is not installed or not in PATH.")
if not check_tool_installed("parallel-fastq-dump"):
    raise EnvironmentError("parallel-fastq-dump is not installed or not in PATH.")

# Get input from command line arguments
input_type = sys.argv[1]

if input_type == "manual":
    file_type = sys.argv[2]
    study_number = sys.argv[3]
    sra_numbers = sys.argv[4:]
    studies = [(file_type, study_number, sra_numbers)]
elif input_type == "csv":
    csv_file = sys.argv[2]
    df = pd.read_csv(csv_file, delimiter=',')
    studies = []
    for study_number, group in df.groupby('SRA_study'):
        sra_numbers = group['Run'].tolist()
        file_type = group['LibraryLayout'].iloc[0]
        studies.append((file_type, study_number, sra_numbers))
else:
    raise ValueError("Invalid input type. Must be 'manual' or 'csv'.")

# Download and process each SRA file
for file_type, study_number, sra_numbers in studies:
    # Define directories for each study
    study_directory = os.path.join("Data", study_number)
    os.makedirs(study_directory, exist_ok=True)
    
    sra_directory = os.path.join(study_directory, "SRA")
    os.makedirs(sra_directory, exist_ok=True)
    
    fastq_directory = os.path.join(study_directory, "fastq")
    os.makedirs(fastq_directory, exist_ok=True)
    
    for sra_id in sra_numbers:
        if download_sra(sra_id, sra_directory):
            extract_fastq(sra_id, sra_directory, fastq_directory)

# Unzip .fastq.gz files
for _, study_number, _ in studies:
    fastq_directory = os.path.join("Data", study_number, "fastq")
    for fastq_file in os.listdir(fastq_directory):
        if fastq_file.endswith(".fastq.gz"):
            fastq_file_path = os.path.join(fastq_directory, fastq_file)
            unzipped_file_path = os.path.join(fastq_directory, fastq_file[:-3])  # Remove .gz extension
            with gzip.open(fastq_file_path, 'rb') as f_in:
                with open(unzipped_file_path, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            os.remove(fastq_file_path)  # Remove the original .fastq.gz file
