#!/bin/bash

# =============================================================================
# Script Name: fastq_processing_pipeline.sh
# Description: Automates the preprocessing and quality control of paired-end
#              FASTQ data using FastQC, MultiQC, and fastp. It creates the 
#              necessary directory structure, allows environment selection,
#              logs read statistics, and generates a read statistics graph.
# Author: Daniel João
# Date: 09-12-2024
# =============================================================================

# =============================================================================
# Exit immediately if a command exits with a non-zero status.
# This ensures that the script stops executing upon encountering an error.
# =============================================================================
set -e

# ---------------------------
# Bash Version Check
# ---------------------------
# Ensures the script is run with Bash version 4 or higher for associative array support.
if (( BASH_VERSINFO[0] < 4 )); then
    echo "❌  Error: This script requires Bash version 4 or higher."
    echo "    Current version: $BASH_VERSION"
    exit 1
fi

# ---------------------------
# Step Title Definitions
# ---------------------------
# These variables define the titles for each major step in the pipeline.
# They are used to display clear and consistent headers in the script's output.
TITLE_PIPELINE="FASTQ Processing Pipeline"
TITLE_CONDA="Conda Environment Selection"
TITLE_CHECK_TOOLS="Checking Required Tools"
TITLE_CREATE_DIR="Creating Directory Structure"
TITLE_MOVE_FILES="Organizing FASTQ Files"
TITLE_IDENTIFY_PAIRS="Identifying Paired-End Files"
TITLE_RUN_FASTQC_RAW="Running FastQC on Raw Data"
TITLE_RUN_MULTIQC_RAW="Running MultiQC on Raw Data"
TITLE_CONFIG_FASTP="Configuring Fastp Parameters"
TITLE_RUN_FASTP="Running Fastp on Raw Data"
TITLE_RUN_FASTQC_PROCESSED="Running FastQC on Processed Data"
TITLE_RUN_MULTIQC_PROCESSED="Running MultiQC on Processed Data"
TITLE_EXTRACT_STAT_PROCESSED="Extracting Fastp Read Statistics and Generating Graph"
TITLE_COMPLETED="Pipeline Completed"

# ---------------------------
# Function: display_step
# Description: Displays a formatted title for each step in the pipeline.
# Parameters:
#   $1 - The title of the step to display.
# ---------------------------
display_step() {
    echo ""
    echo "===== $1 ====="
    echo ""
}

# ---------------------------
# Function: check_tools
# Description: Verifies that all required tools are installed and accessible.
# Parameters:
#   $@ - A list of tool names to check.
# ---------------------------
check_tools() {
    echo ""
    for tool in "$@"; do
        if command -v "$tool" &> /dev/null; then
            echo "✔️  Tool '$tool' is installed and accessible."
        else
            echo "❌  Error: Tool '$tool' is not installed or not in PATH."
            exit 1
        fi
    done
    echo ""
    echo "All required tools are present."
    echo ""
}

# ---------------------------
# Function: source_conda
# Description: Sources the Conda initialization script to enable environment activation.
#              It checks common installation paths for Conda.
# ---------------------------
source_conda() {
    echo ""
    # Attempt to source conda.sh from common installation paths
    if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/miniconda3/etc/profile.d/conda.sh"
        echo "✔️  Sourced Conda from Miniconda3."
    elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
        source "$HOME/anaconda3/etc/profile.d/conda.sh"
        echo "✔️  Sourced Conda from Anaconda3."
    elif [ -f "/opt/conda/etc/profile.d/conda.sh" ]; then
        source "/opt/conda/etc/profile.d/conda.sh"
        echo "✔️  Sourced Conda from /opt/conda."
    else
        echo "❌  Error: conda.sh not found. Please ensure Conda is installed."
        exit 1
    fi
    echo ""
}

# ---------------------------
# Function: ask_yes_no
# Description: Prompts the user with a yes/no question and validates the input.
# Parameters:
#   $1 - The prompt message to display.
# Returns:
#   0 if the user answers yes, 1 if no.
# ---------------------------
ask_yes_no() {
    while true; do
        read -rp "$1 (y/n): " yn
        case "$yn" in
            [Yy]*) return 0 ;;
            [Nn]*) return 1 ;;
            *) echo "❌  Please answer yes (y) or no (n)." >&2 ;;
        esac
    done
}

# ---------------------------
# Function: ask_number
# Description: Prompts the user to enter a number within a specified range and validates the input.
# Parameters:
#   $1 - The prompt message to display.
#   $2 - The default value if the user provides no input.
#   $3 - The minimum acceptable value (inclusive).
#   $4 - The maximum acceptable value (optional).
# Returns:
#   The validated number entered by the user.
# ---------------------------
ask_number() {
    local prompt="$1"
    local default="$2"
    local min="$3"
    local max="$4"
    local input
    while true; do
        read -rp "$prompt " input
        # Use default if input is empty
        input="${input:-$default}"
        # Check if input is a valid number (only digits)
        if [[ "$input" =~ ^[0-9]+$ ]]; then
            if [[ -n "$max" ]]; then
                if (( input >= min && input <= max )); then
                    echo "$input"
                    return 0
                fi
            else
                if (( input >= min )); then
                    echo "$input"
                    return 0
                fi
            fi
        fi
        # Redirect error messages to stderr to prevent them from being captured
        echo "❌  Invalid input. Please enter a number" >&2
        if [[ -n "$max" ]]; then
            echo "    between $min and $max." >&2
        else
            echo "    greater than or equal to $min." >&2
        fi
    done
}

# ---------------------------
# Initialize Associative Array Globally
# ---------------------------
declare -gA fastp_options=()

# ---------------------------
# Function: configure_fastp
# Description: Prompts the user to configure fastp parameters with input validation.
#              Includes Quality Filtering, Length Filtering, Low Complexity Filtering,
#              Adapter Trimming, PolyG Tail Trimming, PolyX Tail Trimming,
#              and Deduplication.
# ---------------------------
configure_fastp() {
    echo ""
    
    # Quality Filtering
    if ask_yes_no "Enable Quality Filtering?"; then
        qual_phred=$(ask_number "Set Qualified Quality PHRED (0-40, default=15):" 15 0 40)
        fastp_options["--qualified_quality_phred"]="${qual_phred}"
        
        unqual_percent=$(ask_number "Set Unqualified Percent Limit (0-100, default=40):" 40 0 100)
        fastp_options["--unqualified_percent_limit"]="${unqual_percent}"
    else
        fastp_options["--disable_quality_filtering"]="true"
    fi
    echo ""
    
    # Length Filtering
    if ask_yes_no "Enable Length Filtering?"; then
        min_len=$(ask_number "Set Minimum Length Required (default=15):" 15 1)
        fastp_options["--length_required"]="${min_len}"
        
        max_len=$(ask_number "Set Maximum Length Limit (0 = no limit, default=0):" 0 0)
        fastp_options["--length_limit"]="${max_len}"
    else
        fastp_options["--disable_length_filtering"]="true"
    fi
    echo ""
    
    # Low Complexity Filtering
    if ask_yes_no "Enable Low Complexity Filtering?"; then
        fastp_options["--low_complexity_filter"]="true"
        complexity_thresh=$(ask_number "Set Complexity Threshold (0-100, default=30):" 30 0 100)
        fastp_options["--complexity_threshold"]="${complexity_thresh}"
    else
        fastp_options["--low_complexity_filter"]="false"
    fi
    echo ""
    
    # Adapter Trimming
    if ask_yes_no "Enable Adapter Trimming?"; then
        # Adapter trimming is enabled by default; no additional options needed
        :
    else
        fastp_options["--disable_adapter_trimming"]="true"
    fi
    echo ""
    
    # PolyG Tail Trimming
    if ask_yes_no "Enable PolyG Tail Trimming?"; then
        poly_g_min_len=$(ask_number "Set PolyG minimum length (default=10):" 10 1)
        fastp_options["--trim_poly_g"]="${poly_g_min_len}"
    else
        fastp_options["--disable_trim_poly_g"]="true"
    fi
    echo ""
    
    # PolyX Tail Trimming
    if ask_yes_no "Enable PolyX Tail Trimming?"; then
        poly_x_min_len=$(ask_number "Set PolyX minimum length (default=10):" 10 1)
        fastp_options["--trim_poly_x"]="${poly_x_min_len}"
    else
        # PolyX tail trimming is disabled by default; no action needed
        :
    fi
    echo ""
    
    # Deduplication
    if ask_yes_no "Enable Deduplication?"; then
        fastp_options["--dedup"]="true"
    fi
    echo ""
    
    # Debug: Display configured fastp options
    echo "🔧 Configured Fastp Options:"
    for key in "${!fastp_options[@]}"; do
        echo "$key: ${fastp_options[$key]}"
    done
    echo ""
    
    echo "✅ Fastp configuration complete."
    echo ""
}

# ---------------------------
# Function: plot_read_statistics
# Description: Generates a multiple bar plot of read statistics using python3 plotext.
# ---------------------------
plot_read_statistics() {
    echo ""
    
    # Check if plotext is installed
    if ! python3 -c "import plotext" &> /dev/null; then
        echo "❌  plotext is not installed. Installing plotext..."
        pip install plotext
        echo "✔️  plotext installed successfully."
    fi
    
    # Generate the plot using Python and plotext
    python3 - <<END
import plotext as plt
import sys

# Read the readLOG.txt file
samples = []
before_reads = []
after_reads = []
discarded_reads = []

try:
    with open('readLOG.txt', 'r') as f:
        next(f)  # Skip header
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) != 4:
                continue
            sample, before, after, discarded = parts
            samples.append(sample)
            before_reads.append(int(before))
            after_reads.append(int(after))
            discarded_reads.append(int(discarded))
except FileNotFoundError:
    print("❌  readLOG.txt not found. Cannot generate graph.")
    sys.exit(1)
except ValueError:
    print("❌  Invalid data format in readLOG.txt.")
    sys.exit(1)

# Check if there is data to plot
if not samples:
    print("❌  No valid data found in readLOG.txt. Cannot generate graph.")
    sys.exit(1)

# Plot Before, After, and Discarded Reads as Multiple Bar Plot
plt.simple_multiple_bar(
    samples,
    [before_reads, after_reads, discarded_reads],
    width=100,
    labels=["Before Filtering", "After Filtering", "Discarded Reads"],
    title="Fastp Read Statistics"
)
plt.show()
END

    echo ""
    echo "✅ Read Statistics Graph generated successfully."
    echo ""
}

# ---------------------------
# ---------------------------
# START
# ---------------------------
# Initialize Conda and Display Pipeline Title
# ---------------------------
display_step "$TITLE_PIPELINE"
source_conda

# ---------------------------
# Conda Environment Selection
# ---------------------------
display_step "$TITLE_CONDA"

echo "Available Conda Environments:"
# Retrieve a list of Conda environments, excluding commented lines and the 'base' environment
mapfile -t envs < <(conda env list | awk '{print $1}' | grep -v "^#\|^\$" | grep -v "^base$")

# Check if any environments are available
if [ ${#envs[@]} -eq 0 ]; then
    echo "❌  No Conda environments found excluding 'base'. Please create one or adjust the script."
    exit 1
fi

# Display the list of available environments with numbering
for i in "${!envs[@]}"; do
    printf "%d) %s\n" $((i + 1)) "${envs[$i]}"
done

# Prompt the user to select a Conda environment by number
while true; do
    read -rp "Select a Conda environment by number: " env_num
    # Validate the input to ensure it's a number within the valid range
    if [[ "$env_num" =~ ^[0-9]+$ ]] && (( env_num >= 1 && env_num <= ${#envs[@]} )); then
        selected_env="${envs[$((env_num - 1))]}"
        break
    else
        echo "❌  Invalid selection. Please enter a number between 1 and ${#envs[@]}." >&2
    fi
done

# Activate the selected Conda environment
echo "🔄 Activating Conda environment: ${selected_env}"
conda activate "${selected_env}"
echo ""

# ---------------------------
# Check for Required Tools
# ---------------------------
display_step "$TITLE_CHECK_TOOLS"
# List of required tools: fastqc, multiqc, fastp, jq
check_tools fastqc multiqc fastp jq

# ---------------------------
# Create Directory Structure
# ---------------------------
display_step "$TITLE_CREATE_DIR"

# Define the names of the raw and processed data directories
RAW_DIR="rawdata"
PROCESSED_DIR="processed_data"

# Ask the user if directories are already created
if ask_yes_no "Are the necessary directories already created?"; then
    echo ""
    echo "✅ Assuming directories already exist."
else
    # Create the necessary directory structure for reports, including a dedicated fastp/ folder
    echo ""
    echo "📁 Creating directory structure..."
    mkdir -p "${RAW_DIR}/reports/fastqc"
    mkdir -p "${RAW_DIR}/reports/multiqc"
    mkdir -p "${PROCESSED_DIR}/reports/fastqc"
    mkdir -p "${PROCESSED_DIR}/reports/multiqc"
    mkdir -p "${PROCESSED_DIR}/reports/fastp"  # New directory for fastp reports
    
    # Inform the user about the created directories
    echo "📁 Directory structure created:"
    echo "- ${RAW_DIR}/"
    echo "  - reports/"
    echo "    - fastqc/"
    echo "    - multiqc/"
    echo "- ${PROCESSED_DIR}/"
    echo "  - reports/"
    echo "    - fastqc/"
    echo "    - multiqc/"
    echo "    - fastp/"
fi
echo ""

# ---------------------------
# Move FASTQ Files to Raw Data Directory
# ---------------------------
display_step "$TITLE_MOVE_FILES"

# Check if any FASTQ files exist in the current directory
if ls *.fastq.gz 1> /dev/null 2>&1; then
    # Move the files to the rawdata directory
    echo "📦 Moving FASTQ files to ${RAW_DIR}/..."
    mv *.fastq.gz "${RAW_DIR}/"
    echo "✅ FASTQ files moved to ${RAW_DIR}/"
else
    echo "⚠️  No FASTQ files found to move."
fi
echo ""

# ---------------------------
# Identify Paired-End Files
# ---------------------------
display_step "$TITLE_IDENTIFY_PAIRS"

# Change directory to rawdata
cd "${RAW_DIR}" || exit

# Find all R1 (forward) FASTQ files based on the naming pattern
read1_files=(*_plus_1_aaa.fastq.gz)

# Check if any paired-end files exist
if [ ${#read1_files[@]} -eq 0 ]; then
    echo "❌ No paired-end FASTQ files found with pattern *_plus_1_aaa.fastq.gz"
    exit 1
fi

# Initialize an array to store sample names by stripping the R1 suffix
declare -a samples
for file in "${read1_files[@]}"; do
    # Extract the sample name by removing the '_plus_1_aaa.fastq.gz' suffix
    sample="${file/_plus_1_aaa.fastq.gz/}"
    samples+=("$sample")
done

# Inform the user about the number of paired-end samples found
echo "🔍 Found ${#samples[@]} paired-end samples."
echo ""

# ---------------------------
# Run FastQC on Raw Data
# ---------------------------
display_step "$TITLE_RUN_FASTQC_RAW"

# Iterate over each sample and run FastQC on both R1 and R2 FASTQ files
for sample in "${samples[@]}"; do
    read1="${sample}_plus_1_aaa.fastq.gz"
    read2="${sample}_plus_2_aaa.fastq.gz"

    # Check if both R1 and R2 files exist
    if [[ -f "$read1" && -f "$read2" ]]; then
        echo "🔧 Running FastQC for sample: $sample"
        # Run FastQC and output the results to the designated fastqc directory
        fastqc -o reports/fastqc/ "$read1" "$read2" > /dev/null 2>&1
        echo "✔️  FastQC completed for sample: $sample"
    else
        echo "⚠️  Warning: Pair not found for sample $sample. Skipping FastQC."
    fi
done
echo ""

# ---------------------------
# Run MultiQC on Raw Data
# ---------------------------
display_step "$TITLE_RUN_MULTIQC_RAW"

echo "📊 Running MultiQC on FastQC reports..."
# Aggregate FastQC reports into a single MultiQC report
multiqc reports/fastqc/ -o reports/multiqc/ > /dev/null 2>&1
echo "✔️  MultiQC report generated in ${RAW_DIR}/reports/multiqc/"
echo ""

# ---------------------------
# Configure Fastp Parameters
# ---------------------------
display_step "$TITLE_CONFIG_FASTP"

# Call the function to configure fastp based on user input
configure_fastp

# ---------------------------
# Run Fastp on Raw Data
# ---------------------------
display_step "$TITLE_RUN_FASTP"

# Iterate over each sample and run fastp for data trimming and filtering
for sample in "${samples[@]}"; do
    read1="${sample}_plus_1_aaa.fastq.gz"
    read2="${sample}_plus_2_aaa.fastq.gz"

    # Check if both R1 and R2 files exist before processing
    if [[ -f "$read1" && -f "$read2" ]]; then
        echo "✂️  Running fastp for sample: $sample"
        # Define output filenames and report paths within the dedicated fastp/ directory
        out1="../${PROCESSED_DIR}/${sample}_plus_1_aaa.fastq.gz"
        out2="../${PROCESSED_DIR}/${sample}_plus_2_aaa.fastq.gz"
        report="../${PROCESSED_DIR}/reports/fastp/${sample}_fastp.json"
        html_report="../${PROCESSED_DIR}/reports/fastp/${sample}_fastp.html"
        log_file="../${PROCESSED_DIR}/reports/fastp/${sample}_fastp.log"

        # Build the fastp command using a Bash array for safer argument handling
        fastp_cmd=(fastp -i "$read1" -I "$read2" -o "$out1" -O "$out2" --json "$report" --html "$html_report")

        # Append user-selected fastp options to the command
        for option in "${!fastp_options[@]}"; do
            value=${fastp_options[$option]}
            if [[ "$value" == "true" ]]; then
                fastp_cmd+=("$option")
            elif [[ "$value" != "false" && -n "$value" ]]; then
                fastp_cmd+=("$option" "$value")
            fi
        done

        # Debug: Display the constructed fastp command
        echo "📋 Running command: ${fastp_cmd[@]}"

        # Execute the fastp command and redirect output to a log file
        if "${fastp_cmd[@]}" > "$log_file" 2>&1; then
            echo "✔️  Fastp completed for sample: $sample"
        else
            echo "❌  Error: fastp failed for sample: $sample. Check the log file at $log_file"
            exit 1
        fi
    else
        echo "⚠️  Warning: Pair not found for sample $sample. Skipping fastp."
    fi
done
echo ""

# ---------------------------
# Run FastQC on Processed Data
# ---------------------------
display_step "$TITLE_RUN_FASTQC_PROCESSED"

echo "🔧 Running FastQC for processed FASTQ files..."
# Navigate back to the parent directory to access processed_data
cd ../

# Iterate over each sample and run FastQC on the processed (trimmed) FASTQ files
for sample in "${samples[@]}"; do
    read1="${PROCESSED_DIR}/${sample}_plus_1_aaa.fastq.gz"
    read2="${PROCESSED_DIR}/${sample}_plus_2_aaa.fastq.gz"

    # Check if both processed R1 and R2 files exist
    if [[ -f "$read1" && -f "$read2" ]]; then
        echo "🔧 Running FastQC for processed sample: $sample"
        # Run FastQC and output the results to the designated fastqc directory
        fastqc -o "${PROCESSED_DIR}/reports/fastqc/" "$read1" "$read2" > /dev/null 2>&1
        echo "✔️  FastQC completed for processed sample: $sample"
    else
        echo "⚠️  Warning: Processed pair not found for sample $sample. Skipping FastQC."
    fi
done
echo ""

# ---------------------------
# Run MultiQC on Processed Data
# ---------------------------
display_step "$TITLE_RUN_MULTIQC_PROCESSED"

echo "📊 Running MultiQC on processed FastQC reports..."
# Change directory to processed_data to run MultiQC
cd "${PROCESSED_DIR}" || exit
# Aggregate FastQC reports into a single MultiQC report
multiqc reports/fastqc/ -o reports/multiqc/ > /dev/null 2>&1
echo "✔️  MultiQC report generated in ${PROCESSED_DIR}/reports/multiqc/"
echo ""

# ---------------------------
# Extract and Display Fastp Read Statistics
# ---------------------------
display_step "$TITLE_EXTRACT_STAT_PROCESSED"

echo "📈 Extracting Fastp read statistics..."
# Change directory to the fastp report directory
cd reports/fastp/ || exit

# Initialize the read statistics log file with headers in the current directory
echo -e "Sample\tBefore_Reads\tAfter_Reads\tDiscarded_Reads" > readLOG.txt

# Iterate over each sample to extract read statistics from fastp JSON reports
for sample in "${samples[@]}"; do
    json_report="${sample}_fastp.json"

    # Check if the fastp JSON report exists
    if [[ -f "$json_report" ]]; then
        # Extract read counts using jq
        before_reads=$(jq '.summary.before_filtering.total_reads' "$json_report")
        after_reads=$(jq '.summary.after_filtering.total_reads' "$json_report")
        # Calculate discarded reads
        discarded_reads=$((before_reads - after_reads))

        # Append the statistics to the readLOG.txt file
        echo -e "${sample}\t${before_reads}\t${after_reads}\t${discarded_reads}" >> readLOG.txt
    else
        echo "⚠️  Warning: fastp report not found for sample $sample. Skipping statistics extraction."
    fi
done

# Display the read statistics in a tabular format
echo "📊 Read Statistics (Processed Data):"
column -t readLOG.txt
echo ""

# ---------------------------
# Generate Read Statistics Graph with Plotext
# ---------------------------
plot_read_statistics

# ---------------------------
# Script Completion
# ---------------------------
display_step "$TITLE_COMPLETED"

echo "✅ FASTQ processing pipeline has been successfully executed."
echo "📄 Check 'readLOG.txt' for read statistics and the 'reports' directories for detailed reports."
echo "📈 A graph of the read statistics has been generated in the terminal."
