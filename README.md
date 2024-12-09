# FASTQ Processing Pipeline

![Pipeline Workflow](https://img.shields.io/badge/Pipeline-FASTQ%20Processing-brightgreen)

Automate the preprocessing and quality control of paired-end FASTQ data using **FastQC**, **MultiQC**, and **fastp**. This Bash script streamlines the entire workflow by creating necessary directory structures, allowing environment selection, logging read statistics, and generating informative read statistics graphs.

## Table of Contents

- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Overview](#pipeline-overview)
- [Configuration](#configuration)
- [Output](#output)
- [Troubleshooting](#troubleshooting)
- [Contributing](#contributing)
- [License](#license)
- [Contact](#contact)

## Features

- **Automated Workflow:** Seamlessly integrates FastQC, MultiQC, and fastp for comprehensive FASTQ data preprocessing and quality control.
- **Environment Selection:** Allows users to select from available Conda environments to ensure the correct tool versions and dependencies.
- **Dynamic Configuration:** Interactive prompts enable users to customize fastp parameters based on their specific requirements.
- **Comprehensive Reporting:** Generates detailed reports in both HTML and JSON formats, including aggregated MultiQC reports.
- **Read Statistics Logging:** Extracts and logs read statistics before and after filtering, providing insights into data quality and processing effectiveness.
- **Visualization:** Creates informative read statistics graphs using Python's `plotext` library, visualized directly in the terminal.
- **Error Handling:** Robust checks ensure all necessary tools are installed and valid user inputs are provided.

## Prerequisites

Before running the script, ensure that the following tools and dependencies are installed and accessible:

- **Bash:** Version 4 or higher.
- **Conda:** For environment management.
- **FastQC:** Quality control tool for high throughput sequence data.
- **MultiQC:** Aggregates results from multiple tools into a single report.
- **fastp:** Fast all-in-one FASTQ preprocessor.
- **jq:** Command-line JSON processor.
- **Python 3:** With `plotext` library installed (the script will install it if missing).

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/yourusername/fastq-processing-pipeline.git
   cd fastq-processing-pipeline
   ```

2. **Ensure Execution Permissions:**

   ```bash
   chmod +x fastq_processing_pipeline.sh
   ```

3. **Install Required Conda Environments:**

   The script will prompt you to select a Conda environment. Ensure you have the necessary environments created with the required tools installed.

   Example to create a new environment:

   ```bash
   conda create -n fastq_env fastqc multiqc fastp jq python=3.8
   conda activate fastq_env
   ```

## Usage

Execute the script from your terminal:

```bash
./fastq_processing_pipeline.sh
```

### Interactive Prompts

The script will guide you through several interactive prompts:

1. **Conda Environment Selection:**
   - Lists available Conda environments (excluding `base`).
   - Prompts you to select an environment by entering its corresponding number.

2. **Directory Structure:**
   - Asks if the necessary directories (`rawdata/` and `processed_data/`) are already created.
   - If not, it creates the required directory structure.

3. **FASTQ File Organization:**
   - Moves all `.fastq.gz` files from the current directory to the `rawdata/` directory.

4. **fastp Configuration:**
   - Interactive prompts to enable/disable and set parameters for:
     - Quality Filtering
     - Length Filtering
     - Low Complexity Filtering
     - Adapter Trimming
     - PolyG Tail Trimming
     - PolyX Tail Trimming
     - Deduplication

## Pipeline Overview

1. **Environment Setup:**
   - Sources the appropriate Conda initialization script.
   - Activates the user-selected Conda environment.

2. **Tool Verification:**
   - Checks for the presence of required tools: `fastqc`, `multiqc`, `fastp`, and `jq`.

3. **Directory Structure Creation:**
   - Sets up `rawdata/` for input FASTQ files and `processed_data/` for output files.
   - Creates subdirectories for reports generated by FastQC, MultiQC, and fastp.

4. **FASTQ File Organization:**
   - Moves all `.fastq.gz` files to the `rawdata/` directory.

5. **Paired-End File Identification:**
   - Identifies paired-end FASTQ files based on the naming pattern `*_plus_1_aaa.fastq.gz`.

6. **Quality Control with FastQC:**
   - Runs FastQC on raw FASTQ files.
   - Aggregates FastQC reports using MultiQC.

7. **Preprocessing with fastp:**
   - Configures fastp parameters based on user input.
   - Executes fastp to perform filtering, trimming, and optional deduplication.
   - Logs read statistics before and after processing.

8. **Post-Processing Quality Control:**
   - Runs FastQC on processed FASTQ files.
   - Aggregates FastQC reports using MultiQC.

9. **Statistics Extraction and Visualization:**
   - Extracts read statistics from fastp's JSON reports using `jq`.
   - Logs statistics into `readLOG.txt`.
   - Generates a read statistics graph displayed in the terminal using Python's `plotext`.

10. **Completion:**
    - Provides a summary of the pipeline execution and points to the generated reports.

## Configuration

During the pipeline execution, you'll be prompted to configure various fastp parameters:

- **Quality Filtering:**
  - **Qualified Quality PHRED:** Minimum PHRED score for a base to be considered high quality.
  - **Unqualified Percent Limit:** Maximum percentage of low-quality bases allowed per read.

- **Length Filtering:**
  - **Minimum Length Required:** Reads shorter than this length will be discarded.
  - **Maximum Length Limit:** Reads longer than this length will be discarded (0 means no limit).

- **Low Complexity Filtering:**
  - **Complexity Threshold:** Percentage of base diversity required to pass.

- **Adapter Trimming:**
  - Option to enable or disable adapter trimming.

- **PolyG and PolyX Tail Trimming:**
  - **PolyG Minimum Length:** Minimum length to detect and trim PolyG tails.
  - **PolyX Minimum Length:** Minimum length to detect and trim PolyX tails.

- **Deduplication:**
  - Option to enable deduplication to remove duplicated reads.

These configurations allow you to tailor the preprocessing steps to the specific requirements of your sequencing data.

## Output

### Directory Structure

After execution, the following directory structure will be created:

```
rawdata/
└── reports/
    ├── fastqc/
    └── multiqc/

processed_data/
└── reports/
    ├── fastqc/
    ├── multiqc/
    └── fastp/
```

### Reports

- **FastQC Reports:**
  - Located in `rawdata/reports/fastqc/` for raw data.
  - Located in `processed_data/reports/fastqc/` for processed data.

- **MultiQC Reports:**
  - Aggregated reports in `rawdata/reports/multiqc/` for raw data.
  - Aggregated reports in `processed_data/reports/multiqc/` for processed data.

- **fastp Reports:**
  - JSON and HTML reports for each sample in `processed_data/reports/fastp/`.
  - Log files capturing the fastp execution details.

### Read Statistics

- **readLOG.txt:**
  - Located in `processed_data/reports/fastp/`.
  - Logs the number of reads before and after filtering, as well as the number of discarded reads for each sample.

  **Sample Content:**

  ```
  Sample      Before_Reads    After_Reads    Discarded_Reads
  0_mM_NOD    1000000         950000          50000
  400_mM_NOD  1000000         960000          40000
  ```

- **Read Statistics Graph:**
  - Displayed directly in the terminal after execution.
  - Visual representation of read statistics using a multiple bar plot.

## Troubleshooting

- **Conda Activation Issues:**
  - Ensure that Conda is correctly installed and that the script is sourcing the correct `conda.sh` path.
  - Verify that the selected Conda environment has all required tools installed.

- **Missing Tools:**
  - The script checks for `fastqc`, `multiqc`, `fastp`, and `jq`. Ensure these are installed within the selected Conda environment.

- **Invalid FASTQ Files:**
  - Ensure that your FASTQ files follow the naming convention `*_plus_1_aaa.fastq.gz` and `*_plus_2_aaa.fastq.gz` for paired-end data.

- **Python Dependencies:**
  - The script uses `plotext` for generating graphs. If installation fails, manually install it using:

    ```bash
    pip install plotext
    ```

- **Insufficient Permissions:**
  - Ensure you have the necessary read/write permissions for the directories involved.

- **Error Messages During Execution:**
  - Review the corresponding log files in `processed_data/reports/fastp/` for detailed error information.

## Contributing

Contributions are welcome! If you have suggestions, bug reports, or enhancements, please open an issue or submit a pull request.

1. **Fork the Repository**

2. **Create a Feature Branch:**

   ```bash
   git checkout -b feature/YourFeature
   ```

3. **Commit Your Changes:**

   ```bash
   git commit -m "Add your feature"
   ```

4. **Push to the Branch:**

   ```bash
   git push origin feature/YourFeature
   ```

5. **Open a Pull Request**

## License

This project is licensed under the [MIT License](LICENSE).

## Contact

For any questions, issues, or suggestions, please contact [Daniel João](mailto:daniel.joao@example.com).

