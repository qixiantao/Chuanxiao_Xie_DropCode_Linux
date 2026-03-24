#!/bin/bash
set -euo pipefail

ENV_NAME="dropcode"
CONDA_BASE=$(conda info --base 2>/dev/null || echo "")

if [ -z "$CONDA_BASE" ]; then
    echo "Error: Conda is not installed or not in PATH."
    echo "Please install Miniconda first: https://docs.conda.io/en/latest/miniconda.html"
    exit 1
fi

echo "Using Conda at: $CONDA_BASE"

# Create environment if not exists
if ! conda env list | grep -q "^$ENV_NAME "; then
    echo "Creating conda environment: $ENV_NAME"
    conda create -y -n "$ENV_NAME" python=3.9
fi

# Activate environment
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

# Function to check if tool is installed
check_tool() {
    command -v "$1" >/dev/null 2>&1
}

# Function to install via apt if possible
install_with_apt() {
    if command -v apt-get >/dev/null 2>&1; then
        echo "Installing $1 via apt..."
        sudo apt-get update -qq
        sudo apt-get install -y "$1"
        return 0
    else
        return 1
    fi
}

# Define required tools
TOOLS="bwa samtools bcftools openjdk-11-jre"
# Note: fastp may not be available in default repos; will be handled separately

# Try apt first for standard tools
MISSING_APT=""
for tool in $TOOLS; do
    if check_tool "$tool"; then
        echo "$tool already installed."
    else
        if install_with_apt "$tool"; then
            echo "$tool installed via apt."
        else
            echo "apt not available or installation failed; will use conda later."
            MISSING_APT="$MISSING_APT $tool"
        fi
    fi
done

# Handle fastp specially (may need conda)
if check_tool fastp; then
    echo "fastp already installed."
else
    if command -v apt-get >/dev/null 2>&1 && sudo apt-get install -y fastp 2>/dev/null; then
        echo "fastp installed via apt."
    else
        echo "fastp not found in apt; will use conda."
        MISSING_APT="$MISSING_APT fastp"
    fi
fi

# Install any remaining tools with conda
if [ -n "$MISSING_APT" ]; then
    echo "Installing missing tools via conda:$MISSING_APT"
    conda install -y -c bioconda $MISSING_APT
    # Ensure java is installed (if openjdk-11-jre not available, use conda)
    if ! check_tool java; then
        conda install -y -c conda-forge openjdk=11
    fi
else
    # Ensure java is available (if not, install via conda)
    if ! check_tool java; then
        echo "Java not found; installing via conda..."
        conda install -y -c conda-forge openjdk=11
    fi
fi

echo "Installing Python dependencies..."
pip install --upgrade pip
pip install pandas openpyxl biopython vcfpy

# Verify Python packages
echo "Verifying Python packages..."
python -c "import Bio; print('Biopython OK')" || { echo "Biopython installation failed"; exit 1; }
python -c "import pandas; print('pandas OK')" || { echo "pandas installation failed"; exit 1; }
python -c "import openpyxl; print('openpyxl OK')" || { echo "openpyxl installation failed"; exit 1; }
python -c "import vcfpy; print('vcfpy OK')" || { echo "vcfpy installation failed"; exit 1; }

# ================== Download Picard JAR ==================
echo "Downloading Picard JAR..."
mkdir -p src
if [ ! -f "src/picard.jar" ]; then
    echo "Attempting to download from GitHub release..."
    # curl 
    curl -L --retry 5 --connect-timeout 30 -o src/picard.jar \
        https://github.com/qixiantao/Chuanxiao_Xie_DropCode_Linux/releases/download/DropCode_picard.jar/picard.jar || {
        echo "GitHub release download failed. Trying official Picard release..."
        # wget
        wget --tries=5 --timeout=30 -O src/picard.jar \
            https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar || {
            echo "ERROR: Could not download Picard JAR. Please download manually and place in src/"
            exit 1
        }
    }
    echo "Picard JAR downloaded."
else
    echo "Picard JAR already exists, skipping download."
fi

# ================== Download DEMO data ==================
echo "Setting up DEMO data..."
DEMO_DIR="input_file/DEMO"
mkdir -p "$DEMO_DIR"

# Define expected demo files (we will check for any file ending with 1.fq.gz and 2.fq.gz)
if [ -z "$(ls -A "$DEMO_DIR"/*1.fq.gz 2>/dev/null)" ] || [ -z "$(ls -A "$DEMO_DIR"/*2.fq.gz 2>/dev/null)" ]; then
    echo "Downloading demo FASTQ files from GitHub release..."
    # Download forward read
    curl -L --retry 5 --connect-timeout 30 -o "$DEMO_DIR/demo_1.fq.gz" \
        https://github.com/qixiantao/Chuanxiao_Xie_DropCode_Linux/releases/download/DropCode_picard.jar/202603091635_AE01-250302016_4P251123103US293267A2_A_20260309_MWBMWB0309_L02_R1.fq.gz || {
        echo "Failed to download forward demo data. Please check network."
        exit 1
    }
    # Download reverse read
    curl -L --retry 5 --connect-timeout 30 -o "$DEMO_DIR/demo_2.fq.gz" \
        https://github.com/qixiantao/Chuanxiao_Xie_DropCode_Linux/releases/download/DropCode_picard.jar/202603091635_AE01-250302016_4P251123103US293267A2_A_20260309_MWBMWB0309_L02_R2.fq.gz || {
        echo "Failed to download reverse demo data. Please check network."
        exit 1
    }
    echo "Demo FASTQ files downloaded."
else
    echo "Demo FASTQ files already exist in $DEMO_DIR, skipping download."
fi

# Note: The demo folder should also contain reference.fasta, target.fasta, barcode.xlsx, etc.
# These are not provided in the release; they are expected to be present in the repository or copied manually.
# We'll print a reminder.
echo "Note: Demo data must include reference.fasta, target.fasta, barcode.xlsx, and optionally library.xlsx."
echo "If they are missing, please copy them from the repository's 'demo' folder (if exists) or create them manually."

# ================== Finalize ==================
echo "Setting script permissions..."
chmod +x batch_run.sh run_sample.sh 2>/dev/null || true

echo "Setup completed successfully!"
echo "To manually activate the environment, run: conda activate $ENV_NAME"
echo "You can now run the pipeline using: bash batch_run.sh --t 8 --ram 4 --q 20"