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

echo "Downloading Picard JAR..."
mkdir -p src
if [ ! -f "src/picard.jar" ]; then
    wget -O src/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
fi

echo "Setting script permissions..."
chmod +x batch_run.sh run_sample.sh 2>/dev/null || true

echo "Setup completed successfully!"
echo "To manually activate the environment, run: conda activate $ENV_NAME"
echo "You can now run the pipeline using: bash batch_run.sh --t 8 --ram 4 --q 20"