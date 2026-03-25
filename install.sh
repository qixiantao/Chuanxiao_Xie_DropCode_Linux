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

echo "Installing required tools via conda..."
conda install -y -c bioconda bwa samtools bcftools fastp
conda install -y -c conda-forge openjdk=11

echo "Installing Python dependencies..."
pip install --upgrade pip
pip install pandas openpyxl biopython vcfpy

# Verify Python packages
echo "Verifying Python packages..."
python -c "import Bio; print('Biopython OK')" || { echo "Biopython installation failed"; exit 1; }
python -c "import pandas; print('pandas OK')" || { echo "pandas installation failed"; exit 1; }
python -c "import openpyxl; print('openpyxl OK')" || { echo "openpyxl installation failed"; exit 1; }
python -c "import vcfpy; print('vcfpy OK')" || { echo "vcfpy installation failed"; exit 1; }


# Reminder about other required files
echo "Note: Demo data must include reference.fasta, target.fasta, barcode.xlsx, and optionally library.xlsx."
echo "If they are missing, please copy them from the repository's 'demo' folder (if exists) or create them manually."

# ================== Finalize ==================
echo "Setting script permissions..."
chmod +x batch_run.sh run_sample.sh 2>/dev/null || true

echo "Setup completed successfully!"
echo "To manually activate the environment, run: conda activate $ENV_NAME"
echo "You can now run the pipeline using: bash batch_run.sh --t 8 --ram 4 --q 20"
