#!/bin/bash
set -euo pipefail

# Activate conda environment
CONDA_BASE=$(conda info --base 2>/dev/null || echo "")
if [ -z "$CONDA_BASE" ]; then
    echo "Error: Conda not found. Please install Conda first."
    exit 1
fi
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate dropcode || { echo "Failed to activate dropcode environment"; exit 1; }

# Ensure library path is set (avoid missing libcrypto issues)
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH

# Default values
THREADS=1
RAM=1
QUAL=20
SAMPLE_NAME=""
INPUT_DIR=""
OUTPUT_DIR=""
LOG_DIR="./logs"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --name) SAMPLE_NAME="$2"; shift 2 ;;
        --input) INPUT_DIR="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        --t) THREADS="$2"; shift 2 ;;
        --ram) RAM="$2"; shift 2 ;;
        --q) QUAL="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

if [ -z "$SAMPLE_NAME" ] || [ -z "$INPUT_DIR" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing required arguments."
    exit 1
fi

# Setup logging for this sample
mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/run_sample_${SAMPLE_NAME}.log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "=== Starting sample $SAMPLE_NAME at $(date) ==="

# Find forward and reverse FASTQ files (ending with 1.fq.gz and 2.fq.gz)
FWD_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*1.fq.gz" | head -n 1)
REV_FILE=$(find "$INPUT_DIR" -maxdepth 1 -type f -name "*2.fq.gz" | head -n 1)

if [ -z "$FWD_FILE" ] || [ -z "$REV_FILE" ]; then
    echo "Error: Could not find files ending with 1.fq.gz and 2.fq.gz in $INPUT_DIR"
    exit 1
fi

echo "Forward file: $FWD_FILE"
echo "Reverse file: $REV_FILE"

# Check other required files
if [ ! -f "$INPUT_DIR/reference.fasta" ]; then
    echo "Error: reference.fasta not found in $INPUT_DIR"
    exit 1
fi
if [ ! -f "$INPUT_DIR/target.fasta" ]; then
    echo "Error: target.fasta not found in $INPUT_DIR"
    exit 1
fi
if [ ! -f "$INPUT_DIR/barcode.xlsx" ]; then
    echo "Error: barcode.xlsx not found in $INPUT_DIR"
    exit 1
fi

# Create temp and output dirs
TEMP_DIR="$OUTPUT_DIR/temp"
mkdir -p "$TEMP_DIR"

# Step 1: Preprocess
echo "Step 1: Extracting upstream sequence..."
python src/preprocess.py \
    --reference "$INPUT_DIR/reference.fasta" \
    --target "$INPUT_DIR/target.fasta" \
    --output "$TEMP_DIR/index.fasta" \
    --upstream 10 \
    --log-level INFO || { echo "Preprocess failed"; exit 1; }

# Step 2: QC filtering
echo "Step 2: Quality filtering..."
python src/qc_filter.py \
    --forward "$FWD_FILE" \
    --reverse "$REV_FILE" \
    --output_fwd "$TEMP_DIR/filtered_1.fastq" \
    --output_rev "$TEMP_DIR/filtered_2.fastq" \
    --quality "$QUAL" \
    --threads "$THREADS" \
    --log-level INFO || { echo "QC filtering failed"; exit 1; }

# Step 3: Demultiplexing
echo "Step 3: Demultiplexing..."
# Check for library.xlsx in input dir or src
LIBRARY_FILE=""
if [ -f "$INPUT_DIR/library.xlsx" ]; then
    LIBRARY_FILE="$INPUT_DIR/library.xlsx"
elif [ -f "src/library.xlsx" ]; then
    LIBRARY_FILE="src/library.xlsx"
else
    echo "Warning: library.xlsx not found, skipping name-to-password replacement."
fi

python src/demultiplex.py \
    --forward "$TEMP_DIR/filtered_1.fastq" \
    --reverse "$TEMP_DIR/filtered_2.fastq" \
    --barcode "$INPUT_DIR/barcode.xlsx" \
    ${LIBRARY_FILE:+--library "$LIBRARY_FILE"} \
    --output_dir "$TEMP_DIR/samples" \
    --log-level INFO || { echo "Demultiplexing failed"; exit 1; }

# Step 4: Alignment
echo "Step 4: Aligning to reference..."
python src/align.py \
    --samples_dir "$TEMP_DIR/samples" \
    --reference "$INPUT_DIR/reference.fasta" \
    --output_bam_dir "$OUTPUT_DIR/BAM" \
    --threads "$THREADS" \
    --ram "$RAM" \
    --picard_jar "src/picard.jar" \
    --log-level INFO || { echo "Alignment failed (some samples may have been skipped)"; }

# Check if any BAM files were generated
if [ ! -d "$OUTPUT_DIR/BAM" ] || [ -z "$(ls -A $OUTPUT_DIR/BAM/*.sorted.bam 2>/dev/null)" ]; then
    echo "Error: No BAM files were generated. Alignment may have failed for all samples."
    exit 1
fi

# Step 5: Variant calling
echo "Step 5: Calling variants and generating allele table..."
python src/variant_call.py \
    --bam_dir "$OUTPUT_DIR/BAM" \
    --reference "$INPUT_DIR/reference.fasta" \
    --target_fasta "$TEMP_DIR/index.fasta" \
    --output_xlsx "$OUTPUT_DIR/GEAnalysis_result.xlsx" \
    --threads "$THREADS" \
    --qual "$QUAL" \
    --log-level INFO || { echo "Variant calling failed"; exit 1; }

echo "=== Finished sample $SAMPLE_NAME at $(date) ==="