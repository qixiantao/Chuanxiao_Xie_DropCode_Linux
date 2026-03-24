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

# Default parameters
THREADS=1
RAM=1
QUAL=20
INPUT_DIR="./input_file"
OUTPUT_DIR="./output_file"
LOG_DIR="./logs"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --t) THREADS="$2"; shift 2 ;;
        --ram) RAM="$2"; shift 2 ;;
        --q) QUAL="$2"; shift 2 ;;
        --input) INPUT_DIR="$2"; shift 2 ;;
        --output) OUTPUT_DIR="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/batch_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo " GEAtool Batch Processing"
echo "=========================================="
echo "Threads: $THREADS"
echo "Memory (GB): $RAM"
echo "Quality threshold: $QUAL"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
echo "=========================================="

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist!"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

for sample_dir in "$INPUT_DIR"/*/; do
    if [ -d "$sample_dir" ]; then
        sample_name=$(basename "$sample_dir")
        echo ""
        echo ">>> Processing sample: $sample_name"

        sample_output="$OUTPUT_DIR/$sample_name"
        mkdir -p "$sample_output"

        # Use bash to execute run_sample.sh to avoid permission issues
        bash ./run_sample.sh \
            --name "$sample_name" \
            --input "$sample_dir" \
            --output "$sample_output" \
            --t "$THREADS" \
            --ram "$RAM" \
            --q "$QUAL"

        echo ">>> Finished processing $sample_name"
    fi
done

echo ""
echo "All samples processed successfully!"