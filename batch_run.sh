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
FILTER_FREQ=0
INPUT_DIR="./input_file"
OUTPUT_DIR="./output_file"
LOG_DIR="./logs"

# Help message
usage() {
    echo "Usage: bash batch_run.sh [options]"
    echo ""
    echo "Options:"
    echo "  --t <int>         Number of threads (default: 1)"
    echo "  --ram <int>       Memory in GB (default: 1)"
    echo "  --q <int>         Quality threshold (default: 20)"
    echo "  --f <int>         Allele frequency filter threshold, 0-20 (default: 0)"
    echo "  --input <dir>     Input directory (default: ./input_file)"
    echo "  --output <dir>    Output directory (default: ./output_file)"
    echo "  -h, --help        Show this help message"
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --t)
            THREADS="$2"
            shift 2
            ;;
        --ram)
            RAM="$2"
            shift 2
            ;;
        --q)
            QUAL="$2"
            shift 2
            ;;
        --f)
            FILTER_FREQ="$2"
            shift 2
            ;;
        --input)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Validate numeric arguments
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [ "$THREADS" -lt 1 ]; then
    echo "Error: --t must be a positive integer."
    exit 1
fi

if ! [[ "$RAM" =~ ^[0-9]+$ ]] || [ "$RAM" -lt 1 ]; then
    echo "Error: --ram must be a positive integer."
    exit 1
fi

if ! [[ "$QUAL" =~ ^[0-9]+$ ]] || [ "$QUAL" -lt 0 ]; then
    echo "Error: --q must be a non-negative integer."
    exit 1
fi

if ! [[ "$FILTER_FREQ" =~ ^[0-9]+$ ]] || [ "$FILTER_FREQ" -lt 0 ] || [ "$FILTER_FREQ" -gt 20 ]; then
    echo "Error: --f must be an integer between 0 and 20."
    exit 1
fi

mkdir -p "$LOG_DIR"
LOG_FILE="$LOG_DIR/batch_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "=========================================="
echo " GEAtool Batch Processing"
echo "=========================================="
echo "Threads: $THREADS"
echo "Memory (GB): $RAM"
echo "Quality threshold: $QUAL"
echo "Allele frequency filter (--f): $FILTER_FREQ"
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Log file: $LOG_FILE"
echo "=========================================="

if [ ! -d "$INPUT_DIR" ]; then
    echo "Error: Input directory $INPUT_DIR does not exist!"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"

shopt -s nullglob
sample_found=0

for sample_dir in "$INPUT_DIR"/*/; do
    if [ -d "$sample_dir" ]; then
        sample_found=1
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
            --q "$QUAL" \
            --f "$FILTER_FREQ"

        echo ">>> Finished processing $sample_name"
    fi
done

if [ "$sample_found" -eq 0 ]; then
    echo "Error: No sample subdirectories found under $INPUT_DIR"
    exit 1
fi

echo ""
echo "All samples processed successfully!"