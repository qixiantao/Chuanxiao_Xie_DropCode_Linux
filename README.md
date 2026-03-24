# DropCode: A Pipeline for Gene Editing Site Variation Detection

DropCode is an automated analysis pipeline for next-generation sequencing (NGS) data, specifically designed for amplicon sequencing with dual-end barcodes. It integrates data preprocessing, quality control, demultiplexing, alignment, and variant calling, finally generating an allele frequency table for each sample.
<img width="1867" height="897" alt="f1589bf8-b758-424a-87f5-b084e02f435e" src="https://github.com/user-attachments/assets/0a53e97a-881e-4654-b9eb-846cb7ef1139" />

## Features

- **Preprocessing**: Automatically decompresses gzipped FASTQ files and extracts the upstream region of the target editing site.
- **Quality Filtering**: Uses `fastp` to trim low-quality bases (user-defined quality threshold, default Q20).
- **Demultiplexing**: Splits mixed sequencing data into individual samples based on dual-end barcodes provided in `barcode.xlsx` and an optional `library.xlsx` for name‑to‑password replacement.
- **Alignment**: Maps reads to the reference genome with `bwa mem`, filters with `samtools`, sorts, adds read groups with `Picard`, and generates sorted BAM files with indices.
- **Variant Calling**: Calls variants using `bcftools mpileup` and `bcftools call` (or a direct mpileup‑based method), calculates allele frequencies in the target region, and outputs an Excel table.
- **Batch Processing**: Supports automatic processing of multiple sample folders under a single input directory.
- **Error Handling & Logging**: Detailed logs are written to the `logs/` directory for easy troubleshooting.

## System Requirements

- Linux (Ubuntu 18.04+ recommended)
- Python 3.6+
- The following tools must be available in `PATH`:
  - `bwa`
  - `samtools`
  - `fastp`
  - `bcftools`
  - `java` (for Picard)
- Conda (Miniconda or Anaconda) for environment management

## Installation

1. **Clone the repository**
   ```bash
   git clone https://github.com/qixiantao/Chuanxiao_Xie_DropCode_Linux.git
   cd Chuanxiao_Xie_DropCode_Linux/DropCode
   ```

2. **Run the installation script**
   ```bash
   bash install.sh
   ```
   This script will:
   - Create a Conda environment named `dropcode` (if not present)
   - Install required tools (`bwa`, `samtools`, `fastp`, `bcftools`, `openjdk`)
   - Install Python dependencies (`pandas`, `openpyxl`, `biopython`, `vcfpy`)
   - Download Picard JAR into the `src/` directory
   - Make the main scripts executable

   **Note**: If you prefer to install system tools via `apt`, the script will attempt to use `apt` first and fall back to Conda.

3. **Verify the installation**
   ```bash
   conda activate dropcode
   bwa --version
   samtools --version
   fastp --version
   bcftools --version
   java -version
   ```

## Input File Preparation

Place each sample's data in a separate subfolder under `input_file` (the folder name will be used as the sample name). Each subfolder must contain:

- **`*1.fq.gz`** and **`*2.fq.gz`**: Forward and reverse reads (gzipped). The file names must end with `1.fq.gz` and `2.fq.gz`.
- **`reference.fasta`**: Reference genome sequence (single contig containing the target region).
- **`target.fasta`**: Target editing site sequence (single sequence; the script extracts its upstream/downstream region).
- **`barcode.xlsx`**: Barcode information. Format (first row is header):
  | Sample | Index1 | Index2 |
  |--------|--------|--------|
  | A      | ATCG   | GCTA   |
  | B      | CGAT   | TAGC   |
- **`library.xlsx`** (optional): Mapping from sample names to passwords, used to replace names in `barcode.xlsx` with passwords. Format:
  | Name | Password |
  |------|----------|
  | A    | pwd1     |
  | B    | pwd2     |

Directory structure example:

```
input_file/
├── sample1/
│   ├── sample1_1.fq.gz
│   ├── sample1_2.fq.gz
│   ├── reference.fasta
│   ├── target.fasta
│   ├── barcode.xlsx
│   └── library.xlsx (optional)
├── sample2/
│   └── ...
```

## Usage

### Process all samples in batch

```bash
bash batch_run.sh [--t THREADS] [--ram RAM_GB] [--q QUALITY] [--input INPUT_DIR] [--output OUTPUT_DIR]
```

Options:
- `--t`: Number of CPU threads (default: 1)
- `--ram`: Total memory in GB for sorting (default: 1)
- `--q`: Quality threshold for fastp and mpileup (10/20/30, default: 20)
- `--input`: Input folder path (default: `./input_file`)
- `--output`: Output folder path (default: `./output_file`)

Example:
```bash
bash batch_run.sh --t 8 --ram 4 --q 20
```

Results for each sample will be saved under `output_file/<sample_name>/`.

### Process a single sample

```bash
bash run_sample.sh --name SAMPLE_NAME --input SAMPLE_DIR --output OUTPUT_DIR [--t THREADS] [--ram RAM_GB] [--q QUALITY]
```

Example:
```bash
bash run_sample.sh --name test --input ./input_file/test --output ./output_file/test --t 8 --ram 4 --q 20
```

## Output

For each sample, the directory `output_file/<sample_name>/` contains:

- **`BAM/`**: Aligned BAM files (`*.sorted.bam`) and indices (`*.bai`).
- **`GEAnalysis_result.xlsx`**: Excel file with allele frequency statistics. The top 5 alleles and their percentages are listed.

Example Excel content:

| Sample | Allele1_Seq_Percent | Allele2_Seq_Percent | ... |
|--------|---------------------|---------------------|-----|
| test   | A: 85.23%           | C: 12.45%           | ... |

## Logs

All logs are stored in the `logs/` directory:
- `batch_run_YYYYMMDD_HHMMSS.log`: Master log for batch processing
- `run_sample_<sample_name>.log`: Detailed log for each sample
- `preprocess.log`, `qc_filter.log`, `demultiplex.log`, `align.log`, `variant_call.log`: Logs for individual steps

## FAQ

**Q: Tools like `bwa` are not found.**  
A: Ensure they are installed and in `PATH`. You can use `install.sh` to install them via Conda.

**Q: `picard.jar` not found.**  
A: Place `picard.jar` in the `src/` directory, or download it manually:
```bash
wget -O src/picard.jar https://github.com/broadinstitute/picard/releases/download/2.27.5/picard.jar
```

**Q: Many reads are unmatched during demultiplexing.**  
A: Check that the barcode sequences in `barcode.xlsx` are correct (10 bp for each end) and that the mapping in `library.xlsx` is accurate.

**Q: No variants are detected.**  
A: Possibly the target region is not covered or quality filtering is too strict. Verify that `target.fasta` exists in the reference sequence, or lower the `--q` parameter. Check the logs for coverage information.

**Q: `samtools sort` fails due to memory.**  
A: Reduce the `--ram` parameter or the number of threads (`--t`). The pipeline will automatically fall back to single‑thread sorting if multi‑thread fails.

**Q: How to see detailed error messages?**  
A: Check the relevant log file in the `logs/` directory or the console output; errors are logged at ERROR level.

## Citation

If you use DropCode in your research, please cite this repository:

```
DropCode: a pipeline for gene editing site variation detection. https://github.com/yourusername/DropCode
```
