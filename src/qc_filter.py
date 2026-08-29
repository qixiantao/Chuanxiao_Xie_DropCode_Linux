#!/usr/bin/env python3
"""
Decompress raw gzipped FASTQ and run fastp for quality trimming.
"""
import argparse
import gzip
import shutil
import subprocess
import sys
import os
import logging
from datetime import datetime

def setup_logger(name, log_dir="./logs", level=logging.INFO):
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    log_file = os.path.join(log_dir, f"{name}.log")
    fh = logging.FileHandler(log_file)
    fh.setLevel(level)
    ch = logging.StreamHandler()
    ch.setLevel(level)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def decompress_gz(gz_path, output_path, logger):
    """Decompress gzip file."""
    try:
        with gzip.open(gz_path, 'rb') as fin, open(output_path, 'wb') as fout:
            shutil.copyfileobj(fin, fout)
        logger.info(f"Decompressed {gz_path} -> {output_path}")
    except Exception as e:
        logger.error(f"Failed to decompress {gz_path}: {e}")
        raise

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--forward", required=True, help="Forward gzipped FASTQ")
    parser.add_argument("--reverse", required=True, help="Reverse gzipped FASTQ")
    parser.add_argument("--output_fwd", required=True, help="Output filtered forward FASTQ")
    parser.add_argument("--output_rev", required=True, help="Output filtered reverse FASTQ")
    parser.add_argument("--quality", type=int, default=20, help="Quality threshold for fastp")
    parser.add_argument("--threads", type=int, default=1, help="Number of threads")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger("qc_filter", level=log_level)

    try:
        # Create temp directory for decompressed files
        temp_dir = os.path.dirname(args.output_fwd)
        os.makedirs(temp_dir, exist_ok=True)

        fwd_fq = os.path.join(temp_dir, "raw_1.fastq")
        rev_fq = os.path.join(temp_dir, "raw_2.fastq")

        logger.info("Decompressing forward reads...")
        decompress_gz(args.forward, fwd_fq, logger)
        logger.info("Decompressing reverse reads...")
        decompress_gz(args.reverse, rev_fq, logger)

        # Run fastp
        cmd = [
            "fastp",
            "-i", fwd_fq,
            "-I", rev_fq,
            "-o", args.output_fwd,
            "-O", args.output_rev,
            "-q", str(args.quality),
            "-l", "70",
            "-h", os.path.join(temp_dir, "fastp_report.html"),
            "-j", os.path.join(temp_dir, "fastp_report.json"),
            "-w", str(args.threads)
        ]
        logger.info(f"Running fastp: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"fastp failed with code {result.returncode}: {result.stderr}")
            sys.exit(1)
        logger.info("fastp completed successfully.")

        # Clean up temporary raw files
        os.remove(fwd_fq)
        os.remove(rev_fq)
        logger.info("Removed temporary decompressed files.")
        logger.info("QC filtering completed.")

    except Exception as e:
        logger.error(f"QC filtering failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()