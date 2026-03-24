#!/usr/bin/env python3
"""
Extract upstream sequence from reference based on target.
"""
import argparse
import logging
import os
import sys
from Bio import SeqIO
from datetime import datetime

def setup_logger(name, log_dir="./logs", level=logging.INFO):
    """Setup logger with file and console handlers."""
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # File handler
    log_file = os.path.join(log_dir, f"{name}.log")
    fh = logging.FileHandler(log_file)
    fh.setLevel(level)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(level)

    # Formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

def convert_to_uppercase(filepath, logger):
    """Convert sequence lines to uppercase in place."""
    try:
        tmp = filepath + ".tmp"
        with open(filepath, 'r') as fin, open(tmp, 'w') as fout:
            for line in fin:
                if line.startswith('>'):
                    fout.write(line)
                else:
                    fout.write(line.upper())
        os.replace(tmp, filepath)
        logger.info(f"Converted {filepath} to uppercase.")
    except Exception as e:
        logger.error(f"Failed to convert {filepath}: {e}")
        raise

def extract_upstream(reference_fasta, target_fasta, output_fasta, upstream_len, logger):
    """Extract upstream (or reverse complement downstream) sequence."""
    try:
        logger.info("Converting reference and target to uppercase...")
        convert_to_uppercase(reference_fasta, logger)
        convert_to_uppercase(target_fasta, logger)

        logger.info("Reading sequences...")
        target = str(next(SeqIO.parse(target_fasta, "fasta")).seq)
        ref = str(next(SeqIO.parse(reference_fasta, "fasta")).seq)

        pos = ref.find(target)
        if pos == -1:
            raise ValueError("Target sequence not found in reference!")

        logger.info(f"Target found at position {pos}.")
        if pos < len(ref) // 2:
            # Upstream
            start = pos - upstream_len
            if start < 0:
                raise ValueError(f"Upstream region extends beyond reference start (needs {upstream_len} bases).")
            extracted = ref[start:pos]
            logger.info("Extracted upstream sequence.")
        else:
            # Downstream, reverse complement
            start = pos + len(target)
            end = start + upstream_len
            if end > len(ref):
                raise ValueError(f"Downstream region extends beyond reference end (needs {upstream_len} bases).")
            extracted = ref[start:end]
            extracted = str(SeqIO.Seq(extracted).reverse_complement())
            logger.info("Extracted downstream sequence and reverse complemented.")

        os.makedirs(os.path.dirname(output_fasta), exist_ok=True)
        with open(output_fasta, 'w') as f:
            f.write(f">upstream_seq\n{extracted}\n")
        logger.info(f"Extracted sequence written to {output_fasta}")

    except Exception as e:
        logger.error(f"Extraction failed: {e}")
        raise

def main():
    parser = argparse.ArgumentParser(description="Extract upstream sequence from reference.")
    parser.add_argument("--reference", required=True, help="Reference FASTA file")
    parser.add_argument("--target", required=True, help="Target FASTA file")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    parser.add_argument("--upstream", type=int, default=10, help="Number of bases to extract")
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"], help="Logging level")
    args = parser.parse_args()

    # Setup logger
    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger("preprocess", level=log_level)

    try:
        extract_upstream(args.reference, args.target, args.output, args.upstream, logger)
        logger.info("Preprocessing completed successfully.")
    except Exception:
        sys.exit(1)

if __name__ == "__main__":
    main()