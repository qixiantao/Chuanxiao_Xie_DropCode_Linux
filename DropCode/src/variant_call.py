#!/usr/bin/env python3
"""
Call variants using samtools mpileup, extract alleles after target sequence.
"""
import argparse
import subprocess
import pandas as pd
import os
import sys
import logging
from pathlib import Path
from Bio import SeqIO
from collections import Counter

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

def get_reference_name(reference_fasta):
    """Extract chromosome name from FASTA header."""
    with open(reference_fasta) as f:
        for line in f:
            if line.startswith('>'):
                # header format: >ref_name [optional description]
                parts = line[1:].strip().split()
                return parts[0]
    raise ValueError("No header found in reference FASTA")

def get_target_region(reference_fasta, target_fasta, logger):
    """Return (chrom, start, end) of the target site in reference coordinates."""
    try:
        target = str(next(SeqIO.parse(target_fasta, "fasta")).seq).upper()
        ref = str(next(SeqIO.parse(reference_fasta, "fasta")).seq).upper()
        pos = ref.find(target)
        if pos == -1:
            raise ValueError("Target sequence not found in reference.")
        start = max(0, pos - 20)
        end = min(len(ref), pos + len(target) + 20)
        chrom = get_reference_name(reference_fasta)
        logger.info(f"Target region: {chrom}:{start}-{end}")
        return (chrom, start, end)
    except Exception as e:
        logger.error(f"Failed to locate target region: {e}")
        raise

def ensure_bam_index(bam_path, logger):
    if not os.path.exists(f"{bam_path}.bai"):
        logger.info(f"Index not found, generating for {bam_path}")
        cmd = ["samtools", "index", bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to index {bam_path}: {result.stderr}")
            return False
        logger.info("Index created.")
    return True

def extract_alleles_from_pileup(mpileup_output, target_seq, logger):
    """
    Parse mpileup output and extract the 20bp sequence following target.
    mpileup output: each line corresponds to a position; we need to walk through
    the target region and reconstruct reads.
    However, mpileup gives per-base depth, not read sequences.
    Alternative: use samtools view to fetch reads covering target region.
    """
    # We'll use a simpler method: extract read sequences via samtools view
    # For speed, we may skip this and just count base frequencies as before.
    # But your requirement is to get the 20bp allele sequence after target.
    # Since mpileup doesn't give read sequences, we need to fetch reads.
    pass

def fetch_alleles_from_bam(bam_path, chrom, start, end, target_seq, extract_len=20, logger=None):
    """
    Fetch all reads covering the target region, extract the sequence after target.
    """
    # Use samtools view to get reads overlapping the region
    region = f"{chrom}:{start}-{end}"
    cmd = ["samtools", "view", "-F", "0x904", bam_path, region]  # filter out unmapped, secondary, supplementary
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"Failed to fetch reads: {result.stderr}")
        return []
    alleles = []
    # Each line: QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL
    for line in result.stdout.splitlines():
        fields = line.split('\t')
        if len(fields) < 10:
            continue
        seq = fields[9]
        pos = int(fields[3])  # 1-based
        # Find target in read sequence (approximate, may need to account for indels)
        # We'll search for the target sequence in the read
        if target_seq in seq:
            idx = seq.find(target_seq)
            allele = seq[idx + len(target_seq): idx + len(target_seq) + extract_len]
            alleles.append(allele)
    return alleles

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--target_fasta", required=True)
    parser.add_argument("--output_xlsx", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--qual", type=int, default=20)
    parser.add_argument("--extract_len", type=int, default=20, help="Length of allele to extract after target")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger("variant_call", level=log_level)

    try:
        chrom, start, end = get_target_region(args.reference, args.target_fasta, logger)
        target_seq = str(next(SeqIO.parse(args.target_fasta, "fasta")).seq).upper()
        logger.info(f"Target sequence: {target_seq}")

        bam_files = list(Path(args.bam_dir).glob("*.sorted.bam"))
        if not bam_files:
            raise FileNotFoundError(f"No BAM files found in {args.bam_dir}")

        results = []
        for bam in bam_files:
            sample = bam.stem.replace(".sorted", "")
            logger.info(f"Processing {sample}...")

            if not ensure_bam_index(str(bam), logger):
                logger.error(f"Skipping {sample} due to missing index.")
                continue

            # Fetch alleles from BAM
            alleles = fetch_alleles_from_bam(str(bam), chrom, start, end, target_seq, args.extract_len, logger)
            if not alleles:
                logger.warning(f"No reads containing target found for {sample}")
                top = [("", 0)] * 5
            else:
                counter = Counter(alleles)
                total = len(alleles)
                items = sorted(counter.items(), key=lambda x: x[1], reverse=True)
                top = [(seq, cnt/total*100) for seq, cnt in items[:5]]
                top += [("", 0)] * (5 - len(top))
                logger.info(f"Top alleles: {top}")

            row = [sample] + [f"{seq}: {pct:.2f}%" for seq, pct in top]
            results.append(row)

        if not results:
            logger.warning("No samples processed successfully.")
        else:
            columns = ["Sample"] + [f"Allele{i+1}_Seq_Percent" for i in range(5)]
            df = pd.DataFrame(results, columns=columns)
            os.makedirs(os.path.dirname(args.output_xlsx), exist_ok=True)
            df.to_excel(args.output_xlsx, index=False)
            logger.info(f"Results written to {args.output_xlsx}")

    except Exception as e:
        logger.error(f"Variant calling failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()