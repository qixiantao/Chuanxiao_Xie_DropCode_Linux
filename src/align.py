#!/usr/bin/env python3
"""
Align each sample to reference, sort, add RG, index.
"""
import argparse
import subprocess
import sys
import os
import logging
from pathlib import Path

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

def run_cmd(cmd, desc, logger, check=True):
    """Run a command, log output, raise on error if check=True."""
    logger.info(f"Running: {desc}")
    logger.debug(f"Command: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=check)
        if result.stdout:
            logger.debug(result.stdout)
        return result
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with code {e.returncode}: {e.stderr}")
        if check:
            raise
        return e

def bwa_index(ref_fasta, logger):
    idx_files = [f"{ref_fasta}.bwt", f"{ref_fasta}.sa", f"{ref_fasta}.pac"]
    if all(os.path.exists(f) for f in idx_files):
        logger.info("BWA index already exists, skipping.")
        return True
    logger.info("Building BWA index...")
    run_cmd(["bwa", "index", ref_fasta], "bwa index", logger)

def samtools_faidx(ref_fasta, logger):
    if os.path.exists(f"{ref_fasta}.fai"):
        logger.info("Fasta index already exists, skipping.")
        return True
    run_cmd(["samtools", "faidx", ref_fasta], "samtools faidx", logger)

def sort_bam(bam_path, sorted_bam_path, threads, total_mem_gb, temp_dir, logger):
    """Sort BAM with memory calculation and fallback to single-thread if needed."""
    # Calculate per-thread memory
    total_mem_mb = total_mem_gb * 1024
    per_thread_mem_mb = max(500, total_mem_mb // threads)
    per_thread_mem = f"{per_thread_mem_mb}M"
    logger.info(f"Trying multi-thread sort: threads={threads}, per-thread memory={per_thread_mem}")

    # Multi-thread command
    cmd = [
        "samtools", "sort",
        "-@", str(threads),
        "-m", per_thread_mem,
        "-T", str(temp_dir / "tmp"),
        str(bam_path),
        "-o", str(sorted_bam_path)
    ]
    result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
    if result.returncode == 0:
        logger.info("Multi-thread sort succeeded.")
        return True
    else:
        logger.warning(f"Multi-thread sort failed: {result.stderr}")
        logger.info("Falling back to single-thread sort (using all memory)...")
        # Single-thread command (use total memory)
        cmd = [
            "samtools", "sort",
            "-@", "1",
            "-m", f"{total_mem_mb}M",
            "-T", str(temp_dir / "tmp"),
            str(bam_path),
            "-o", str(sorted_bam_path)
        ]
        result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
        if result.returncode == 0:
            logger.info("Single-thread sort succeeded.")
            return True
        else:
            logger.error(f"Single-thread sort also failed: {result.stderr}")
            return False

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples_dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--output_bam_dir", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--ram", type=int, default=1, help="Total memory in GB for sorting")
    parser.add_argument("--picard_jar", required=True)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG","INFO","WARNING","ERROR"])
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper())
    logger = setup_logger("align", level=log_level)

    try:
        if not os.path.exists(args.reference):
            raise FileNotFoundError(f"Reference not found: {args.reference}")
        if not os.path.exists(args.picard_jar):
            raise FileNotFoundError(f"Picard JAR not found: {args.picard_jar}")

        bwa_index(args.reference, logger)
        samtools_faidx(args.reference, logger)

        samples_dir = Path(args.samples_dir)
        fwd_files = list(samples_dir.glob("*_1.fastq"))
        if not fwd_files:
            logger.warning(f"No *_1.fastq files found in {samples_dir}")
            return

        out_dir = Path(args.output_bam_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        temp_dir = out_dir / "tmp"
        temp_dir.mkdir(exist_ok=True)

        any_success = False
        for fwd in fwd_files:
            prefix = fwd.stem.replace("_1", "")
            rev = samples_dir / f"{prefix}_2.fastq"
            if not rev.exists():
                logger.warning(f"Missing reverse file for {prefix}, skipping.")
                continue
            if fwd.stat().st_size == 0 or rev.stat().st_size == 0:
                logger.warning(f"FASTQ file(s) for {prefix} are empty, skipping.")
                continue

            logger.info(f"Processing sample: {prefix}")

            sam = out_dir / f"{prefix}.sam"
            bam = out_dir / f"{prefix}.bam"
            sorted_bam = out_dir / f"{prefix}.sorted.bam"

            # BWA MEM
            cmd = ["bwa", "mem", "-t", str(args.threads), "-k", "15", "-M", args.reference, str(fwd), str(rev)]
            with open(sam, "w") as sam_out:
                result = subprocess.run(cmd, stdout=sam_out, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    logger.error(f"BWA failed for {prefix}: {result.stderr}")
                    continue
                if sam.stat().st_size == 0:
                    logger.error(f"SAM file is empty for {prefix}, alignment may have failed.")
                    continue
            logger.info(f"Alignment complete for {prefix}")

            # SAM -> BAM and filter
            view_opts = ["-b", "-f", "2", "-F", "4", "-F", "256", "-F", "2048", "-q", "20"]
            cmd = ["samtools", "view"] + view_opts + [str(sam)]
            with open(bam, "w") as bam_out:
                result = subprocess.run(cmd, stdout=bam_out, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    logger.error(f"samtools view failed for {prefix}: {result.stderr}")
                    continue
            logger.info(f"Converted to BAM and filtered for {prefix}")

            # Sort (with fallback)
            if not sort_bam(bam, sorted_bam, args.threads, args.ram, temp_dir, logger):
                logger.error(f"Sorting failed for {prefix}, skipping.")
                continue

            # Clean intermediate
            sam.unlink()
            bam.unlink()
            logger.debug(f"Removed intermediate files for {prefix}")

            # Add read groups
            temp_bam = out_dir / f"{prefix}.sorted.rg.bam"
            cmd = [
                "java", "-jar", args.picard_jar, "AddOrReplaceReadGroups",
                f"I={sorted_bam}", f"O={temp_bam}",
                f"RGID={prefix}", "RGLB=lib1", "RGPL=ILLUMINA", "RGPU=unit1", f"RGSM={prefix}"
            ]
            result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                logger.error(f"Picard failed for {prefix}: {result.stderr}")
                continue
            sorted_bam.unlink()
            temp_bam.rename(sorted_bam)

            # Index
            cmd = ["samtools", "index", str(sorted_bam)]
            result = subprocess.run(cmd, stderr=subprocess.PIPE, text=True)
            if result.returncode != 0:
                logger.error(f"samtools index failed for {prefix}: {result.stderr}")
                continue
            logger.info(f"Finished processing {prefix}")
            any_success = True

        if not any_success:
            logger.error("No sample was successfully aligned. Check logs for details.")
            sys.exit(1)
        logger.info("All samples aligned successfully.")

    except Exception as e:
        logger.error(f"Alignment pipeline failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()