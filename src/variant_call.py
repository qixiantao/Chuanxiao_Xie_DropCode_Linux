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
    """Create and return a logger."""
    os.makedirs(log_dir, exist_ok=True)
    logger = logging.getLogger(name)
    logger.setLevel(level)

    # Avoid duplicated handlers when the script is imported or re-run in the same process
    if logger.handlers:
        return logger

    log_file = os.path.join(log_dir, f"{name}.log")
    fh = logging.FileHandler(log_file)
    fh.setLevel(level)

    ch = logging.StreamHandler()
    ch.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)

    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger


def get_reference_name(reference_fasta):
    """Extract the reference/chromosome name from the FASTA header."""
    with open(reference_fasta) as f:
        for line in f:
            if line.startswith(">"):
                parts = line[1:].strip().split()
                return parts[0]
    raise ValueError("No header found in reference FASTA.")


def get_target_region(reference_fasta, target_fasta, logger):
    """
    Return (chrom, start, end, target_start, target_end, reference_sequence).

    start and end are expanded by 20 bp around the target region for BAM fetching.
    target_start and target_end are 0-based coordinates in the reference sequence.
    """
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
        return chrom, start, end, pos, pos + len(target), ref
    except Exception as e:
        logger.error(f"Failed to locate target region: {e}")
        raise


def get_wildtype_allele(reference_seq, target_end, extract_len):
    """
    Extract the expected wild-type allele sequence immediately downstream of target.
    """
    return reference_seq[target_end: target_end + extract_len].upper()


def ensure_bam_index(bam_path, logger):
    """Ensure BAM index exists."""
    if not os.path.exists(f"{bam_path}.bai"):
        logger.info(f"Index not found, generating for {bam_path}")
        cmd = ["samtools", "index", bam_path]
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"Failed to index {bam_path}: {result.stderr}")
            return False
        logger.info("Index created.")
    return True


def fetch_alleles_from_bam(bam_path, chrom, start, end, target_seq, extract_len=20, logger=None):
    """
    Fetch all reads overlapping the target region and extract the sequence
    immediately after the target sequence within each read.

    Notes:
    - This follows the original logic and searches the target sequence directly
      in the read sequence.
    - Reads not containing the exact target sequence are ignored.
    """
    region = f"{chrom}:{start}-{end}"
    cmd = ["samtools", "view", "-F", "0x904", bam_path, region]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        if logger:
            logger.error(f"Failed to fetch reads: {result.stderr}")
        return []

    alleles = []

    for line in result.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) < 10:
            continue

        seq = fields[9].upper()

        if target_seq in seq:
            idx = seq.find(target_seq)
            allele = seq[idx + len(target_seq): idx + len(target_seq) + extract_len]
            alleles.append(allele)

    return alleles


def validate_filter_threshold(value):
    """Validate filter threshold for --f."""
    ivalue = int(value)
    if ivalue < 0 or ivalue > 20:
        raise argparse.ArgumentTypeError("The --f value must be an integer between 0 and 20.")
    return ivalue


def get_mutation_category(allele_seq, wt_seq):
    """
    Infer a simple mutation category relative to the wild-type allele.

    Categories:
    - WT
    - Insertion
    - Deletion
    - Substitution
    - Complex
    """
    if allele_seq == wt_seq:
        return "WT"

    if len(allele_seq) > len(wt_seq):
        return "Insertion"
    elif len(allele_seq) < len(wt_seq):
        return "Deletion"
    else:
        mismatch_count = sum(1 for a, b in zip(allele_seq, wt_seq) if a != b)
        if mismatch_count == 1:
            return "Substitution"
        else:
            return "Complex"


def apply_allele_filter(sorted_items, total_reads, filter_threshold):
    """
    Apply allele frequency filtering based on the original total target reads.

    Parameters
    ----------
    sorted_items : list of (allele_seq, count)
        Alleles sorted by descending count.
    total_reads : int
        Total number of target-containing reads before filtering.
    filter_threshold : int
        Minimum percentage threshold (0-20).

    Returns
    -------
    filtered_items : list of dict
        Each item contains:
        {
            "seq": allele sequence,
            "count": count,
            "raw_pct": percentage based on total_reads before filtering
        }
    """
    filtered_items = []

    for seq, cnt in sorted_items:
        raw_pct = (cnt / total_reads * 100.0) if total_reads > 0 else 0.0
        if raw_pct >= filter_threshold:
            filtered_items.append({
                "seq": seq,
                "count": cnt,
                "raw_pct": raw_pct
            })

    return filtered_items


def recalculate_percentages(filtered_items):
    """
    Recalculate allele percentages using only the filtered reads.

    Parameters
    ----------
    filtered_items : list of dict
        Output from apply_allele_filter().

    Returns
    -------
    recalculated_items : list of dict
        Each item contains:
        {
            "seq": allele sequence,
            "count": count,
            "raw_pct": percentage before filtering,
            "filtered_pct": percentage after filtering
        }
    filtered_total_reads : int
        Sum of counts for retained alleles.
    """
    filtered_total_reads = sum(item["count"] for item in filtered_items)

    recalculated_items = []
    for item in filtered_items:
        filtered_pct = (
            item["count"] / filtered_total_reads * 100.0
            if filtered_total_reads > 0 else 0.0
        )
        recalculated_items.append({
            "seq": item["seq"],
            "count": item["count"],
            "raw_pct": item["raw_pct"],
            "filtered_pct": filtered_pct
        })

    return recalculated_items, filtered_total_reads


def infer_genotype_from_filtered_alleles(filtered_items, wt_allele):
    """
    Infer genotype based on the user-defined rules, using only the alleles retained
    after filtering.

    Rules implemented:
    1. If no allele remains after filtering:
       - No_Alleles_After_Filter

    2. If only one allele remains:
       - If reads1 == WT -> WildType
       - Else -> Homozygous_Mutant

    3. If exactly two alleles remain:
       - If one of reads1 or reads2 is WT:
         - If reads2 > 10% -> Heterozygous_Mutant
         - Else -> Chimeric_or_WT_with_PCR_Contamination

       - If both are non-WT:
         - If reads2 > 10%:
           - If both have the same mutation category -> HomozygousMutant
           - Else -> Biallelic_Mutant
         - Else -> Chimeric_Or_Homozygous_Mutant_with_PCR_Contamination

    4. If three or more alleles remain:
       - Chimeric

    Notes
    -----
    - reads2 percentage is evaluated based on the filtered allele total,
      as requested for filter-aware statistics.
    """
    if not filtered_items:
        return "No_Alleles_After_Filter"

    allele_count = len(filtered_items)

    if allele_count >= 3:
        return "Chimeric"

    seq1 = filtered_items[0]["seq"]
    cnt1 = filtered_items[0]["count"]
    pct1 = filtered_items[0]["filtered_pct"]
    is_wt1 = (seq1 == wt_allele)
    cat1 = get_mutation_category(seq1, wt_allele)

    if allele_count == 1:
        if is_wt1:
            return "Wild Type"
        else:
            return "Homozygous Mutant"

    seq2 = filtered_items[1]["seq"]
    cnt2 = filtered_items[1]["count"]
    pct2 = filtered_items[1]["filtered_pct"]
    is_wt2 = (seq2 == wt_allele)
    cat2 = get_mutation_category(seq2, wt_allele)

    # Exactly two alleles remain after filtering
    if is_wt1 or is_wt2:
        if pct2 > 10.0:
            return "Heterozygous Mutant"
        else:
            return "Chimeric or WT with PCR Contamination"
    else:
        if pct2 > 10.0:
            if cat1 == cat2:
                return "Homozygous Mutant"
            else:
                return "Biallelic Mutant"
        else:
            return "Chimeric Or Homozygous with PCR Contamination"


def format_allele_cell(seq, count, pct):
    """Format one allele cell for Excel output."""
    return f"{seq}: {count} reads, {pct:.2f}%"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bam_dir", required=True)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--target_fasta", required=True)
    parser.add_argument("--output_xlsx", required=True)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--qual", type=int, default=20)
    parser.add_argument(
        "--extract_len",
        type=int,
        default=20,
        help="Length of allele to extract after target"
    )
    parser.add_argument(
        "--f",
        type=validate_filter_threshold,
        default=0,
        help="Optional allele frequency filter threshold in percent (0-20). Default: 0"
    )
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args()

    log_level = getattr(logging, args.log_level.upper(), logging.INFO)
    logger = setup_logger("variant_call", level=log_level)

    try:
        chrom, start, end, target_start, target_end, reference_seq = get_target_region(
            args.reference, args.target_fasta, logger
        )
        target_seq = str(next(SeqIO.parse(args.target_fasta, "fasta")).seq).upper()
        wt_allele = get_wildtype_allele(reference_seq, target_end, args.extract_len)

        logger.info(f"Target sequence: {target_seq}")
        logger.info(f"Expected wild-type downstream allele: {wt_allele}")
        logger.info(f"Filter threshold: {args.f}%")

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

            alleles = fetch_alleles_from_bam(
                str(bam), chrom, start, end, target_seq, args.extract_len, logger
            )

            if not alleles:
                logger.warning(f"No reads containing target found for {sample}")

                total_target_reads = 0
                filtered_target_reads = 0
                filtered_reads_percent_of_total = 0.0
                genotype = "NoTargetReads"
                top_display = [("", 0, 0.0)] * 5

            else:
                counter = Counter(alleles)
                total_target_reads = len(alleles)

                # Sort alleles by descending count
                sorted_items = sorted(counter.items(), key=lambda x: x[1], reverse=True)

                # Step 1: apply filter based on original total reads
                filtered_items = apply_allele_filter(
                    sorted_items=sorted_items,
                    total_reads=total_target_reads,
                    filter_threshold=args.f
                )

                # Step 2: recalculate percentages based on filtered reads only
                recalculated_items, filtered_target_reads = recalculate_percentages(filtered_items)

                # Step 3: genotype inference based on filtered allele set
                genotype = infer_genotype_from_filtered_alleles(
                    filtered_items=recalculated_items,
                    wt_allele=wt_allele
                )

                filtered_reads_percent_of_total = (
                    filtered_target_reads / total_target_reads * 100.0
                    if total_target_reads > 0 else 0.0
                )

                # Only display the top 5 alleles after filtering
                top_display = [
                    (item["seq"], item["count"], item["filtered_pct"])
                    for item in recalculated_items[:5]
                ]
                top_display += [("", 0, 0.0)] * (5 - len(top_display))

                logger.info(f"Genotype: {genotype}")
                logger.info(
                    "Top alleles before filter: "
                    f"{[(seq, cnt, cnt / total_target_reads * 100.0) for seq, cnt in sorted_items[:10]]}"
                )
                logger.info(
                    "Alleles retained after filter and recalculated: "
                    f"{[(item['seq'], item['count'], item['filtered_pct']) for item in recalculated_items[:10]]}"
                )

            row = [
                sample,
                genotype,
                wt_allele,
                total_target_reads,
                args.f,
                filtered_target_reads,
                f"{filtered_reads_percent_of_total:.2f}%"
            ]

            for seq, cnt, pct in top_display:
                if seq:
                    row.append(format_allele_cell(seq, cnt, pct))
                else:
                    row.append("")

            results.append(row)

        if not results:
            logger.warning("No samples processed successfully.")
        else:
            columns = [
                "Sample",
                "Genotype",
                "WildType_Allele",
                "Total_Target_Reads",
                "Filter_Threshold_Percent",
                "Reads_Used_After_Filter",
                "Reads_Used_Percent_Of_Total",
                "Allele1_Seq_Count_Percent",
                "Allele2_Seq_Count_Percent",
                "Allele3_Seq_Count_Percent",
                "Allele4_Seq_Count_Percent",
                "Allele5_Seq_Count_Percent"
            ]

            df = pd.DataFrame(results, columns=columns)

            output_dir = os.path.dirname(args.output_xlsx)
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)

            df.to_excel(args.output_xlsx, index=False)
            logger.info(f"Results written to {args.output_xlsx}")

    except Exception as e:
        logger.error(f"Variant calling failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()