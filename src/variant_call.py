#!/usr/bin/env python3
"""
Call variants using samtools mpileup, extract alleles after target sequence.
"""
import argparse
import logging
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def setup_logger(name, log_dir="./logs", level=logging.INFO):
    """Create a logger with console and file output."""
    os.makedirs(log_dir, exist_ok=True)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.propagate = False

    if logger.handlers:
        logger.handlers.clear()

    log_file = os.path.join(log_dir, f"{name}.log")

    file_handler = logging.FileHandler(log_file, encoding="utf-8")
    file_handler.setLevel(level)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(level)

    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    return logger


def validate_filter_threshold(value):
    """Validate --f as an integer from 0 to 20."""
    try:
        integer_value = int(value)
    except ValueError as error:
        raise argparse.ArgumentTypeError(
            "--f must be an integer."
        ) from error

    if integer_value < 0 or integer_value > 20:
        raise argparse.ArgumentTypeError(
            "--f must be an integer between 0 and 20."
        )

    return integer_value


def read_reference(reference_fasta):
    """Read all reference FASTA records into a dictionary."""
    if not os.path.isfile(reference_fasta):
        raise FileNotFoundError(
            f"Reference FASTA not found: {reference_fasta}"
        )

    reference_dict = {}

    for record in SeqIO.parse(reference_fasta, "fasta"):
        reference_dict[record.id] = str(record.seq).upper()

    if not reference_dict:
        raise ValueError(
            f"No FASTA records found in reference: {reference_fasta}"
        )

    return reference_dict


def parse_fasta_metadata(description):
    """
    Parse key=value metadata from a FASTA description.

    Example:
        upstream_seq chrom=SH2 target_start=100 target_end=120
        orientation=forward anchor_start=90 anchor_end=100
    """
    metadata = {}

    for token in description.split()[1:]:
        if "=" not in token:
            continue

        key, value = token.split("=", 1)
        metadata[key] = value

    return metadata


def read_anchor_fasta(anchor_fasta):
    """
    Read the first anchor FASTA record.

    Returns
    -------
    anchor_sequence : str
    metadata : dict
    """
    if not os.path.isfile(anchor_fasta):
        raise FileNotFoundError(
            f"Anchor FASTA not found: {anchor_fasta}"
        )

    records = list(SeqIO.parse(anchor_fasta, "fasta"))

    if not records:
        raise ValueError(
            f"No FASTA record found in anchor file: {anchor_fasta}"
        )

    record = records[0]
    anchor_sequence = str(record.seq).upper()
    metadata = parse_fasta_metadata(record.description)

    if not anchor_sequence:
        raise ValueError("Upstream anchor sequence is empty.")

    return anchor_sequence, metadata


def find_all_occurrences(sequence, query):
    """Return all overlapping occurrences of a query."""
    positions = []
    search_start = 0

    while True:
        position = sequence.find(query, search_start)
        if position == -1:
            break

        positions.append(position)
        search_start = position + 1

    return positions


def determine_target_interval_from_metadata(
    reference_dict,
    metadata,
    total_len,
    logger
):
    """
    Determine the genomic target interval using preprocess metadata.

    Returns
    -------
    chrom : str
    interval_start : int
        0-based inclusive.
    interval_end : int
        0-based exclusive.
    orientation : str
        forward or reverse.
    """
    required_keys = {
        "chrom",
        "target_start",
        "target_end",
        "orientation"
    }

    if not required_keys.issubset(metadata):
        return None

    chrom = metadata["chrom"]

    if chrom not in reference_dict:
        raise ValueError(
            f"Chromosome/reference record '{chrom}' from anchor metadata "
            f"was not found in the reference FASTA."
        )

    try:
        target_start = int(metadata["target_start"])
        target_end = int(metadata["target_end"])
    except ValueError as error:
        raise ValueError(
            "Invalid target coordinates in anchor FASTA metadata."
        ) from error

    orientation = metadata["orientation"].lower()
    reference_length = len(reference_dict[chrom])

    target_len = target_end - target_start

    if target_len <= 0:
        raise ValueError(
            "Invalid target coordinates: target_end must be greater "
            f"than target_start (got {target_start}-{target_end})."
        )

    # If no total length is requested, report the exact target (no flanking).
    if total_len is None:
        total_len = target_len

    if total_len < target_len:
        raise ValueError(
            f"--l ({total_len}) is shorter than the target length "
            f"({target_len} bp). The total region must be at least as "
            f"long as the target."
        )

    # Symmetric flanking around the target.
    flank_total = total_len - target_len
    flank_left = flank_total // 2
    flank_right = flank_total - flank_left

    if orientation == "forward":
        interval_start = target_start - flank_left
        interval_end = target_end + flank_right

    elif orientation in {"reverse", "reverse_complement"}:
        # For a reverse target the biological 5' end maps to the reference
        # target_end, so the left/right flanking sides are swapped.
        interval_start = target_start - flank_right
        interval_end = target_end + flank_left
        orientation = "reverse"

    else:
        raise ValueError(
            f"Invalid target orientation in metadata: {orientation}"
        )

    if interval_start < 0 or interval_end > reference_length:
        raise ValueError(
            "Requested allele interval extends beyond the reference: "
            f"{chrom}:{interval_start}-{interval_end}, "
            f"reference length={reference_length}. "
            "Try a smaller --l."
        )

    logger.info(
        "Using target metadata: chrom=%s, interval=%d-%d, "
        "target_len=%d, total_len=%d, orientation=%s",
        chrom,
        interval_start + 1,
        interval_end,
        target_len,
        total_len,
        orientation
    )

    return chrom, interval_start, interval_end, orientation


def determine_target_interval_from_anchor(
    reference_dict,
    anchor_sequence,
    total_len,
    logger
):
    """
    Fallback for old index.fasta files without metadata.

    Locate the anchor sequence or its reverse complement in the reference.
    """
    anchor_sequence = anchor_sequence.upper()
    anchor_rc = str(
        Seq(anchor_sequence).reverse_complement()
    ).upper()

    matches = []

    for chrom, reference_sequence in reference_dict.items():
        for position in find_all_occurrences(
            reference_sequence,
            anchor_sequence
        ):
            matches.append({
                "chrom": chrom,
                "position": position,
                "orientation": "forward"
            })

        if anchor_rc != anchor_sequence:
            for position in find_all_occurrences(
                reference_sequence,
                anchor_rc
            ):
                matches.append({
                    "chrom": chrom,
                    "position": position,
                    "orientation": "reverse"
                })

    if not matches:
        raise ValueError(
            "The upstream anchor was not found in the reference."
        )

    if len(matches) > 1:
        match_text = "; ".join(
            f"{item['chrom']}:{item['position'] + 1}"
            f"({item['orientation']})"
            for item in matches
        )
        raise ValueError(
            f"The short upstream anchor has {len(matches)} reference "
            f"matches: {match_text}. Regenerate index.fasta using the "
            f"corrected preprocess.py so that locus metadata is included."
        )

    match = matches[0]
    chrom = match["chrom"]
    position = match["position"]
    orientation = match["orientation"]

    # Legacy fallback: without target metadata we cannot center the target,
    # so the requested length is taken directly downstream of the anchor.
    if total_len is None:
        total_len = 20

    if orientation == "forward":
        interval_start = position + len(anchor_sequence)
        interval_end = interval_start + total_len
    else:
        interval_end = position
        interval_start = interval_end - total_len

    if (
        interval_start < 0
        or interval_end > len(reference_dict[chrom])
    ):
        raise ValueError(
            "Target interval inferred from anchor extends beyond reference."
        )

    logger.warning(
        "Anchor FASTA does not contain target metadata. "
        "Using anchor sequence fallback."
    )
    logger.info(
        "Anchor-based target interval: %s:%d-%d, orientation=%s",
        chrom,
        interval_start + 1,
        interval_end,
        orientation
    )

    return chrom, interval_start, interval_end, orientation


def determine_target_interval(
    reference_dict,
    anchor_sequence,
    metadata,
    total_len,
    logger
):
    """Determine target interval from metadata or anchor sequence."""
    metadata_result = determine_target_interval_from_metadata(
        reference_dict=reference_dict,
        metadata=metadata,
        total_len=total_len,
        logger=logger
    )

    if metadata_result is not None:
        return metadata_result

    return determine_target_interval_from_anchor(
        reference_dict=reference_dict,
        anchor_sequence=anchor_sequence,
        total_len=total_len,
        logger=logger
    )


def get_wildtype_allele(
    reference_sequence,
    interval_start,
    interval_end,
    orientation
):
    """Extract expected wild-type allele in biological target orientation."""
    wildtype = reference_sequence[interval_start:interval_end].upper()

    if orientation == "reverse":
        wildtype = str(
            Seq(wildtype).reverse_complement()
        ).upper()

    return wildtype


def ensure_bam_index(bam_path, threads, logger):
    """Create BAM index if it does not already exist."""
    possible_indexes = [
        f"{bam_path}.bai",
        str(Path(bam_path).with_suffix(".bai"))
    ]

    if any(os.path.isfile(index) for index in possible_indexes):
        return True

    logger.info("BAM index not found; indexing %s", bam_path)

    command = [
        "samtools",
        "index",
        "-@",
        str(max(1, threads)),
        bam_path
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        logger.error(
            "Failed to index BAM %s: %s",
            bam_path,
            result.stderr.strip()
        )
        return False

    logger.info("BAM index generated successfully.")
    return True


def parse_cigar(cigar):
    """Parse a CIGAR string."""
    if not cigar or cigar == "*":
        return []

    operations = [
        (int(length), operation)
        for length, operation in CIGAR_PATTERN.findall(cigar)
    ]

    reconstructed = "".join(
        f"{length}{operation}"
        for length, operation in operations
    )

    if reconstructed != cigar:
        raise ValueError(f"Unsupported or invalid CIGAR string: {cigar}")

    return operations


def get_reference_end(alignment_start, cigar_operations):
    """Calculate 0-based exclusive reference end from CIGAR."""
    reference_position = alignment_start

    for length, operation in cigar_operations:
        if operation in {"M", "D", "N", "=", "X"}:
            reference_position += length

    return reference_position


def extract_allele_from_alignment(
    query_sequence,
    cigar,
    alignment_start,
    interval_start,
    interval_end
):
    """
    Reconstruct an allele for a reference interval using CIGAR.

    Insertions inside the interval are included.
    Deleted reference bases are omitted.

    Returns
    -------
    str or None
        None means the read does not fully span the target interval.
        An empty reconstructed allele is represented as "-".
    """
    query_sequence = query_sequence.upper()
    cigar_operations = parse_cigar(cigar)

    if not cigar_operations:
        return None

    alignment_end = get_reference_end(
        alignment_start,
        cigar_operations
    )

    # 只统计完整跨越整个目标区域的 reads，避免部分 reads 产生截短 allele
    if (
        alignment_start > interval_start
        or alignment_end < interval_end
    ):
        return None

    query_position = 0
    reference_position = alignment_start
    allele_parts = []

    for length, operation in cigar_operations:
        if operation in {"M", "=", "X"}:
            operation_ref_start = reference_position
            operation_ref_end = reference_position + length

            overlap_start = max(
                operation_ref_start,
                interval_start
            )
            overlap_end = min(
                operation_ref_end,
                interval_end
            )

            if overlap_start < overlap_end:
                query_offset_start = (
                    query_position
                    + overlap_start
                    - operation_ref_start
                )
                query_offset_end = (
                    query_position
                    + overlap_end
                    - operation_ref_start
                )

                allele_parts.append(
                    query_sequence[
                        query_offset_start:query_offset_end
                    ]
                )

            query_position += length
            reference_position += length

        elif operation == "I":
            # 插入发生在当前 reference_position 之前。
            # 只保留位于目标区间内部的插入。
            if interval_start <= reference_position < interval_end:
                allele_parts.append(
                    query_sequence[
                        query_position:query_position + length
                    ]
                )

            query_position += length

        elif operation in {"D", "N"}:
            # 缺失/跳跃消耗参考序列，但不消耗 query。
            reference_position += length

        elif operation == "S":
            query_position += length

        elif operation in {"H", "P"}:
            # Hard clip/padding 不消耗 query_sequence。
            continue

        else:
            raise ValueError(
                f"Unsupported CIGAR operation: {operation}"
            )

    allele = "".join(allele_parts).upper()

    if not allele:
        return "-"

    return allele


def fetch_alleles_from_bam(
    bam_path,
    chrom,
    interval_start,
    interval_end,
    orientation,
    threads=1,
    mapq=0,
    logger=None
):
    """
    Fetch alignments and reconstruct target alleles.

    SAM sequences are represented in alignment/reference orientation.
    For a reverse target, the reconstructed allele is reverse-complemented
    into the biological orientation of the input target.
    """
    region = f"{chrom}:{interval_start + 1}-{interval_end}"

    command = [
        "samtools",
        "view",
        "-@",
        str(max(1, threads)),
        "-F",
        "0x904",
        "-q",
        str(max(0, mapq)),
        bam_path,
        region
    ]

    result = subprocess.run(
        command,
        capture_output=True,
        text=True
    )

    if result.returncode != 0:
        if logger:
            logger.error(
                "Failed to read BAM %s: %s",
                bam_path,
                result.stderr.strip()
            )
        return []

    alleles = []

    for line_number, line in enumerate(
        result.stdout.splitlines(),
        start=1
    ):
        if not line:
            continue

        fields = line.split("\t")

        if len(fields) < 11:
            continue

        try:
            alignment_start = int(fields[3]) - 1
            cigar = fields[5]
            query_sequence = fields[9]

            if query_sequence == "*" or cigar == "*":
                continue

            allele = extract_allele_from_alignment(
                query_sequence=query_sequence,
                cigar=cigar,
                alignment_start=alignment_start,
                interval_start=interval_start,
                interval_end=interval_end
            )

            if allele is None:
                continue

            if orientation == "reverse":
                if allele == "-":
                    oriented_allele = "-"
                else:
                    oriented_allele = str(
                        Seq(allele).reverse_complement()
                    ).upper()
            else:
                oriented_allele = allele

            alleles.append(oriented_allele)

        except Exception as error:
            if logger:
                logger.debug(
                    "Skipping malformed alignment line %d from %s: %s",
                    line_number,
                    bam_path,
                    error
                )

    return alleles


def get_mutation_category(allele_seq, wt_seq):
    """Classify an allele relative to wild type."""
    if allele_seq == wt_seq:
        return "WT"

    if allele_seq == "-":
        return "Deletion"

    if len(allele_seq) > len(wt_seq):
        return "Insertion"

    if len(allele_seq) < len(wt_seq):
        return "Deletion"

    mismatch_count = sum(
        allele_base != wt_base
        for allele_base, wt_base in zip(allele_seq, wt_seq)
    )

    if mismatch_count == 1:
        return "Substitution"

    return "Complex"


def apply_allele_filter(
    sorted_items,
    total_reads,
    filter_threshold
):
    """
    Retain alleles whose raw frequency is at least filter_threshold.
    """
    filtered_items = []

    for sequence, count in sorted_items:
        raw_percent = (
            count / total_reads * 100.0
            if total_reads > 0
            else 0.0
        )

        if raw_percent >= filter_threshold:
            filtered_items.append({
                "seq": sequence,
                "count": count,
                "raw_pct": raw_percent
            })

    return filtered_items


def recalculate_percentages(filtered_items):
    """Recalculate frequencies using retained reads only."""
    filtered_total = sum(
        item["count"]
        for item in filtered_items
    )

    recalculated = []

    for item in filtered_items:
        filtered_percent = (
            item["count"] / filtered_total * 100.0
            if filtered_total > 0
            else 0.0
        )

        recalculated.append({
            "seq": item["seq"],
            "count": item["count"],
            "raw_pct": item["raw_pct"],
            "filtered_pct": filtered_percent
        })

    return recalculated, filtered_total


def infer_genotype_from_filtered_alleles(
    filtered_items,
    wt_allele
):
    """
    Infer a simple genotype from filtered alleles.

    Rules
    -----
    1 allele:
        WT allele       -> Wild Type
        mutant allele   -> Homozygous Mutant

    2 alleles:
        WT + mutant with minor allele >10% -> Heterozygous Mutant
        two different mutants >10%         -> Biallelic Mutant
        minor allele <=10%                  -> possible chimera/contamination

    >=3 alleles:
        Chimeric
    """
    if not filtered_items:
        return "No_Alleles_After_Filter"

    if len(filtered_items) >= 3:
        return "Chimeric"

    first = filtered_items[0]
    first_is_wt = first["seq"] == wt_allele

    if len(filtered_items) == 1:
        if first_is_wt:
            return "Wild Type"
        return "Homozygous Mutant"

    second = filtered_items[1]
    second_is_wt = second["seq"] == wt_allele
    minor_percent = second["filtered_pct"]

    if minor_percent <= 10.0:
        if first_is_wt:
            return "WT with Low-Frequency Variant"
        return "Homozygous Mutant with Low-Frequency Variant"

    if first_is_wt or second_is_wt:
        return "Heterozygous Mutant"

    # 两个不同的非 WT allele 应判定为双等位突变，
    # 不能仅因二者同为 deletion/insertion 就判为纯合。
    return "Biallelic Mutant"


def format_allele_cell(sequence, count, percent, wt_allele):
    """Format an allele cell for Excel."""
    category = get_mutation_category(sequence, wt_allele)

    return (
        f"{sequence}: {count} reads, "
        f"{percent:.2f}%, {category}"
    )


def write_results_excel(results, output_xlsx, res=5):
    """Write final results to Excel."""
    columns = [
        "Sample",
        "Genotype",
        "WildType_Allele",
        "Target_Orientation",
        "Target_Region",
        "Total_Target_Reads",
        "Filter_Threshold_Percent",
        "Reads_Used_After_Filter",
        "Reads_Used_Percent_Of_Total",
    ]

    columns += [
        f"Allele{i}_Seq_Count_Percent_Type"
        for i in range(1, res + 1)
    ]

    dataframe = pd.DataFrame(results, columns=columns)

    output_dir = os.path.dirname(output_xlsx)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    dataframe.to_excel(output_xlsx, index=False)


def write_top_reads_fasta(entries, output_fasta):
    """Write top read sequences to a FASTA file.

    Parameters
    ----------
    entries : list of tuples
        (sample, rank, sequence, count, percent, category)
    """
    if not entries:
        return

    output_dir = os.path.dirname(output_fasta)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    with open(output_fasta, "w", encoding="utf-8") as handle:
        for sample, rank, sequence, count, percent, category in entries:
            handle.write(
                f">sample={sample} rank={rank} count={count} "
                f"percent={percent:.2f}% type={category}\n"
            )
            handle.write(f"{sequence}\n")


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct target alleles from BAM CIGAR alignments."
        )
    )

    parser.add_argument(
        "--bam_dir",
        required=True,
        help="Directory containing *.sorted.bam files"
    )
    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA"
    )
    parser.add_argument(
        "--target_fasta",
        "--anchor_fasta",
        dest="anchor_fasta",
        required=True,
        help=(
            "Upstream anchor FASTA generated by preprocess.py. "
            "--target_fasta is retained for compatibility."
        )
    )
    parser.add_argument(
        "--output_xlsx",
        required=True,
        help="Output Excel file"
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of samtools threads"
    )
    parser.add_argument(
        "--qual",
        type=int,
        default=20,
        help=(
            "Compatibility parameter inherited from the original pipeline. "
            "FASTQ quality filtering should be performed in qc_filter.py."
        )
    )
    parser.add_argument(
        "--mapq",
        type=int,
        default=0,
        help="Minimum BAM mapping quality; default: 0"
    )
    parser.add_argument(
        "--l",
        "--extract_len",
        dest="total_len",
        type=int,
        default=None,
        help=(
            "Total length of the reported target region (bp), centered on "
            "the target with symmetric left/right flanking. For example, "
            "with a 20 bp target and --l 60, 20 bp of flanking is added to "
            "each side (60 bp total). If omitted or shorter than the target, "
            "the exact target length is used. "
            "--extract_len is retained as an alias."
        )
    )
    parser.add_argument(
        "--res",
        type=int,
        default=5,
        help=(
            "Number of top allele/read sequences to export, ordered by read "
            "count. If fewer distinct alleles exist, all are exported. "
            "default: 5"
        )
    )
    parser.add_argument(
        "--f",
        type=validate_filter_threshold,
        default=0,
        help=(
            "Allele-frequency filter threshold in percent, 0-20; "
            "default: 0"
        )
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )

    args = parser.parse_args()

    logger = setup_logger(
        "variant_call",
        level=getattr(logging, args.log_level.upper())
    )

    try:
        if args.threads < 1:
            raise ValueError("--threads must be a positive integer.")

        if args.total_len is not None and args.total_len < 1:
            raise ValueError("--l must be a positive integer.")

        if args.res < 1:
            raise ValueError("--res must be a positive integer.")

        if args.mapq < 0:
            raise ValueError("--mapq must be a non-negative integer.")

        reference_dict = read_reference(args.reference)

        anchor_sequence, metadata = read_anchor_fasta(
            args.anchor_fasta
        )

        (
            chrom,
            interval_start,
            interval_end,
            orientation
        ) = determine_target_interval(
            reference_dict=reference_dict,
            anchor_sequence=anchor_sequence,
            metadata=metadata,
            total_len=args.total_len,
            logger=logger
        )

        wildtype_allele = get_wildtype_allele(
            reference_sequence=reference_dict[chrom],
            interval_start=interval_start,
            interval_end=interval_end,
            orientation=orientation
        )

        target_region_text = (
            f"{chrom}:{interval_start + 1}-{interval_end}"
        )

        logger.info("Upstream anchor: %s", anchor_sequence)
        logger.info("Target region: %s", target_region_text)
        logger.info("Target orientation: %s", orientation)
        logger.info(
            "Wild-type allele (%d bp): %s",
            len(wildtype_allele),
            wildtype_allele
        )
        logger.info("Allele filter threshold: %d%%", args.f)
        logger.info("Minimum mapping quality: %d", args.mapq)
        logger.info("Reported target length (--l): %s", args.total_len)
        logger.info("Top reads to export (--res): %d", args.res)
        logger.info(
            "--qual=%d is used by the earlier FASTQ QC step; "
            "variant_call.py does not reinterpret it as base quality.",
            args.qual
        )

        bam_files = sorted(
            Path(args.bam_dir).glob("*.sorted.bam")
        )

        if not bam_files:
            raise FileNotFoundError(
                f"No *.sorted.bam files found in {args.bam_dir}"
            )

        results = []
        top_reads_entries = []

        for bam_path in bam_files:
            sample = bam_path.name

            if sample.endswith(".sorted.bam"):
                sample = sample[:-len(".sorted.bam")]

            logger.info("Processing sample: %s", sample)

            if not ensure_bam_index(
                str(bam_path),
                args.threads,
                logger
            ):
                logger.error(
                    "Skipping sample %s because BAM indexing failed.",
                    sample
                )
                continue

            alleles = fetch_alleles_from_bam(
                bam_path=str(bam_path),
                chrom=chrom,
                interval_start=interval_start,
                interval_end=interval_end,
                orientation=orientation,
                threads=args.threads,
                mapq=args.mapq,
                logger=logger
            )

            if not alleles:
                logger.warning(
                    "No reads fully spanning the target were found for %s",
                    sample
                )

                total_target_reads = 0
                filtered_target_reads = 0
                filtered_percent_of_total = 0.0
                genotype = "NoTargetReads"
                top_display = []

            else:
                allele_counter = Counter(alleles)
                total_target_reads = len(alleles)

                sorted_items = sorted(
                    allele_counter.items(),
                    key=lambda item: (-item[1], item[0])
                )

                filtered_items = apply_allele_filter(
                    sorted_items=sorted_items,
                    total_reads=total_target_reads,
                    filter_threshold=args.f
                )

                (
                    recalculated_items,
                    filtered_target_reads
                ) = recalculate_percentages(filtered_items)

                genotype = infer_genotype_from_filtered_alleles(
                    filtered_items=recalculated_items,
                    wt_allele=wildtype_allele
                )

                filtered_percent_of_total = (
                    filtered_target_reads
                    / total_target_reads
                    * 100.0
                    if total_target_reads > 0
                    else 0.0
                )

                top_display = recalculated_items[:args.res]

                for rank, (sequence, count) in enumerate(
                    sorted_items[:args.res],
                    start=1
                ):
                    top_reads_entries.append((
                        sample,
                        rank,
                        sequence,
                        count,
                        count / total_target_reads * 100.0,
                        get_mutation_category(sequence, wildtype_allele)
                    ))

                logger.info(
                    "Sample=%s, total spanning reads=%d, "
                    "filtered reads=%d, genotype=%s",
                    sample,
                    total_target_reads,
                    filtered_target_reads,
                    genotype
                )

                logger.info(
                    "Raw top alleles: %s",
                    [
                        {
                            "sequence": sequence,
                            "count": count,
                            "percent": (
                                count / total_target_reads * 100.0
                            ),
                            "type": get_mutation_category(
                                sequence,
                                wildtype_allele
                            )
                        }
                        for sequence, count in sorted_items[:10]
                    ]
                )

            row = [
                sample,
                genotype,
                wildtype_allele,
                orientation,
                target_region_text,
                total_target_reads,
                args.f,
                filtered_target_reads,
                f"{filtered_percent_of_total:.2f}%"
            ]

            for item in top_display:
                row.append(
                    format_allele_cell(
                        sequence=item["seq"],
                        count=item["count"],
                        percent=item["filtered_pct"],
                        wt_allele=wildtype_allele
                    )
                )

            while len(row) < 9 + args.res:
                row.append("")

            results.append(row)

        if not results:
            raise RuntimeError(
                "No BAM samples were processed successfully."
            )

        write_results_excel(
            results=results,
            output_xlsx=args.output_xlsx,
            res=args.res
        )

        top_reads_fasta = (
            os.path.splitext(args.output_xlsx)[0] + "_top_reads.fa"
        )
        write_top_reads_fasta(top_reads_entries, top_reads_fasta)

        logger.info(
            "Results written successfully to %s",
            args.output_xlsx
        )
        logger.info(
            "Top reads (%d per sample) written to %s",
            args.res,
            top_reads_fasta
        )

    except Exception as error:
        logger.exception("Variant calling failed: %s", error)
        sys.exit(1)


if __name__ == "__main__":
    main()