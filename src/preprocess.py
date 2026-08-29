#!/usr/bin/env python3
"""
Extract upstream sequence from reference based on target.
"""
import argparse
import logging
import os
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq


ALLOWED_BASES = set("ACGTN")


def setup_logger(name, log_dir="./logs", level=logging.INFO):
    """Create a logger with console and file handlers."""
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


def validate_sequence(sequence, name):
    """Validate a nucleotide sequence."""
    if not sequence:
        raise ValueError(f"{name} is empty.")

    invalid = set(sequence.upper()) - ALLOWED_BASES
    if invalid:
        raise ValueError(
            f"{name} contains invalid characters: "
            f"{''.join(sorted(invalid))}. "
            f"Only A/C/G/T/N are supported."
        )


def read_reference_records(filepath, logger):
    """
    Read reference sequences.

    Supports:
    1. FASTA files with one or more records.
    2. Plain sequence text without a FASTA header.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Reference file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8") as handle:
        content = handle.read().strip()

    if not content:
        raise ValueError(f"Reference file is empty: {filepath}")

    records = []

    if content.startswith(">"):
        for record in SeqIO.parse(filepath, "fasta"):
            sequence = str(record.seq).replace(" ", "").upper()
            validate_sequence(sequence, f"Reference record '{record.id}'")
            records.append((record.id, sequence))
    else:
        sequence = "".join(
            line.strip()
            for line in content.splitlines()
            if line.strip() and not line.startswith(">")
        )
        sequence = sequence.replace(" ", "").upper()

        record_id = Path(filepath).stem
        validate_sequence(sequence, f"Reference record '{record_id}'")
        records.append((record_id, sequence))

    if not records:
        raise ValueError(f"No valid reference records found in {filepath}")

    logger.info(
        "Loaded %d reference record(s): %s",
        len(records),
        ", ".join(f"{record_id}({len(seq)} bp)" for record_id, seq in records)
    )

    return records


def read_target_sequence(filepath, logger):
    """
    Read the first target sequence.

    Supports FASTA and plain sequence text.
    """
    if not os.path.isfile(filepath):
        raise FileNotFoundError(f"Target file not found: {filepath}")

    with open(filepath, "r", encoding="utf-8") as handle:
        content = handle.read().strip()

    if not content:
        raise ValueError(f"Target file is empty: {filepath}")

    if content.startswith(">"):
        records = list(SeqIO.parse(filepath, "fasta"))
        if not records:
            raise ValueError(f"No FASTA record found in {filepath}")

        target_id = records[0].id
        target_seq = str(records[0].seq).replace(" ", "").upper()
    else:
        target_id = Path(filepath).stem
        target_seq = "".join(
            line.strip()
            for line in content.splitlines()
            if line.strip() and not line.startswith(">")
        )
        target_seq = target_seq.replace(" ", "").upper()

    validate_sequence(target_seq, "Target sequence")

    logger.info(
        "Loaded target '%s', length=%d bp: %s",
        target_id,
        len(target_seq),
        target_seq
    )

    return target_id, target_seq


def find_all_occurrences(sequence, query):
    """Return all overlapping occurrences of query in sequence."""
    positions = []
    search_start = 0

    while True:
        position = sequence.find(query, search_start)
        if position == -1:
            break

        positions.append(position)
        search_start = position + 1

    return positions


def locate_target(reference_records, target_seq, logger):
    """
    Locate target or its reverse complement in all reference records.

    Returns a dictionary containing:
        chrom
        reference_seq
        target_start
        target_end
        orientation
        matched_sequence
    """
    target_seq = target_seq.upper()
    target_rc = str(Seq(target_seq).reverse_complement()).upper()

    matches = []

    for chrom, reference_seq in reference_records:
        for position in find_all_occurrences(reference_seq, target_seq):
            matches.append({
                "chrom": chrom,
                "reference_seq": reference_seq,
                "target_start": position,
                "target_end": position + len(target_seq),
                "orientation": "forward",
                "matched_sequence": target_seq
            })

        # 避免回文序列被重复记录
        if target_rc != target_seq:
            for position in find_all_occurrences(reference_seq, target_rc):
                matches.append({
                    "chrom": chrom,
                    "reference_seq": reference_seq,
                    "target_start": position,
                    "target_end": position + len(target_rc),
                    "orientation": "reverse",
                    "matched_sequence": target_rc
                })

    if not matches:
        raise ValueError(
            "Target sequence was not found in the reference, either in "
            "forward orientation or as a reverse complement."
        )

    if len(matches) > 1:
        match_text = "; ".join(
            (
                f"{item['chrom']}:{item['target_start'] + 1}-"
                f"{item['target_end']}({item['orientation']})"
            )
            for item in matches
        )

        raise ValueError(
            f"Target sequence has {len(matches)} matches in the reference: "
            f"{match_text}. A unique target sequence is required."
        )

    match = matches[0]

    logger.info(
        "Target located at %s:%d-%d, orientation=%s",
        match["chrom"],
        match["target_start"] + 1,
        match["target_end"],
        match["orientation"]
    )

    return match


def extract_upstream_anchor(match, upstream_len, logger):
    """
    Extract an upstream anchor in the biological orientation of the target.
    """
    if upstream_len <= 0:
        raise ValueError("--upstream must be a positive integer.")

    reference_seq = match["reference_seq"]
    target_start = match["target_start"]
    target_end = match["target_end"]
    orientation = match["orientation"]

    if orientation == "forward":
        anchor_start = target_start - upstream_len
        anchor_end = target_start

        if anchor_start < 0:
            raise ValueError(
                "Not enough upstream sequence before the forward target: "
                f"target_start={target_start}, requested={upstream_len}."
            )

        anchor_seq = reference_seq[anchor_start:anchor_end]

    elif orientation == "reverse":
        # 对于反向靶点，其生物学上游位于参考序列坐标的 target_end 之后
        anchor_start = target_end
        anchor_end = target_end + upstream_len

        if anchor_end > len(reference_seq):
            raise ValueError(
                "Not enough upstream sequence for the reverse target: "
                f"target_end={target_end}, requested={upstream_len}, "
                f"reference_length={len(reference_seq)}."
            )

        reference_fragment = reference_seq[anchor_start:anchor_end]
        anchor_seq = str(Seq(reference_fragment).reverse_complement()).upper()

    else:
        raise ValueError(f"Unsupported target orientation: {orientation}")

    logger.info(
        "Extracted upstream anchor: genomic interval %s:%d-%d, "
        "target orientation=%s, anchor=%s",
        match["chrom"],
        anchor_start + 1,
        anchor_end,
        orientation,
        anchor_seq
    )

    return anchor_seq, anchor_start, anchor_end


def write_anchor_fasta(
    output_fasta,
    anchor_seq,
    match,
    anchor_start,
    anchor_end
):
    """
    Write anchor FASTA with target locus metadata.

    Coordinates in the metadata are 0-based half-open coordinates.
    """
    output_dir = os.path.dirname(output_fasta)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)

    header = (
        f"upstream_seq "
        f"chrom={match['chrom']} "
        f"target_start={match['target_start']} "
        f"target_end={match['target_end']} "
        f"orientation={match['orientation']} "
        f"anchor_start={anchor_start} "
        f"anchor_end={anchor_end}"
    )

    with open(output_fasta, "w", encoding="utf-8") as handle:
        handle.write(f">{header}\n")
        handle.write(f"{anchor_seq}\n")


def extract_upstream(
    reference_fasta,
    target_fasta,
    output_fasta,
    upstream_len,
    logger
):
    """Main preprocessing procedure."""
    reference_records = read_reference_records(reference_fasta, logger)
    _, target_seq = read_target_sequence(target_fasta, logger)

    match = locate_target(reference_records, target_seq, logger)

    anchor_seq, anchor_start, anchor_end = extract_upstream_anchor(
        match=match,
        upstream_len=upstream_len,
        logger=logger
    )

    write_anchor_fasta(
        output_fasta=output_fasta,
        anchor_seq=anchor_seq,
        match=match,
        anchor_start=anchor_start,
        anchor_end=anchor_end
    )

    logger.info("Upstream anchor written to %s", output_fasta)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Locate a target in a reference and extract its upstream anchor."
        )
    )

    parser.add_argument(
        "--reference",
        required=True,
        help="Reference FASTA or plain sequence file"
    )
    parser.add_argument(
        "--target",
        required=True,
        help="Target FASTA or plain sequence file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output upstream-anchor FASTA"
    )
    parser.add_argument(
        "--upstream",
        type=int,
        default=10,
        help="Number of upstream bases to extract; default: 10"
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level"
    )

    args = parser.parse_args()

    logger = setup_logger(
        "preprocess",
        level=getattr(logging, args.log_level.upper())
    )

    try:
        extract_upstream(
            reference_fasta=args.reference,
            target_fasta=args.target,
            output_fasta=args.output,
            upstream_len=args.upstream,
            logger=logger
        )
        logger.info("Preprocessing completed successfully.")
    except Exception as error:
        logger.exception("Preprocessing failed: %s", error)
        sys.exit(1)


if __name__ == "__main__":
    main()