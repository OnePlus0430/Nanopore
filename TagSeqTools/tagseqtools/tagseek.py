#!/usr/bin/env python3
"""
tagseek.py - Tag Sequence Identification Module

This module provides functionality to identify and separate tagged (NAD-capped)
and non-tagged RNA reads from Nanopore sequencing FASTQ files. It searches for
synthetic tag sequences in the first N bases of each read.

The module uses a pattern-based matching approach where multiple overlapping
patterns are generated from the tag sequence based on the similarity threshold.

Functions:
    - generate_patterns: Generate search patterns from a tag sequence
    - run_tagseek: Main function to process FASTQ and separate tagged/non-tagged reads

Output Files:
    - {prefix}.tag.fastq: Reads containing the tag sequence
    - {prefix}.nontag.fastq: Reads without the tag sequence
    - Tag_statistics.txt: Summary statistics of the tagging process

Command-line usage:
    $ tagseqtools tagseek --fastq sample --tag 'CCTGAA...' --similarity 12

Author: Huan ZHONG
License: Apache 2.0
"""

from __future__ import annotations

from collections import Counter
from pathlib import Path
from typing import Optional, TextIO, Dict, Tuple

try:
    import regex
except ImportError:
    raise ImportError(
        "Required package 'regex' is not installed. "
        "Please install it with: pip install regex"
    )

from tagseqtools.utils import log_msg, ensure_exists


def generate_patterns(tag_sequence: str, similarity: int) -> Tuple[str, ...]:
    """
    Generate search patterns from a tag sequence based on similarity threshold.

    Creates overlapping subsequences of the tag sequence that will be used
    for pattern matching. RNA sequences (containing U) are automatically
    converted to DNA (U -> T).

    Args:
        tag_sequence: The synthetic tag RNA/DNA sequence to search for.
                     Can contain either U (RNA) or T (DNA).
        similarity: The length of each pattern (number of consecutive bases).
                   Also determines the number of patterns generated.

    Returns:
        A tuple of unique pattern strings to search for.

    Example:
        >>> patterns = generate_patterns("CCTGAACCTGAA", similarity=6)
        >>> print(patterns)
        ('CCTGAA', 'CTGAAC', 'TGAACC', 'GAACCT', 'AACCTG', 'ACCTGA', 'CCTGAA')

    Note:
        The number of patterns generated is approximately (similarity + 1),
        minus any duplicates that may occur in repetitive sequences.
    """
    # Convert RNA to DNA (U -> T)
    tag_dna = tag_sequence.upper().replace('U', 'T')
    
    # Generate overlapping patterns
    patterns = set()
    for i in range(similarity + 1):
        if i + similarity <= len(tag_dna):
            patterns.add(tag_dna[i:similarity + i])
    
    return tuple(patterns)


def run_tagseek(
    input_fastq_prefix: str,
    tag_sequence: str,
    similarity: int,
    search_window: int = 40,
    output_dir: Path = Path("."),
    logfp: Optional[TextIO] = None
) -> Dict:
    """
    Identify and separate tagged and non-tagged reads from a FASTQ file.

    This function reads a FASTQ file, searches for tag patterns in the first
    N bases of each read, and separates reads into tagged and non-tagged
    output files. It handles corrupted FASTQ records gracefully by skipping
    malformed entries.

    Args:
        input_fastq_prefix: Prefix of the input FASTQ file (without .fastq extension).
                           The function expects a file named "{prefix}.fastq".
        tag_sequence: The synthetic tag sequence to search for (RNA or DNA format).
        similarity: Number of consecutive matching bases required for tag detection.
        search_window: Number of bases from the start of each read to search.
                      Default is 40 (first 40 bp).
        output_dir: Directory for output files. Default is current directory.
        logfp: Optional file handle for logging. If None, only prints to console.

    Returns:
        Dictionary containing:
            - 'total': Total number of valid reads processed
            - 'tagged': Number of reads identified as tagged
            - 'skipped': Number of corrupted reads skipped
            - 'tag_fastq': Path to the tagged reads output file
            - 'nontag_fastq': Path to the non-tagged reads output file
            - 'stats_file': Path to the statistics file

    Raises:
        FileNotFoundError: If the input FASTQ file does not exist.

    Example:
        >>> from pathlib import Path
        >>> result = run_tagseek(
        ...     input_fastq_prefix="sample",
        ...     tag_sequence="CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT",
        ...     similarity=12,
        ...     search_window=40,
        ...     output_dir=Path("results")
        ... )
        >>> print(f"Found {result['tagged']} tagged reads out of {result['total']}")

    Output Files:
        - {prefix}.tag.fastq: FASTQ file containing tagged reads
        - {prefix}.nontag.fastq: FASTQ file containing non-tagged reads
        - Tag_statistics.txt: Summary statistics including:
            - Total reads analyzed
            - Number of tagged reads (with percentage)
            - Number of non-tagged reads
            - Number of corrupted reads skipped
    """
    log_msg(logfp, "=" * 60)
    log_msg(logfp, "TagSeek: Identifying tagged reads")
    log_msg(logfp, "=" * 60)

    # Define input and output paths
    input_fastq = Path(f"{input_fastq_prefix}.fastq")
    output_dir.mkdir(parents=True, exist_ok=True)

    output_tag = output_dir / f"{input_fastq_prefix}.tag.fastq"
    output_nontag = output_dir / f"{input_fastq_prefix}.nontag.fastq"
    output_stats = output_dir / "Tag_statistics.txt"

    # Validate input file exists
    ensure_exists(input_fastq, "Input FASTQ")

    # Generate search patterns from tag sequence
    patterns = generate_patterns(tag_sequence, similarity)
    log_msg(logfp, f"Generated {len(patterns)} search patterns")

    # Initialize counters
    total_reads = 0
    tagged_reads = 0
    skipped_reads = 0
    pattern_counts: Counter = Counter()

    log_msg(logfp, f"Processing file: {input_fastq}")

    # Process FASTQ file
    with open(output_tag, 'w') as tag_h, open(output_nontag, 'w') as nontag_h:
        # Manual FASTQ parsing to handle corrupted records
        with open(input_fastq, 'r') as fq_handle:
            while True:
                # Read 4 lines (FASTQ format)
                header = fq_handle.readline()
                if not header:
                    break  # End of file

                seq = fq_handle.readline().strip()
                plus = fq_handle.readline()
                qual = fq_handle.readline().strip()

                # Validate FASTQ format
                if not header.startswith('@') or not plus.startswith('+'):
                    skipped_reads += 1
                    continue

                # Check sequence and quality length match
                if len(seq) != len(qual):
                    skipped_reads += 1
                    continue

                total_reads += 1
                search_region = seq[:search_window]

                # Search for tag patterns
                found = False
                for pattern in patterns:
                    if regex.search(pattern, search_region):
                        found = True
                        pattern_counts[pattern] += 1
                        break

                # Write to appropriate output file
                if found:
                    tag_h.write(f"{header}{seq}\n+\n{qual}\n")
                    tagged_reads += 1
                else:
                    nontag_h.write(f"{header}{seq}\n+\n{qual}\n")

                # Progress reporting (every 100,000 reads)
                if total_reads % 100000 == 0:
                    log_msg(logfp, f"  Processed {total_reads:,} reads...")

    # Calculate statistics
    tag_ratio = tagged_reads / total_reads * 100 if total_reads > 0 else 0

    # Write statistics file
    with open(output_stats, 'w') as f:
        f.write("TagSeek Statistics Report\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Input file: {input_fastq}\n")
        f.write(f"Tag sequence: {tag_sequence}\n")
        f.write(f"Similarity threshold: {similarity}\n")
        f.write(f"Search window: {search_window} bp\n\n")
        f.write("-" * 50 + "\n")
        f.write(f"Total reads processed: {total_reads:,}\n")
        f.write(f"Corrupted reads skipped: {skipped_reads:,}\n")
        f.write(f"Tagged reads: {tagged_reads:,} ({tag_ratio:.2f}%)\n")
        f.write(f"Non-tagged reads: {total_reads - tagged_reads:,} ({100 - tag_ratio:.2f}%)\n")
        f.write("-" * 50 + "\n\n")
        
        if pattern_counts:
            f.write("Pattern detection counts:\n")
            for pattern, count in sorted(pattern_counts.items(), key=lambda x: -x[1]):
                f.write(f"  {pattern}: {count:,}\n")

    # Log summary
    log_msg(logfp, f"Complete! Tagged: {tagged_reads:,}/{total_reads:,} ({tag_ratio:.2f}%)")
    if skipped_reads > 0:
        log_msg(logfp, f"Skipped corrupted reads: {skipped_reads:,}")
    log_msg(logfp, f"Output files: {output_tag}, {output_nontag}")

    return {
        'total': total_reads,
        'tagged': tagged_reads,
        'skipped': skipped_reads,
        'tag_fastq': output_tag,
        'nontag_fastq': output_nontag,
        'stats_file': output_stats
    }
