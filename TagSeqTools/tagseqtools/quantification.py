#!/usr/bin/env python3
"""
quantification.py - Quantification Module for TagSeqTools

This module provides functionality to quantify NAD-capped RNAs from PAF alignment
files. It replaces the original R script (TagSeqQuant.r) with a pure Python
implementation.

Functions:
    - read_paf: Read and parse PAF alignment files
    - get_gene_name: Extract gene name from transcript ID
    - count_reads: Count reads per target
    - run_quantification: Main quantification pipeline

Input Files:
    - trans.NAD.paf: Transcriptome alignments of tagged reads
    - trans.nonNAD.paf: Transcriptome alignments of non-tagged reads

Output Files:
    - NAD_total_counts.txt: Gene-level counts
    - NAD_total_isoform_counts.txt: Isoform-level counts
    - Counting_statistics.txt: Summary statistics

Author: Huan ZHONG
License: Apache 2.0
"""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import Optional, TextIO, Dict, List, Tuple

from tagseqtools.utils import log_msg, ensure_exists


# PAF column indices
PAF_QUERY_NAME = 0
PAF_QUERY_LEN = 1
PAF_QUERY_START = 2
PAF_QUERY_END = 3
PAF_STRAND = 4
PAF_TARGET_NAME = 5
PAF_TARGET_LEN = 6
PAF_TARGET_START = 7
PAF_TARGET_END = 8
PAF_MATCHES = 9
PAF_BLOCK_LEN = 10
PAF_MAPQ = 11


def read_paf(paf_file: Path) -> List[Dict]:
    """
    Read PAF file and extract alignment information.

    PAF (Pairwise mApping Format) is a lightweight format for storing
    sequence alignments. This function parses the essential columns.

    Args:
        paf_file: Path to the PAF file.

    Returns:
        List of dictionaries, each containing alignment information with keys:
            - query_name: Read/query sequence name
            - query_len: Length of query sequence
            - query_start: Start position on query
            - query_end: End position on query
            - strand: Mapping strand (+ or -)
            - target_name: Reference sequence name (transcript/gene ID)

    Raises:
        FileNotFoundError: If the PAF file does not exist.

    Example:
        >>> alignments = read_paf(Path("trans.NAD.paf"))
        >>> print(f"Found {len(alignments)} alignments")
    """
    ensure_exists(paf_file, "PAF file")

    alignments = []
    with open(paf_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) < 6:
                continue

            alignment = {
                'query_name': fields[PAF_QUERY_NAME],
                'query_len': int(fields[PAF_QUERY_LEN]) if fields[PAF_QUERY_LEN].isdigit() else 0,
                'query_start': int(fields[PAF_QUERY_START]) if fields[PAF_QUERY_START].isdigit() else 0,
                'query_end': int(fields[PAF_QUERY_END]) if fields[PAF_QUERY_END].isdigit() else 0,
                'strand': fields[PAF_STRAND],
                'target_name': fields[PAF_TARGET_NAME],
            }
            alignments.append(alignment)

    return alignments


def get_gene_name(transcript_id: str) -> str:
    """
    Extract gene name from transcript ID by removing isoform suffix.

    Removes trailing version/isoform numbers (e.g., ".1", ".2") from
    transcript identifiers to aggregate counts at the gene level.

    Args:
        transcript_id: Transcript identifier (e.g., "AT1G01100.1").

    Returns:
        Gene name without isoform suffix (e.g., "AT1G01100").

    Example:
        >>> get_gene_name("AT1G01100.1")
        'AT1G01100'
        >>> get_gene_name("ENST00000456328.2")
        'ENST00000456328'
        >>> get_gene_name("gene_without_suffix")
        'gene_without_suffix'
    """
    return re.sub(r'\.[0-9]+$', '', transcript_id)


def count_reads(alignments: List[Dict], level: str = 'isoform') -> Counter:
    """
    Count reads per target (gene or isoform level).

    Args:
        alignments: List of alignment dictionaries from read_paf().
        level: Counting level - 'isoform' for transcript-level or
               'gene' for gene-level (aggregates isoforms).

    Returns:
        Counter object with target names as keys and read counts as values.

    Example:
        >>> alignments = read_paf(Path("trans.NAD.paf"))
        >>> isoform_counts = count_reads(alignments, level='isoform')
        >>> gene_counts = count_reads(alignments, level='gene')
    """
    counts = Counter()

    for aln in alignments:
        target = aln['target_name']
        if level == 'gene':
            target = get_gene_name(target)
        counts[target] += 1

    return counts


def merge_counts(
    nad_counts: Counter,
    total_counts: Counter
) -> List[Tuple[str, int, int]]:
    """
    Merge NAD and total counts into a combined table.

    Args:
        nad_counts: Counter with NAD-tagged read counts.
        total_counts: Counter with total read counts.

    Returns:
        List of tuples (target_name, nad_count, total_count) sorted by
        NAD count in descending order.
    """
    all_targets = set(nad_counts.keys()) | set(total_counts.keys())

    merged = []
    for target in all_targets:
        nad = nad_counts.get(target, 0)
        total = total_counts.get(target, 0)
        merged.append((target, nad, total))

    # Sort by NAD count (descending), then by target name
    merged.sort(key=lambda x: (-x[1], x[0]))

    return merged


def write_counts_table(
    counts: List[Tuple[str, int, int]],
    output_file: Path,
    header: Tuple[str, str, str] = ("Gene", "NAD.count", "total.count")
) -> None:
    """
    Write counts table to a tab-separated file.

    Args:
        counts: List of tuples (target, nad_count, total_count).
        output_file: Path to output file.
        header: Column headers for the output file.
    """
    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for target, nad, total in counts:
            f.write(f"{target}\t{nad}\t{total}\n")


def run_quantification(
    nad_paf: Path,
    nonnad_paf: Path,
    output_dir: Path,
    logfp: Optional[TextIO] = None
) -> Dict:
    """
    Run the complete quantification pipeline.

    Reads PAF alignment files, counts reads at both isoform and gene levels,
    and generates output files with quantification results.

    Args:
        nad_paf: Path to NAD-tagged reads PAF file (trans.NAD.paf).
        nonnad_paf: Path to non-tagged reads PAF file (trans.nonNAD.paf).
        output_dir: Directory for output files.
        logfp: Optional file handle for logging.

    Returns:
        Dictionary containing summary statistics:
            - total_alignments: Total number of alignments
            - nad_alignments: Number of NAD alignments
            - nonnad_alignments: Number of non-NAD alignments
            - total_genes: Number of genes detected
            - nad_genes: Number of NAD-capped genes
            - total_isoforms: Number of isoforms detected
            - nad_isoforms: Number of NAD-capped isoforms

    Raises:
        FileNotFoundError: If input PAF files do not exist.

    Output Files:
        - NAD_total_counts.txt: Gene-level counts with columns:
            Gene, NAD.count, total.count
        - NAD_total_isoform_counts.txt: Isoform-level counts with columns:
            Gene, NAD.count, total.count
        - Counting_statistics.txt: Summary statistics

    Example:
        >>> stats = run_quantification(
        ...     nad_paf=Path("trans.NAD.paf"),
        ...     nonnad_paf=Path("trans.nonNAD.paf"),
        ...     output_dir=Path("results")
        ... )
        >>> print(f"Found {stats['nad_genes']} NAD-capped genes")
    """
    log_msg(logfp, "")
    log_msg(logfp, "TagSeqQuant Python Quantification")
    log_msg(logfp, "=" * 40)

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # Read PAF files
    log_msg(logfp, "Reading PAF files...")
    nad_alignments = read_paf(nad_paf)
    nonnad_alignments = read_paf(nonnad_paf)
    all_alignments = nad_alignments + nonnad_alignments

    log_msg(logfp, f"  NAD alignments: {len(nad_alignments):,}")
    log_msg(logfp, f"  nonNAD alignments: {len(nonnad_alignments):,}")

    # Count reads per isoform
    log_msg(logfp, "Counting reads per isoform...")
    nad_isoform_counts = count_reads(nad_alignments, level='isoform')
    total_isoform_counts = count_reads(all_alignments, level='isoform')

    # Merge and write isoform counts
    isoform_merged = merge_counts(nad_isoform_counts, total_isoform_counts)
    isoform_output = output_dir / "NAD_total_isoform_counts.txt"
    write_counts_table(isoform_merged, isoform_output)
    log_msg(logfp, f"  Wrote {isoform_output.name}: {len(isoform_merged):,} isoforms")

    # Count reads per gene (aggregate isoforms)
    log_msg(logfp, "Aggregating to gene level...")
    nad_gene_counts = count_reads(nad_alignments, level='gene')
    total_gene_counts = count_reads(all_alignments, level='gene')

    # Merge and write gene counts
    gene_merged = merge_counts(nad_gene_counts, total_gene_counts)
    gene_output = output_dir / "NAD_total_counts.txt"
    write_counts_table(gene_merged, gene_output)
    log_msg(logfp, f"  Wrote {gene_output.name}: {len(gene_merged):,} genes")

    # Calculate statistics
    stats = {
        'total_alignments': len(all_alignments),
        'nad_alignments': len(nad_alignments),
        'nonnad_alignments': len(nonnad_alignments),
        'total_genes': len(gene_merged),
        'nad_genes': sum(1 for _, nad, _ in gene_merged if nad > 0),
        'total_isoforms': len(isoform_merged),
        'nad_isoforms': sum(1 for _, nad, _ in isoform_merged if nad > 0),
    }

    # Write statistics file
    stats_output = output_dir / "Counting_statistics.txt"
    with open(stats_output, 'w') as f:
        f.write("Metric\tValue\n")
        f.write(f"Total_alignments\t{stats['total_alignments']}\n")
        f.write(f"NAD_alignments\t{stats['nad_alignments']}\n")
        f.write(f"nonNAD_alignments\t{stats['nonnad_alignments']}\n")
        f.write(f"Total_genes\t{stats['total_genes']}\n")
        f.write(f"NAD_genes\t{stats['nad_genes']}\n")
        f.write(f"Total_isoforms\t{stats['total_isoforms']}\n")
        f.write(f"NAD_isoforms\t{stats['nad_isoforms']}\n")

    log_msg(logfp, f"  Wrote {stats_output.name}")

    # Print summary
    log_msg(logfp, "")
    log_msg(logfp, "Summary Statistics:")
    log_msg(logfp, "-" * 30)
    log_msg(logfp, f"  Total alignments: {stats['total_alignments']:,}")
    log_msg(logfp, f"  NAD alignments: {stats['nad_alignments']:,}")
    log_msg(logfp, f"  nonNAD alignments: {stats['nonnad_alignments']:,}")
    log_msg(logfp, f"  Total genes: {stats['total_genes']:,}")
    log_msg(logfp, f"  NAD genes: {stats['nad_genes']:,}")
    log_msg(logfp, f"  Total isoforms: {stats['total_isoforms']:,}")
    log_msg(logfp, f"  NAD isoforms: {stats['nad_isoforms']:,}")
    log_msg(logfp, "")
    log_msg(logfp, "Quantification completed successfully!")

    return stats
