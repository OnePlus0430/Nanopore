#!/usr/bin/env python3
"""
tagseqquant.py - Alignment and Quantification Module

This module handles the alignment of tagged and non-tagged reads to reference
genomes/transcriptomes using minimap2, and generates quantification results.
It integrates with samtools for BAM processing and FastQC for quality control.

Pipeline Steps:
    1. FastQC quality control (optional)
    2. minimap2 alignment to transcriptome (PAF format)
    3. minimap2 alignment to genome (PAF and SAM format)
    4. SAM to sorted BAM conversion with indexing
    5. Mapping statistics (samtools flagstat)
    6. Mapping QC plots (samtools stats + plot-bamstats)
    7. Python-based quantification (replaces R script)

Required External Tools:
    - minimap2: Sequence alignment
    - samtools: BAM file processing
    - fastqc: Quality control (optional)
    - plot-bamstats: Mapping QC visualization (optional)

Command-line usage:
    $ tagseqtools quant --name sample --trans ref.fa --genome genome.fa

"""

from __future__ import annotations
from tagseqtools.utils import *
from tagseqtools.quantification import run_quantification


def _run_minimap2_paf(
    ref: Path,
    fastq: Path,
    output_paf: Path,
    threads: int,
    logfp: Optional[TextIO]
) -> None:
    """
    Run minimap2 alignment and output in PAF format.

    PAF (Pairwise mApping Format) is a lightweight format for storing
    sequence alignments, suitable for quick overlap detection.

    Args:
        ref: Path to reference FASTA file (genome or transcriptome).
        fastq: Path to input FASTQ file.
        output_paf: Path for output PAF file.
        threads: Number of threads for minimap2.
        logfp: Optional file handle for logging.

    Raises:
        RuntimeError: If minimap2 execution fails.
    """
    cmd = [
        "minimap2",
        "-t", str(threads),
        "-cx", "splice",  # Spliced alignment preset
        "-uf",            # Transcript strand
        "-k14",           # k-mer size
        str(ref),
        str(fastq)
    ]
    
    log_msg(logfp, f"[PAF] {fastq.name} -> {output_paf.name}")

    paf_log = output_paf.with_suffix(output_paf.suffix + ".log")
    
    with open(output_paf, "w") as f, open(paf_log, "w") as lf:
        p = subprocess.run(cmd, stdout=f, stderr=lf, text=True)

    if p.returncode != 0:
        raise RuntimeError(f"minimap2 PAF alignment failed: {shlex.join(cmd)}")


def _run_minimap2_sam(
    ref: Path,
    fastq: Path,
    output_sam: Path,
    threads: int,
    logfp: Optional[TextIO]
) -> None:
    """
    Run minimap2 alignment and output in SAM format.

    SAM format is required for downstream BAM conversion and
    visualization in genome browsers like IGV.

    Args:
        ref: Path to reference FASTA file.
        fastq: Path to input FASTQ file.
        output_sam: Path for output SAM file.
        threads: Number of threads for minimap2.
        logfp: Optional file handle for logging.

    Raises:
        RuntimeError: If minimap2 execution fails.
    """
    cmd = [
        "minimap2",
        "-t", str(threads),
        "-ax", "splice",  # Spliced alignment, SAM output
        "-uf",            # Transcript strand
        "-k14",           # k-mer size
        str(ref),
        str(fastq)
    ]
    
    log_msg(logfp, f"[SAM] {fastq.name} -> {output_sam.name}")

    sam_log = output_sam.with_suffix(output_sam.suffix + ".log")
    
    with open(output_sam, "w") as f, open(sam_log, "w") as lf:
        p = subprocess.run(cmd, stdout=f, stderr=lf, text=True)

    if p.returncode != 0:
        raise RuntimeError(f"minimap2 SAM alignment failed: {shlex.join(cmd)}")


def _sam_to_sorted_bam(
    input_sam: Path,
    output_prefix: str,
    output_dir: Path,
    threads: int,
    logfp: Optional[TextIO]
) -> Path:
    """
    Convert SAM to sorted BAM with index.

    Performs three steps:
    1. SAM to BAM conversion (samtools view)
    2. BAM sorting (samtools sort)
    3. BAM indexing (samtools index)

    Args:
        input_sam: Path to input SAM file.
        output_prefix: Prefix for output files (e.g., "NAD" -> "NAD.sort.bam").
        output_dir: Directory for output files.
        threads: Number of threads for samtools sort.
        logfp: Optional file handle for logging.

    Returns:
        Path to the sorted and indexed BAM file.
    """
    bam = output_dir / f"{output_prefix}.bam"
    sorted_bam = output_dir / f"{output_prefix}.sort.bam"

    # SAM -> BAM
    log_msg(logfp, f"Converting SAM to BAM: {input_sam.name}")
    run_cmd(["samtools", "view", "-bS", str(input_sam), "-o", str(bam)], logfp)

    # Sort BAM
    log_msg(logfp, f"Sorting BAM: {bam.name}")
    run_cmd(["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), str(bam)], logfp)

    # Index BAM
    log_msg(logfp, f"Indexing BAM: {sorted_bam.name}")
    run_cmd(["samtools", "index", str(sorted_bam)], logfp)

    # Remove unsorted BAM
    bam.unlink(missing_ok=True)

    return sorted_bam


def _run_mapping_qc(
    bam: Path,
    label: str,
    output_dir: Path,
    logfp: Optional[TextIO]
) -> None:
    """
    Generate mapping QC statistics and plots.

    Uses samtools stats to generate statistics, then plot-bamstats
    to create visualization plots and HTML reports.

    Args:
        bam: Path to input BAM file.
        label: Label for output files (e.g., "NAD" or "nonNAD").
        output_dir: Directory for output files.
        logfp: Optional file handle for logging.
    """
    stats_file = output_dir / f"{label}.stats.txt"
    
    log_msg(logfp, f"Generating {label} mapping statistics...")
    
    # Generate stats
    run_cmd_to_file(
        ["samtools", "stats", str(bam.resolve())],
        stats_file,
        logfp,
        cwd=output_dir
    )
    
    # Generate plots
    run_cmd(
        ["plot-bamstats", "-p", f"{label}_map-", stats_file.name],
        logfp,
        cwd=output_dir
    )


def run_tagseqquant(
    name: str,
    trans_ref: Path,
    genome_ref: Path,
    outdir: Path,
    threads: int = 4,
    skip_fastqc: bool = False,
    skip_mapqc: bool = False,
    skip_flagstat: bool = False,
    keep_sam: bool = False,
    logfp: Optional[TextIO] = None
) -> None:
    """
    Run the complete TagSeqQuant alignment and quantification pipeline.

    This function orchestrates the entire quantification workflow:
    1. Quality control with FastQC
    2. Alignment to transcriptome and genome
    3. BAM file generation and processing
    4. Mapping statistics and QC plots
    5. Python-based quantification

    Args:
        name: Sample name/prefix. Expects files: {name}.fastq, {name}.tag.fastq,
              {name}.nontag.fastq in the current directory.
        trans_ref: Path to transcriptome reference FASTA file.
        genome_ref: Path to genome reference FASTA file.
        outdir: Output directory for all results.
        threads: Number of threads for parallel processing. Default is 4.
        skip_fastqc: If True, skip FastQC quality control. Default is False.
        skip_mapqc: If True, skip mapping QC plots. Default is False.
        skip_flagstat: If True, skip samtools flagstat. Default is False.
        keep_sam: If True, retain intermediate SAM files. Default is False.
        logfp: Optional file handle for logging.

    Raises:
        FileNotFoundError: If required input files or tools are missing.
        RuntimeError: If any pipeline step fails.

    Output Structure:
        outdir/
        ├── QC_results/           # FastQC reports
        ├── Mapping_statistics/   # PAF files and QC plots
        ├── Mapping_results/      # BAM files and mapping stats
        └── Quantification_results/  # Quantification output

    Example:
        >>> from pathlib import Path
        >>> run_tagseqquant(
        ...     name="sample",
        ...     trans_ref=Path("ref/transcriptome.fa"),
        ...     genome_ref=Path("ref/genome.fa"),
        ...     outdir=Path("results"),
        ...     threads=8
        ... )
    """
    log_msg(logfp, "")
    log_msg(logfp, "=" * 60)
    log_msg(logfp, "TagSeqQuant: Alignment and Quantification Pipeline")
    log_msg(logfp, "=" * 60)

    # Define input files
    all_fq = Path(f"{name}.fastq")
    tag_fq = Path(f"{name}.tag.fastq")
    nontag_fq = Path(f"{name}.nontag.fastq")

    # Validate input files
    ensure_exists(all_fq, "Complete FASTQ")
    ensure_exists(tag_fq, "Tagged FASTQ")
    ensure_exists(nontag_fq, "Non-tagged FASTQ")
    ensure_exists(trans_ref, "Transcriptome reference")
    ensure_exists(genome_ref, "Genome reference")

    # Check required tools
    log_msg(logfp, "Checking required tools...")
    ensure_tool_exists("minimap2")
    ensure_tool_exists("samtools")
    if not skip_fastqc:
        ensure_tool_exists("fastqc")
    if not skip_mapqc:
        ensure_tool_exists("plot-bamstats")

    # Create output directories
    qcdir = outdir / "QC_results"
    mapstatdir = outdir / "Mapping_statistics"
    mapresdir = outdir / "Mapping_results"
    quantdir = outdir / "Quantification_results"

    for d in [qcdir, mapstatdir, mapresdir, quantdir]:
        d.mkdir(parents=True, exist_ok=True)

    # =========================================================================
    # Step 1: FastQC
    # =========================================================================
    if not skip_fastqc:
        log_msg(logfp, "\n--- Step 1: FastQC Quality Control ---")
        run_cmd(
            ["fastqc", "-o", str(qcdir), str(all_fq), str(tag_fq), str(nontag_fq)],
            logfp
        )

    # =========================================================================
    # Step 2: minimap2 -> PAF (transcriptome and genome)
    # =========================================================================
    log_msg(logfp, "\n--- Step 2: Generating PAF Alignments ---")

    # Transcriptome PAF
    _run_minimap2_paf(trans_ref, tag_fq, mapstatdir / "trans.NAD.paf", threads, logfp)
    _run_minimap2_paf(trans_ref, all_fq, mapstatdir / "trans.all.paf", threads, logfp)
    _run_minimap2_paf(trans_ref, nontag_fq, mapstatdir / "trans.nonNAD.paf", threads, logfp)

    # Genome PAF
    _run_minimap2_paf(genome_ref, tag_fq, mapstatdir / "genome.NAD.paf", threads, logfp)
    _run_minimap2_paf(genome_ref, all_fq, mapstatdir / "genome.all.paf", threads, logfp)
    _run_minimap2_paf(genome_ref, nontag_fq, mapstatdir / "genome.nonNAD.paf", threads, logfp)

    # =========================================================================
    # Step 3: minimap2 -> SAM (genome and transcriptome)
    # =========================================================================
    log_msg(logfp, "\n--- Step 3: Generating SAM Alignments ---")

    genome_nad_sam = mapresdir / "genome.NAD.sam"
    genome_non_sam = mapresdir / "genome.nonNAD.sam"
    trans_nad_sam = mapresdir / "trans.NAD.sam"
    trans_non_sam = mapresdir / "trans.nonNAD.sam"

    _run_minimap2_sam(genome_ref, tag_fq, genome_nad_sam, threads, logfp)
    _run_minimap2_sam(genome_ref, nontag_fq, genome_non_sam, threads, logfp)
    _run_minimap2_sam(trans_ref, tag_fq, trans_nad_sam, threads, logfp)
    _run_minimap2_sam(trans_ref, nontag_fq, trans_non_sam, threads, logfp)

    # =========================================================================
    # Step 4: SAM -> Sorted BAM + Index
    # =========================================================================
    log_msg(logfp, "\n--- Step 4: Converting SAM to Sorted BAM ---")

    nad_bam = _sam_to_sorted_bam(genome_nad_sam, "NAD", mapresdir, threads, logfp)
    non_bam = _sam_to_sorted_bam(genome_non_sam, "nonNAD", mapresdir, threads, logfp)

    # Remove SAM files unless requested to keep
    if not keep_sam:
        for sam in [genome_nad_sam, genome_non_sam, trans_nad_sam, trans_non_sam]:
            sam.unlink(missing_ok=True)

    # =========================================================================
    # Step 5: Mapping Statistics (flagstat)
    # =========================================================================
    if not skip_flagstat:
        log_msg(logfp, "\n--- Step 5: Generating Mapping Statistics ---")
        run_cmd_to_file(
            ["samtools", "flagstat", str(nad_bam)],
            mapresdir / "NAD_Mapping_statistics.txt",
            logfp
        )
        run_cmd_to_file(
            ["samtools", "flagstat", str(non_bam)],
            mapresdir / "nonNAD_Mapping_statistics.txt",
            logfp
        )

    # =========================================================================
    # Step 6: Mapping QC Plots
    # =========================================================================
    if not skip_mapqc:
        log_msg(logfp, "\n--- Step 6: Generating Mapping QC Plots ---")
        _run_mapping_qc(nad_bam, "NAD", mapstatdir, logfp)
        _run_mapping_qc(non_bam, "nonNAD", mapstatdir, logfp)

    # =========================================================================
    # Step 7: Move Tag_statistics.txt if present
    # =========================================================================
    tagstat = Path("Tag_statistics.txt")
    if tagstat.exists():
        shutil.move(str(tagstat), str(mapresdir / "Tag_statistics.txt"))
        log_msg(logfp, f"Moved Tag_statistics.txt to {mapresdir}")

    # =========================================================================
    # Step 8: Copy PAF files to Quantification_results
    # =========================================================================
    for paf in ["trans.NAD.paf", "trans.nonNAD.paf"]:
        src = mapstatdir / paf
        if src.exists():
            shutil.copy2(src, quantdir / paf)

    # =========================================================================
    # Step 9: Run Python Quantification
    # =========================================================================
    log_msg(logfp, "\n--- Step 9: Running Quantification ---")
    run_quantification(
        nad_paf=quantdir / "trans.NAD.paf",
        nonnad_paf=quantdir / "trans.nonNAD.paf",
        output_dir=quantdir,
        logfp=logfp
    )

    # =========================================================================
    # Complete
    # =========================================================================
    log_msg(logfp, "\n" + "=" * 60)
    log_msg(logfp, "TagSeqQuant pipeline completed successfully!")
    log_msg(logfp, f"Output directory: {outdir}")
    log_msg(logfp, "=" * 60)
