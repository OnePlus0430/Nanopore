#!/usr/bin/env python3
"""
cli.py - Command-Line Interface for TagSeqTools

This module provides the command-line interface for the TagSeqTools pipeline.
It uses argparse to handle subcommands and arguments, orchestrating the
workflow between TagSeek and TagSeqQuant modules.

Subcommands:
    - full: Run complete pipeline (TagSeek + TagSeqQuant)
    - tagseek: Run only tag identification
    - quant: Run only alignment and quantification

Usage:
    $ tagseqtools full --fastq sample --tag 'CCTGAA...' -s 12 --trans ref.fa --genome genome.fa
    $ tagseqtools tagseek --fastq sample --tag 'CCTGAA...' -s 12
    $ tagseqtools quant --name sample --trans ref.fa --genome genome.fa

The CLI can be invoked either through the installed entry point 'tagseqtools'
or directly by running this module:
    $ python -m tagseqtools.cli full --fastq sample ...

"""

from __future__ import annotations

import argparse
import sys
from typing import Optional

from tagseqtools.tagseek import *
from tagseqtools.tagseqquant import *
from tagseqtools.utils import *


def create_parser() -> argparse.ArgumentParser:
    """
    Create and configure the argument parser with all subcommands.

    Returns:
        Configured ArgumentParser instance with all subcommands and arguments.
    """
    parser = argparse.ArgumentParser(
        prog="tagseqtools",
        description="TagSeqTools - A comprehensive pipeline for NAD-capped RNA analysis",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run complete pipeline
  tagseqtools full --fastq sample --tag 'CCTGAACCTGAA...' -s 12 \\
      --trans transcriptome.fa --genome genome.fa --threads 8

  # Run TagSeek only (identify tagged reads)
  tagseqtools tagseek --fastq sample --tag 'CCTGAACCTGAA...' -s 12

  # Run TagSeqQuant only (alignment and quantification)
  tagseqtools quant --name sample --trans transcriptome.fa --genome genome.fa

For more information, visit: https://github.com/dorothyzh/TagSeqTools
        """
    )

    # parser.add_argument(
    #     "--version", "-v",
    #     action="version",
    #     version=f"%(prog)s {__version__}"
    # )

    subparsers = parser.add_subparsers(
        dest="command",
        title="Commands",
        description="Available pipeline commands"
    )

    # =========================================================================
    # 'full' subcommand - Complete pipeline
    # =========================================================================
    full_parser = subparsers.add_parser(
        "full",
        help="Run complete pipeline: TagSeek + TagSeqQuant",
        description="Run the complete TagSeqTools pipeline including tag identification, "
                   "alignment, and quantification."
    )
    
    full_parser.add_argument(
        "--fastq",
        required=True,
        help="FASTQ file prefix (without .fastq extension)"
    )
    full_parser.add_argument(
        "--tag",
        default="CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT",
        help="Synthetic tag RNA/DNA sequence to search for (default: CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT)"
    )
    full_parser.add_argument(
        "--similarity", "-s",
        type=int,
        required=True,
        help="Number of consecutive matching bases required for tag detection"
    )
    full_parser.add_argument(
        "--trans",
        required=True,
        help="Path to transcriptome reference FASTA file"
    )
    full_parser.add_argument(
        "--genome",
        required=True,
        help="Path to genome reference FASTA file"
    )
    full_parser.add_argument(
        "--outdir",
        default="TagSeqTools_out",
        help="Output directory (default: TagSeqTools_out)"
    )
    full_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for alignment (default: 4)"
    )
    full_parser.add_argument(
        "--skip-fastqc",
        action="store_true",
        help="Skip FastQC quality control step"
    )
    full_parser.add_argument(
        "--skip-mapqc",
        action="store_true",
        help="Skip mapping QC plot generation"
    )
    full_parser.add_argument(
        "--keep-sam",
        action="store_true",
        help="Keep intermediate SAM files (default: delete)"
    )

    # =========================================================================
    # 'tagseek' subcommand - Tag identification only
    # =========================================================================
    tagseek_parser = subparsers.add_parser(
        "tagseek",
        help="Run TagSeek only: identify and separate tagged/non-tagged reads",
        description="Identify reads containing the synthetic tag sequence and "
                   "separate them into tagged and non-tagged FASTQ files."
    )
    
    tagseek_parser.add_argument(
        "--fastq", "-f",
        required=True,
        help="FASTQ file prefix (without .fastq extension)"
    )
    tagseek_parser.add_argument(
        "--tag", "-t",
        default="CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT",
        help="Synthetic tag RNA/DNA sequence to search for (default: CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT)"
    )
    tagseek_parser.add_argument(
        "--similarity", "-s",
        type=int,
        required=True,
        help="Number of consecutive matching bases required"
    )
    tagseek_parser.add_argument(
        "--window", "-w",
        type=int,
        default=40,
        help="Search window size in base pairs (default: 40)"
    )
    tagseek_parser.add_argument(
        "--outdir", "-o",
        default=".",
        help="Output directory (default: current directory)"
    )

    # =========================================================================
    # 'quant' subcommand - Alignment and quantification only
    # =========================================================================
    quant_parser = subparsers.add_parser(
        "quant",
        help="Run TagSeqQuant only: alignment and quantification",
        description="Perform alignment to reference sequences and generate "
                   "quantification results. Requires pre-existing .tag.fastq "
                   "and .nontag.fastq files."
    )
    
    quant_parser.add_argument(
        "--name", "-n",
        required=True,
        help="Sample name/prefix (expects {name}.fastq, {name}.tag.fastq, {name}.nontag.fastq)"
    )
    quant_parser.add_argument(
        "--trans",
        required=True,
        help="Path to transcriptome reference FASTA file"
    )
    quant_parser.add_argument(
        "--genome",
        required=True,
        help="Path to genome reference FASTA file"
    )
    quant_parser.add_argument(
        "--outdir",
        default="TagSeqQuant_out",
        help="Output directory (default: TagSeqQuant_out)"
    )
    quant_parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="Number of threads for alignment (default: 4)"
    )
    quant_parser.add_argument(
        "--skip-fastqc",
        action="store_true",
        help="Skip FastQC quality control step"
    )
    quant_parser.add_argument(
        "--skip-mapqc",
        action="store_true",
        help="Skip mapping QC plot generation"
    )
    quant_parser.add_argument(
        "--keep-sam",
        action="store_true",
        help="Keep intermediate SAM files (default: delete)"
    )

    return parser


def main(argv: Optional[list] = None) -> int:
    """
    Main entry point for the TagSeqTools command-line interface.

    Parses command-line arguments and dispatches to the appropriate
    pipeline function based on the subcommand.

    Args:
        argv: Command-line arguments. If None, uses sys.argv.

    Returns:
        Exit code: 0 for success, 1 for error.
    """
    parser = create_parser()
    args = parser.parse_args(argv)

    # Show help if no command provided
    if not args.command:
        parser.print_help()
        return 1

    # Print banner
    print("=" * 60)
    print("NAD-capped RNA Analysis Pipeline")
    print("=" * 60)
    print()

    # Setup output directory and log file
    outdir = Path(args.outdir if hasattr(args, 'outdir') else '.')
    outdir.mkdir(parents=True, exist_ok=True)
    logpath = outdir / "tagseqtools.log"

    # Execute pipeline
    with open(logpath, 'w') as logfp:
        try:
            if args.command == 'tagseek':
                # Run TagSeek only
                run_tagseek(
                    input_fastq_prefix=args.fastq,
                    tag_sequence=args.tag,
                    similarity=args.similarity,
                    search_window=args.window,
                    output_dir=Path(args.outdir),
                    logfp=logfp
                )

            elif args.command == 'quant':
                # Run TagSeqQuant only
                run_tagseqquant(
                    name=args.name,
                    trans_ref=Path(args.trans),
                    genome_ref=Path(args.genome),
                    outdir=outdir,
                    threads=args.threads,
                    skip_fastqc=args.skip_fastqc,
                    skip_mapqc=args.skip_mapqc,
                    skip_flagstat=False,
                    keep_sam=args.keep_sam,
                    logfp=logfp
                )

            elif args.command == 'full':
                # Run complete pipeline
                
                # Step 1: TagSeek
                run_tagseek(
                    input_fastq_prefix=args.fastq,
                    tag_sequence=args.tag,
                    similarity=args.similarity,
                    search_window=40,
                    output_dir=Path('.'),
                    logfp=logfp
                )

                # Step 2: TagSeqQuant
                run_tagseqquant(
                    name=args.fastq,
                    trans_ref=Path(args.trans),
                    genome_ref=Path(args.genome),
                    outdir=outdir,
                    threads=args.threads,
                    skip_fastqc=args.skip_fastqc,
                    skip_mapqc=args.skip_mapqc,
                    skip_flagstat=False,
                    keep_sam=args.keep_sam,
                    logfp=logfp
                )

            log_msg(logfp, f"\nLog file: {logpath}")
            print(f"\n[SUCCESS] Pipeline completed. Log: {logpath}")
            return 0

        except FileNotFoundError as e:
            log_msg(logfp, f"\n[ERROR] File not found: {e}")
            print(f"\n[ERROR] {e}", file=sys.stderr)
            return 1

        except RuntimeError as e:
            log_msg(logfp, f"\n[ERROR] Pipeline failed: {e}")
            print(f"\n[ERROR] {e}", file=sys.stderr)
            return 1

        except Exception as e:
            log_msg(logfp, f"\n[ERROR] Unexpected error: {e}")
            import traceback
            traceback.print_exc()
            return 1


if __name__ == "__main__":
    sys.exit(main())
