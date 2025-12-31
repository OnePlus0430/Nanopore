"""
TagSeqTools - A Comprehensive Pipeline for NAD-capped RNA Analysis

TagSeqTools is a bioinformatics pipeline designed for analyzing NAD tagSeq data
from Nanopore sequencing. It provides tools for identifying tagged (NAD-capped)
RNA reads and performing quantitative analysis.

Modules:
    - tagseek: Identifies and separates tagged/non-tagged reads
    - tagseqquant: Performs alignment and quantification
    - quantification: PAF-based read counting and statistics
    - utils: Shared utility functions
    - cli: Command-line interface

For command-line usage:
    $ tagseqtools full --fastq sample --tag 'CCTGAA...' -s 12 --trans ref.fa --genome genome.fa

"""

from tagseqtools.tagseek import *
from tagseqtools.tagseqquant import *
from tagseqtools.quantification import *
from tagseqtools.utils import *

__all__ = [
    # tagseek
    "run_tagseek",
    "generate_patterns",
    # tagseqquant
    "run_tagseqquant",
    # quantification
    "run_quantification",
    "read_paf",
    "get_gene_name",
    "count_reads",
    # utils
    "log_msg",
    "ensure_exists",
    "ensure_tool_exists",
]
