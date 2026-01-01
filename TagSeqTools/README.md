# TagSeqTools

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-green.svg)](https://opensource.org/licenses/Apache-2.0)

A comprehensive bioinformatics pipeline for analyzing NAD-capped RNA using Nanopore sequencing data. TagSeqTools identifies tagged (NAD-capped) and non-tagged RNA reads, performs alignment, and generates quantification results.

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Full Pipeline](#full-pipeline)
  - [TagSeek Only](#tagseek-only)
  - [TagSeqQuant Only](#tagseqquant-only)
- [Output Structure](#output-structure)
- [Project Structure](#project-structure)
- [Requirements](#requirements)
- [Citation](#citation)
- [License](#license)

## Introduction

TagSeqTools is designed for analyzing NAD tagSeq data from Nanopore sequencing. The pipeline consists of two main modules:

1. **TagSeek**: Differentiates tagged (NAD-capped) and non-tagged RNA reads based on synthetic tag sequence matching
2. **TagSeqQuant**: Performs alignment to reference genomes/transcriptomes and generates quantification results

## Features

- **Tag Identification**: Fast and accurate identification of tagged reads using pattern matching
- **Robust FASTQ Parsing**: Handles corrupted or malformed FASTQ records gracefully
- **Comprehensive Alignment**: Generates both PAF and BAM alignments using minimap2
- **Quality Control**: Integrates FastQC for read quality assessment
- **Mapping Statistics**: Produces detailed mapping statistics and visualization
- **Quantification**: Gene and isoform-level quantification of NAD-capped RNAs
- **Flexible Workflow**: Run the complete pipeline or individual modules

## Installation

### Prerequisites

**System Requirements:**
- Linux-based operating system (Ubuntu 18.04+ recommended)
- Python 3.8 or higher
- R 3.2.1 or higher

**External Tools (must be in PATH):**
- [minimap2](https://github.com/lh3/minimap2) >= 2.12
- [samtools](http://www.htslib.org/) >= 1.7
- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) >= 0.11.4

### Install from Source

```bash
# Clone the repository
git clone https://github.com/yourusername/TagSeqTools.git
cd TagSeqTools

# Install Python dependencies
pip install -r requirements.txt

# Install the package
pip install -e .
```

### Install External Tools

```bash
# Ubuntu/Debian
sudo apt-get install samtools fastqc

# Install minimap2
curl -L https://github.com/lh3/minimap2/releases/download/v2.24/minimap2-2.24_x64-linux.tar.bz2 | tar -jxvf -
export PATH=$PATH:$(pwd)/minimap2-2.24_x64-linux/
```

## Quick Start

```bash
# Run the complete pipeline
tagseqtools full \
    --fastq sample \
    --tag 'CCTGAACCTGAACCTGAACCTGAACCTGAACCTGAACCT' \
    --similarity 12 \
    --trans reference/transcriptome.fa \
    --genome reference/genome.fa \
    --threads 8 \
    --outdir results
```

## Usage

### Full Pipeline

Run both TagSeek and TagSeqQuant in sequence:

```bash
tagseqtools full \
    --fastq <prefix> \
    --tag <tag_sequence> \
    --similarity <threshold> \
    --trans <transcriptome.fa> \
    --genome <genome.fa> \
    [--rscript <TagSeqQuant.r>] \
    [--threads <num>] \
    [--outdir <output_dir>] \
    [--skip-fastqc] \
    [--skip-mapqc] \
    [--keep-sam]
```

**Parameters:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--fastq` | Yes | - | FASTQ file prefix (without .fastq extension) |
| `--tag` | Yes | - | Synthetic tag RNA sequence |
| `--similarity` | Yes | - | Number of consecutive matched bases for tag detection |
| `--trans` | Yes | - | Transcriptome reference FASTA file |
| `--genome` | Yes | - | Genome reference FASTA file |
| `--rscript` | No | TagSeqQuant.r | Path to R quantification script |
| `--threads` | No | 4 | Number of threads for alignment |
| `--outdir` | No | TagSeqTools_out | Output directory |
| `--skip-fastqc` | No | False | Skip FastQC analysis |
| `--skip-mapqc` | No | False | Skip mapping QC plots |
| `--keep-sam` | No | False | Keep intermediate SAM files |

### TagSeek Only

Identify and separate tagged/non-tagged reads:

```bash
tagseqtools tagseek \
    --fastq <prefix> \
    --tag <tag_sequence> \
    --similarity <threshold> \
    [--window <search_window>] \
    [--outdir <output_dir>]
```

**Parameters:**

| Parameter | Required | Default | Description |
|-----------|----------|---------|-------------|
| `--fastq`, `-f` | Yes | - | FASTQ file prefix |
| `--tag`, `-t` | Yes | - | Tag sequence to search for |
| `--similarity`, `-s` | Yes | - | Similarity threshold (consecutive matches) |
| `--window`, `-w` | No | 40 | Search window size (first N bases of read) |
| `--outdir`, `-o` | No | . | Output directory |

### TagSeqQuant Only

Run alignment and quantification (requires pre-existing .tag.fastq and .nontag.fastq):

```bash
tagseqtools quant \
    --name <prefix> \
    --trans <transcriptome.fa> \
    --genome <genome.fa> \
    [--rscript <TagSeqQuant.r>] \
    [--threads <num>] \
    [--outdir <output_dir>]
```

### Module Descriptions

| Module | Description |
|--------|-------------|
| `cli.py` | Command-line interface handling argument parsing and workflow orchestration |
| `tagseek.py` | Tag identification module - separates tagged and non-tagged reads |
| `tagseqquant.py` | Alignment and quantification module using minimap2 and samtools |
| `utils.py` | Shared utility functions for logging, file validation, and command execution |

## Requirements

### Python Packages

- biopython >= 1.78
- regex >= 2021.8.3

### External Tools

| Tool | Version | Purpose |
|------|---------|---------|
| minimap2 | >= 2.12 | Sequence alignment |
| samtools | >= 1.7 | BAM file processing |
| FastQC | >= 0.11.4 | Read quality control |

## Citation

If you use TagSeqTools in your research, please cite:

> Zhang, H.*, Zhong, H.*, Zhang, S., Shao, X., Ni, M., Cai, Z., Chen, X., and Xia, Y. (2019). 
> NAD TagSeq Reveals That NAD+-Capped RNAs Are Mostly Produced from a Large Number of 
> Protein-Coding Genes in Arabidopsis. *Proceedings of the National Academy of Sciences*.
> https://doi.org/10.1073/pnas.1903683116

> Zhong, H., Cai, Z., Yang, Z., and Xia, Y. (2020). 
> TagSeqTools: A flexible and comprehensive analysis pipeline for NAD tagSeq data. 
> *bioRxiv*. https://doi.org/10.1101/2020.03.09.982934

## License

This project is licensed under the Apache License 2.0 - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Author**: Huan ZHONG
- **Email**: zhdorothy5@gmail.com
- **GitHub**: [https://github.com/dorothyzh/TagSeqTools](https://github.com/dorothyzh/TagSeqTools)
