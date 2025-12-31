#!/usr/bin/env python3
"""
setup.py - Installation script for TagSeqTools

This script handles the installation of the TagSeqTools package using setuptools.
It can be used for both development and production installations.

Installation:
    # Standard installation
    pip install .

    # Development installation (editable)
    pip install -e .

    # Install with all dependencies
    pip install -e ".[dev]"

Author: Huan ZHONG
License: Apache 2.0
"""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README for long description
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text(encoding="utf-8")

# Read requirements
requirements = [
    "biopython>=1.78",
    "regex>=2021.8.3",
]

dev_requirements = [
    "pytest>=6.0",
    "pytest-cov>=2.0",
    "black>=21.0",
    "flake8>=3.9",
]

setup(
    name="tagseqtools",
    version="1.0.0",
    author="Huan ZHONG",
    author_email="zhdorothy5@gmail.com",
    description="A comprehensive pipeline for NAD-capped RNA analysis using Nanopore sequencing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/dorothyzh/TagSeqTools",
    project_urls={
        "Bug Tracker": "https://github.com/dorothyzh/TagSeqTools/issues",
        "Documentation": "https://github.com/dorothyzh/TagSeqTools#readme",
        "Source Code": "https://github.com/dorothyzh/TagSeqTools",
    },
    license="Apache-2.0",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    keywords="bioinformatics, nanopore, rna-seq, nad-rna, transcriptomics",
    packages=find_packages(exclude=["tests", "tests.*", "examples", "examples.*"]),
    package_data={
        "tagseqtools": ["r_scripts/*.r", "r_scripts/*.R"],
    },
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=requirements,
    extras_require={
        "dev": dev_requirements,
    },
    entry_points={
        "console_scripts": [
            "tagseqtools=tagseqtools.cli:main",
        ],
    },
    zip_safe=False,
)
