#!/usr/bin/env python3
"""
utils.py - Utility Functions for TagSeqTools

This module provides shared utility functions used across the TagSeqTools pipeline,
including logging, file validation, and command execution helpers.

Functions:
    - log_msg: Write timestamped messages to log file and console
    - ensure_exists: Validate that a file or directory exists
    - ensure_tool_exists: Check if a command-line tool is available in PATH
    - run_cmd: Execute a shell command with logging
    - run_cmd_to_file: Execute a command and redirect stdout to a file

"""

from __future__ import annotations

import shlex
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional, TextIO, List, Union


def log_msg(
    logfp: Optional[TextIO],
    msg: str,
    also_print: bool = True
) -> None:
    """
    Write a timestamped message to log file and optionally print to console.

    Args:
        logfp: File handle for the log file. If None, only prints to console.
        msg: The message to log.
        also_print: If True, also print the message to stdout. Default is True.

    Returns:
        None

    Example:
        >>> with open("run.log", "w") as log:
        ...     log_msg(log, "Processing started")
        ...     log_msg(log, "Step 1 complete", also_print=False)
    """
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    line = f"[{timestamp}] {msg}"
    
    if logfp:
        logfp.write(line + "\n")
        logfp.flush()
    
    if also_print:
        print(msg)


def ensure_exists(p: Path, kind: str) -> None:
    """
    Check if a file or directory exists, raise FileNotFoundError if not.

    Args:
        p: Path object pointing to the file or directory to check.
        kind: A descriptive name for the file type (used in error message).

    Raises:
        FileNotFoundError: If the path does not exist.

    Example:
        >>> from pathlib import Path
        >>> ensure_exists(Path("genome.fa"), "Genome reference")  # Raises if missing
    """
    if not p.exists():
        raise FileNotFoundError(f"{kind} not found: {p}")


def ensure_tool_exists(tool: str) -> None:
    """
    Check if a command-line tool is available in the system PATH.

    Args:
        tool: Name of the command-line tool to check (e.g., "minimap2", "samtools").

    Raises:
        FileNotFoundError: If the tool is not found in PATH.

    Example:
        >>> ensure_tool_exists("minimap2")  # Raises if minimap2 not installed
        >>> ensure_tool_exists("samtools")  # Raises if samtools not installed
    """
    if shutil.which(tool) is None:
        raise FileNotFoundError(
            f"Required tool not found in PATH: {tool}\n"
            f"Please install {tool} and ensure it is available in your PATH."
        )


def run_cmd(
    cmd: List[str],
    logfp: Optional[TextIO],
    cwd: Optional[Path] = None
) -> None:
    """
    Execute a shell command with real-time output logging.

    Runs a command using subprocess, streaming stdout and stderr to the log file.
    Raises an exception if the command fails (non-zero exit code).

    Args:
        cmd: List of command arguments (e.g., ["minimap2", "-t", "4", "ref.fa", "reads.fq"]).
        logfp: File handle for logging. If None, output is not logged.
        cwd: Working directory for command execution. If None, uses current directory.

    Raises:
        RuntimeError: If the command exits with a non-zero status.

    Example:
        >>> with open("run.log", "w") as log:
        ...     run_cmd(["samtools", "view", "-bS", "input.sam"], log)
    """
    if logfp:
        logfp.write("\n" + "=" * 80 + "\n")
        logfp.write(f"[{datetime.now().isoformat(timespec='seconds')}] CMD: {shlex.join(cmd)}\n")
        logfp.flush()

    p = subprocess.Popen(
        cmd,
        cwd=str(cwd) if cwd else None,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    assert p.stdout is not None
    for line in p.stdout:
        if logfp:
            logfp.write(line)
    p.wait()

    if p.returncode != 0:
        raise RuntimeError(f"Command failed (exit {p.returncode}): {shlex.join(cmd)}")


def run_cmd_to_file(
    cmd: List[str],
    outfile: Path,
    logfp: Optional[TextIO],
    cwd: Optional[Path] = None
) -> None:
    """
    Execute a command and redirect stdout to a file.

    Runs a command, capturing stdout to the specified output file while
    logging stderr to the log file. Useful for commands like `samtools stats`
    that output results to stdout.

    Args:
        cmd: List of command arguments.
        outfile: Path to the output file where stdout will be written.
        logfp: File handle for logging stderr. If None, stderr is not logged.
        cwd: Working directory for command execution. If None, uses current directory.

    Raises:
        RuntimeError: If the command exits with a non-zero status.

    Example:
        >>> with open("run.log", "w") as log:
        ...     run_cmd_to_file(
        ...         ["samtools", "flagstat", "aligned.bam"],
        ...         Path("stats.txt"),
        ...         log
        ...     )
    """
    if logfp:
        logfp.write("\n" + "=" * 80 + "\n")
        logfp.write(f"[{datetime.now().isoformat(timespec='seconds')}] CMD: {shlex.join(cmd)}\n")
        logfp.write(f"[{datetime.now().isoformat(timespec='seconds')}] STDOUT -> {outfile}\n")
        logfp.flush()

    outfile.parent.mkdir(parents=True, exist_ok=True)
    
    with open(outfile, "w", encoding="utf-8") as out:
        p = subprocess.run(
            cmd,
            cwd=str(cwd) if cwd else None,
            stdout=out,
            stderr=subprocess.PIPE,
            text=True,
        )

    if p.stderr and logfp:
        logfp.write(p.stderr)
        logfp.flush()

    if p.returncode != 0:
        raise RuntimeError(f"Command failed (exit {p.returncode}): {shlex.join(cmd)}")


def format_number(n: int) -> str:
    """
    Format a number with thousand separators for readable output.

    Args:
        n: Integer to format.

    Returns:
        String representation with comma separators.

    Example:
        >>> format_number(1234567)
        '1,234,567'
    """
    return f"{n:,}"
