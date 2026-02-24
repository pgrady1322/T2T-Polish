#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
T2T-Polish v4.0

Automated Polishing Pipeline v4 — iterative assembly polishing with Merfin and Winnowmap.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

import sys
sys.stdout.flush()

import pysam  # type: ignore
import os
import argparse
import subprocess
from subprocess import PIPE
import time
import json
import concurrent.futures
import re
import shutil
import logging
from tqdm import tqdm
from typing import Tuple
from contextlib import redirect_stdout, redirect_stderr
from collections import defaultdict, Counter
from pathlib import Path
import subprocess

logger = logging.getLogger(__name__)

# Adjust the calls for dependencies as needed
WINNOWMAP = "winnowmap"
FALCONC = "falconc"
MERYL = "meryl"
MERFIN = "merfin"
BCFTOOLS = "bcftools"
JELLYFISH = "jellyfish"
GENOMESCOPE = "genomescope2"
SAMTOOLS = "samtools"
MERQURY_SH = "merqury.sh" 

# Default quality scores for read correctors
DEFAULT_CORRECTOR_Q = {
    "hifiasm": 40,
    "herro": 35,
    "flye": 30,
    "unknown": 30  # fallback
}

DEFAULT_KMER_SIZE = 31

# ---------------------------------------------------------------------
# Sub-command: make_path, validate_dependecies, log_versions, diagnostics, 
# logger command, and needs_kcov
# Helper commands to improve portability and reliability
# ---------------------------------------------------------------------

def make_path(base, *parts, ext=""):
    return f"{base}." + "_".join(parts) + ext

def validate_dependencies():
    for tool in [WINNOWMAP, FALCONC, MERYL, MERFIN, BCFTOOLS, JELLYFISH, GENOMESCOPE, SAMTOOLS, MERQURY_SH]:
        if shutil.which(tool) is None:
            logger.error(f"Required tool not found in PATH: {tool}")
            sys.exit(1)

def needs_kcov(args):
    return args.optimized and (
        not args.fitted_hist or
        (args.ploidy == "haploid" and not args.ideal_dpeak) or
        (args.ploidy == "diploid" and not args.ideal_hpeak)
    )

def conditional_run(outfile, resume_flag, step_desc, command, stream=False, **kwargs):
    """
    Run a command if needed, optionally streaming stdout directly to outfile.
    If stream=True, uses subprocess.Popen to write stdout to file without buffering.
    """
    if resume_flag and os.path.exists(outfile):
        logger.info("Resume: skipping %s (exists: %s)", step_desc, outfile)
        return None

    if stream:
        logger.info("Running: %s (streaming)", step_desc)
        with open(outfile, "w") as fh:
            proc = subprocess.Popen(
                command,
                stdout=fh,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )
            _, stderr = proc.communicate()
            if proc.returncode != 0:
                logger.error("[STDERR]\n%s", stderr.strip() or "No stderr")
                sys.exit(proc.returncode)
            logger.info("Completed: %s (streamed to %s)", step_desc, outfile)
        return None

    return run_command(command, description=step_desc, **kwargs)

def log_tool_versions(prefix):
    version_log = prefix + ".tool_versions.txt"
    tools = [WINNOWMAP, FALCONC, MERYL, MERFIN, BCFTOOLS, JELLYFISH, GENOMESCOPE, SAMTOOLS, MERQURY_SH]
    with open(version_log, "w") as f:
        for tool in tools:
            f.write(f"{tool} version:\n")
            try:
                # Run with a 10-second timeout to avoid hanging
                out = subprocess.run(
                    [tool, "--version"],
                    stdout=PIPE,
                    stderr=PIPE,
                    universal_newlines=True,
                    timeout=10
                )
                f.write(out.stdout.strip() + "\n\n")
            except subprocess.TimeoutExpired:
                f.write(f"Error: {tool} --version timed out after 10 seconds\n\n")
            except Exception as e:
                f.write(f"Error retrieving version: {e}\n\n")
    logger.info("Tool versions written to %s", version_log)

def diagnostic_mode():
    import os, shutil

    logger.info("[Diagnostics] Checking required tools and system status...\n")
    # Check external tools
    for tool in [WINNOWMAP, FALCONC, MERYL, MERFIN, BCFTOOLS, JELLYFISH, GENOMESCOPE, SAMTOOLS, MERQURY_SH]:
        path = shutil.which(tool)
        if path:
            logger.info("[OK] %s found at %s", tool, path)
        else:
            logger.info("[Missing] %s not found in PATH", tool)

    # Check available disk space on current working directory
    total, used, free = shutil.disk_usage(os.getcwd())
    logger.info(
        "Diagnostics: Disk space - Total: %.2f GB, Used: %.2f GB, Free: %.2f GB",
        total/1e9, used/1e9, free/1e9
    )

    # Check presence of key Python packages
    for pkg in ["pysam", "json", "argparse", "subprocess", "concurrent.futures"]:
        try:
            __import__(pkg)
            logger.info("[OK] Python package %s is installed", pkg)
        except ImportError:
            logger.warning("Missing Python package: %s", pkg)

    logger.info("\n[Diagnostics] Check complete.")

def setup_logging(log_file: str = None, quiet: bool = False):
    """
    Configure root logger:
     - Console: INFO (or WARNING if --quiet)
     - File: DEBUG if --log-file given
     - Timestamped entries
    """
    # Remove any existing handlers to avoid duplicate logs on multiple runs
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    fmt    = "%(asctime)s %(levelname)-8s %(message)s"
    datefmt= "%Y-%m-%d %H:%M:%S"
    handlers = []

    # Console handler
    console_level = logging.WARNING if quiet else logging.INFO
    ch = logging.StreamHandler()
    ch.setLevel(console_level)
    ch.setFormatter(logging.Formatter(fmt, datefmt))
    handlers.append(ch)

    # File handler, if requested
    if log_file:
        fh = logging.FileHandler(log_file, mode='w')
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(logging.Formatter(fmt, datefmt))
        handlers.append(fh)

    # Apply configuration
    logging.basicConfig(level=logging.DEBUG, handlers=handlers)
    
# ---------------------------------------------------------------------
# Sub-command: run_command
# A helper function that provides stderr, stdout, and runtimes to the user
# ---------------------------------------------------------------------

class CommandResult:
    def __init__(self, stdout, stderr, returncode, elapsed):
        self.stdout = stdout
        self.stderr = stderr
        self.returncode = returncode
        self.elapsed = elapsed

def run_command(cmd, description=None, use_shell=False, flush=True, logfile=None, max_log_lines=100):
    # Helper: shorten very long stdout/stderr for log display
    def _truncate(text, limit):
        if limit is None:
            return text
        lines = text.splitlines()
        if len(lines) <= limit:
            return text
        return "\n".join(lines[:limit] +
                         [f"... ({len(lines) - limit} more lines truncated)"])
    if description:
        logger.info("Running: %s", description)
    logger.debug("Command: %s", cmd if isinstance(cmd, str) else " ".join(cmd))

    start = time.time()
    try:
        result = subprocess.run(
            cmd,
            shell=use_shell,
            check=True,
            capture_output=True,
            text=True
        )
        elapsed = time.time() - start

        stdout = result.stdout.strip() or "No output"
        stderr = result.stderr.strip() or "No errors"

        logger.info("[STDOUT]\n%s", _truncate(stdout, max_log_lines))
        logger.info("[STDERR]\n%s", _truncate(stderr, max_log_lines))

        if logfile:
            with open(logfile, "a") as f:
                f.write(f"[{description}]\nCommand: {' '.join(cmd) if isinstance(cmd, list) else cmd}\n")
                f.write(f"STDOUT:\n{stdout}\nSTDERR:\n{stderr}\nElapsed: {elapsed:.2f} seconds\n\n")

        return CommandResult(stdout, stderr, result.returncode, elapsed)

    except subprocess.CalledProcessError as e:
        elapsed = time.time() - start
        stdout = e.stdout.strip() if e.stdout else "No output"
        stderr = e.stderr.strip() if e.stderr else "No errors"

        logger.error("[Fatal Error] %s (elapsed: %.2f s)", description if description else cmd, elapsed)
        logger.error("[STDOUT]\n%s", _truncate(stdout, max_log_lines))
        logger.error("[STDERR]\n%s", _truncate(stderr, max_log_lines))

        return CommandResult(stdout, stderr, e.returncode, elapsed)

# ---------------------------------------------------------------------
# Sub-command: computekcov
#   1) Jellyfish count/histo
#   2) GenomeScope2 with ploidy variable
# ---------------------------------------------------------------------

def sub_computekcov(args):
    if not os.path.exists(args.reads):
        logger.error(f"Reads file not found: {args.reads}")
        sys.exit(1)

    # If resume flag is enabled and outputs are already present, skip computation.
    if hasattr(args, 'resume') and args.resume:
        hist_file = f"{args.prefix}.histo"
        gs_lookup = os.path.join(f"{args.prefix}.genomescope_out", "lookup_table.txt")
        if os.path.exists(hist_file) and os.path.exists(gs_lookup):
            logger.info("compute_kcov: Resume enabled; skipping compute-kcov steps as results exist.")
            params_txt = os.path.join(f"{args.prefix}.genomescope_out", "params.txt")
            discovered_kcov = None
            if os.path.exists(params_txt):
                with open(params_txt) as f:
                    for line in f:
                        match = re.match(r"kcov:\s*(\S+)", line)
                        if match:
                            discovered_kcov = match.group(1)
                            break
            fitted_hist = gs_lookup
            logger.info("compute_kcov: GenomeScope2 finished (resumed).")
            logger.info("compute_kcov: Ideal KCov = %s", discovered_kcov if discovered_kcov else "Unknown")
            logger.info("compute_kcov: Fitted histogram location = %s", fitted_hist)
            logger.info("Using previously computed --fitted_hist %s and --peak %s", fitted_hist, discovered_kcov)
            return discovered_kcov, fitted_hist

    # Continue with the full computation if resume isn't enabled or output files aren't found.
    logger.info("compute_kcov: Starting computation of k-mer coverage using Jellyfish and GenomeScope2")

    prefix = args.prefix
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

    jf_out = f"{prefix}.jf"
    jf_histo = f"{prefix}.histo"
    gs_outdir = f"{prefix}.genomescope_out"

    logger.info("compute_kcov: Step 1: Running Jellyfish K-mer counting (this may take a while)...")
    conditional_run(
        outfile=jf_out,
        resume_flag=(hasattr(args, "resume") and args.resume),
        step_desc="[computekcov step 1] Jellyfish K-mer counting",
        command=[
            JELLYFISH, "count",
            "-C",
            "-m", str(args.kmer_size),
            "-s", str(args.jellyfish_hash_size),
            "-t", str(args.threads),
            "-o", jf_out,
            args.reads
        ]
    )

    conditional_run(
        outfile=jf_histo,
        resume_flag=not (args.resume and os.path.exists(jf_histo)),
        step_desc="[computekcov step 2] Jellyfish histogram generation",
        command=[
            JELLYFISH, "histo",
            "-t", str(args.threads),
            jf_out,
            "-o", jf_histo
        ]
    )

    logger.info("compute_kcov: Step 3: Running GenomeScope2 for histogram fitting and coverage estimation...")
    gs_cmd = [
        GENOMESCOPE,
        "-i", jf_histo,
        "-o", gs_outdir,
        "-k", str(args.kmer_size),
        "--fitted_hist"
    ]
    if args.ploidy == "haploid":
        gs_cmd += ["--ploidy", "1"]
    else:
        gs_cmd += ["--ploidy", "2"]

    # Step 4: Run GenomeScope2 and capture kcov from stdout
    result = conditional_run(
        outfile=os.path.join(gs_outdir, "lookup_table.txt"),
        resume_flag=args.resume,
        step_desc="[computekcov step 4] GenomeScope2 fitted histogram file creation and coverage calculation for pipeline",
        command=gs_cmd
    )
    # Try to parse kcov from GenomeScope2 stdout (e.g., "kcov:115")
    discovered_kcov = None
    if result and result.stdout:
        match = re.search(r"kcov:([0-9]+(?:\\.[0-9]+)?)", result.stdout)
        if match:
            discovered_kcov = match.group(1)
    # Fallback: read from params.txt if stdout did not contain kcov
    params_txt = os.path.join(gs_outdir, "params.txt")
    if discovered_kcov is None and os.path.exists(params_txt):
        with open(params_txt) as f:
            for line in f:
                match = re.match(r"kcov:\\s*(\\S+)", line)
                if match:
                    discovered_kcov = match.group(1)
                    break

    fitted_hist = os.path.join(gs_outdir, "lookup_table.txt")

    logger.info("compute_kcov: GenomeScope2 finished.")
    logger.info("compute_kcov: Ideal KCov = %s", discovered_kcov if discovered_kcov else "Unknown")
    logger.info("compute_kcov: Fitted histogram location = %s", fitted_hist)
    logger.info("If running compute_kcov separately, pass --fitted_hist %s and --peak %s to the polish step. Otherwise, the automated pipeline will incorporate them.", fitted_hist, discovered_kcov)
    
    # Persist coverage parameters for reliable resume under original prefix
    meta = {"kcov": discovered_kcov, "fitted_hist": fitted_hist}
    # strip trailing ".kcov" if present
    meta_prefix = args.prefix[:-5] if args.prefix.endswith(".kcov") else args.prefix
    meta_file = f"{meta_prefix}.kcov.json"
    with open(meta_file, "w") as mf:
        json.dump(meta, mf)
    return discovered_kcov, fitted_hist

# ---------------------------------------------------------------------
# Sub-command: run_merfin_eval and run_merqury_eval
# Helper functions for QV evaluations that are called after each step
# ---------------------------------------------------------------------

def run_merfin_eval(consensus, readmers, optimized, ideal_kcov, fitted_hist_location):
    """
    Returns only the QV and completeness lines from Merfin evaluation.
    """
    # Build histogram (“hist”) and completeness (“completeness”) commands
    hist_out = f"{consensus}.merfin_hist.txt"
    comp_out = f"{consensus}.merfin_completeness.txt"

    # Histogram command: sequence → readmers → peak → output
    peak_val = str(ideal_kcov) if optimized and ideal_kcov else "106.7"
    hist_cmd = [
        MERFIN, "-hist",
        "-sequence", consensus,
        "-readmers", readmers,
        "-peak", peak_val,
        "-output", hist_out
    ]
    # Include probability table if available
    if optimized and fitted_hist_location:
        hist_cmd += ["-prob", fitted_hist_location]

    # Completeness command: same ordering
    comp_cmd = [
        MERFIN, "-completeness",
        "-sequence", consensus,
        "-readmers", readmers,
        "-peak", peak_val,
        "-output", comp_out
    ]
    if optimized and fitted_hist_location:
        comp_cmd += ["-prob", fitted_hist_location]

    import subprocess

    # 1) Histogram calculation (Merfin writes output to hist_out file)
    hist_proc = subprocess.run(hist_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    if hist_proc.returncode != 0:
        logger.error("[Merfin Eval] Histogram calculation failed: %s", hist_proc.stderr.strip())
        hist_text = ""
    else:
        # Merfin ≥1.1 writes summary to stdout; fall back to stderr if empty
        hist_text = (hist_proc.stdout or hist_proc.stderr).strip()

    # 2) Completeness calculation (Merfin writes output to comp_out file)
    comp_proc = subprocess.run(comp_cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True)
    if comp_proc.returncode != 0:
        logger.error("[Merfin Eval] Completeness calculation failed: %s", comp_proc.stderr.strip())
        comp_text = ""
    else:
        # Merfin ≥1.1 writes summary to stdout; fall back to stderr if empty
        comp_text = (comp_proc.stdout or comp_proc.stderr).strip()

    # Filter Merfin outputs to only QV and completeness metrics
    hist_lines = [l for l in hist_text.splitlines() if re.search(r'(QV|Mean)', l)]
    comp_lines = [l for l in comp_text.splitlines() if re.search(r'Completeness', l)]
    filtered = "\n".join(hist_lines + comp_lines)
    return filtered

def run_merqury_eval(consensus, readmers, resume=False):
    """
    Returns only the QV and completeness numbers from Merqury output files.
    """
    iter_dir = os.path.dirname(os.path.abspath(consensus))
    base = Path(consensus).stem
    outprefix = os.path.join(iter_dir, f"{base}_merqury")
    # Always launch via bash so we don't rely on merqury.sh having +x or being in PATH
    cmd = ["bash", MERQURY_SH, readmers, consensus, outprefix]

    # Run Merqury (it writes .qv and .completeness.stats files)
    conditional_run(
        outfile=outprefix + ".qv",
        resume_flag=resume,
        step_desc="[Merqury Eval] Running merqury.sh",
        command=cmd
    )
    # If the .qv file is still missing after running, treat as failure
    if not os.path.exists(outprefix + ".qv"):
        logger.error("Merqury failed: %s did not produce %s", MERQURY_SH, outprefix + ".qv")
        return "Merqury QV: ERROR\nMerqury Completeness: ERROR"

    # Extract QV from the .qv file (4th column)
    qv_file = outprefix + ".qv"
    qv_value = None
    if os.path.exists(qv_file):
        with open(qv_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 4:
                    qv_value = parts[3]
                    break
    else:
        logger.warning("Merqury QV file not found: %s", qv_file)

    # Extract completeness from the .completeness.stats file (5th column)
    comp_file = outprefix + ".completeness.stats"
    comp_value = None
    if os.path.exists(comp_file):
        with open(comp_file) as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 5:
                    comp_value = parts[4]
                    break
    else:
        logger.warning("Merqury completeness file not found: %s", comp_file)

    # Build and return a simple two-line summary
    summary_lines = []
    if qv_value is not None:
        summary_lines.append(f"Merqury QV: {qv_value}")
    else:
        summary_lines.append("Merqury QV: N/A")
    if comp_value is not None:
        summary_lines.append(f"Merqury Completeness: {comp_value}")
    else:
        summary_lines.append("Merqury Completeness: N/A")
    return "\n".join(summary_lines)

# ---------------------------------------------------------------------
# Sub-command: compute_readmers_db
# Helper function to automatically calculate a Meryl read DB
# ---------------------------------------------------------------------

def compute_readmers_db(reads, output_readmers, k=31, threads=8):
    cmd = [
        MERYL, f"k={k}",
        f"threads={threads}",
        "count", reads, "output", output_readmers
    ]
    # Safety check: if the intended meryl output path already exists as a
    # *file* (from a prior failed run), delete it so that meryl can create
    # the required directory.
    if os.path.exists(output_readmers) and not os.path.isdir(output_readmers):
        logger.warning("%s exists as a file; removing it so meryl can create a directory.", output_readmers)
        os.remove(output_readmers)
    conditional_run(
        outfile=output_readmers,
        resume_flag=True,  # Or pass a flag if available
        step_desc="[compute_readmers_db] Computing readmers database with Meryl",
        command=cmd,
        max_log_lines=0
    )

# ---------------------------------------------------------------------
# Sub-command: fasta_to_fastq
# A fasta to fastq converter for DeepVariant, based on read correctors
# ---------------------------------------------------------------------

def fasta_to_fastq(fasta_file, fastq_file, quality_score=40):
    qual_char = chr(quality_score + 33)  # FASTQ ASCII encoding
    with open(fasta_file) as fa, open(fastq_file, "w") as fq:
        name, seq = None, []
        for line in fa:
            line = line.strip()
            if line.startswith(">"):
                if name and seq:
                    fq.write(f"@{name}\n{''.join(seq)}\n+\n{qual_char * len(''.join(seq))}\n")
                elif name and not seq:
                    logger.warning("FASTA record %s has no sequence data and was skipped.", name)
                name = line[1:]
                seq = []
            else:
                seq.append(line)
        if name and seq:
            fq.write(f"@{name}\n{''.join(seq)}\n+\n{qual_char * len(''.join(seq))}\n")
        elif name and not seq:
            logger.warning("FASTA record %s has no sequence data and was skipped.", name)

# ---------------------------------------------------------------------
#  Sub-command: run_polishing_iteration
#     1) Meryl DB of repetitive k=15 from draft_fasta
#     2) winnowmap -> samtools view -> samtools sort
#     3) falconc filter => SAM
#     4) VCF Generation
#     5) Merfin => produce .polish.vcf
#     6) bcftools => final FASTA
#     7) Parallel QV: Merfin + Merqury
# Returns: (paths.consensus, merfin_output, merqury_output)
# ---------------------------------------------------------------------

class PolishIteration:
    def __init__(self, folder, prefix):
        self.folder = folder
        self.prefix = prefix

        # Meryl DB lives directly under folder
        self.repet_db = os.path.join(self.folder, "merylDB")

        # All other files should be in the folder as well
        self.repet_k15_txt = os.path.join(self.folder, f"{self.prefix}.repetitive_k15.txt")
        self.sorted_bam     = os.path.join(self.folder, f"{self.prefix}.winnowmap.sorted.bam")
        self.sam            = os.path.join(self.folder, f"{self.prefix}.falconc.sam")
        self.bam            = os.path.join(self.folder, f"{self.prefix}.falconc.sorted.bam")

        # DeepVariant output directory under folder
        self.dv_outdir      = os.path.join(self.folder, "deepvariant_results")
        os.makedirs(self.dv_outdir, exist_ok=True)
        self.dv_vcf         = os.path.join(self.dv_outdir, "deepvariant_output.vcf.gz")
        self.dv_gvcf        = os.path.join(self.dv_outdir, "deepvariant_output.g.vcf.gz")

        # Output Meryl DB from DeepVariant
        self.meryl_db       = os.path.join(self.folder, f"{self.prefix}.dv.meryl")

        # Merfin base and its derived files
        self.merfin_base    = os.path.join(self.folder, f"{self.prefix}.dv.merfin")
        self.merfin_vcf     = os.path.join(self.folder, f"{self.prefix}.dv.merfin.polish.vcf")
        self.merfin_vcfgz   = os.path.join(self.folder, f"{self.prefix}.dv.merfin.polish.vcfgz")

        # Final consensus FASTA
        self.consensus      = os.path.join(self.folder, f"{self.prefix}.consensus.fasta")

def run_polish_iteration(
    iter_folder,
    iteration_id,
    draft_fasta,
    reads,
    readmers,
    k_mer_size,
    num_threads,
    winnowmap_seq_type,
    deepseq_type,
    singularity_sif,
    cleanup=False,
    optimized=False,
    ideal_kcov=None,
    fitted_hist_location=None,
    resume=False,
    resume_from=0,
    meryl_memory=50,
    dv_tmp_dir=None
    ):

    os.makedirs(iter_folder, exist_ok=True)
    # Debug: confirm received coverage parameters
    logger.debug(
        "run_polish_iteration: ideal_kcov=%s, fitted_hist_location=%s",
        ideal_kcov, fitted_hist_location
    )

    logger.info("Polishing pipeline start: %s", time.ctime())

    # keep folder and filename-prefix separate:
    prefix = f"iter_{iteration_id}"
    paths = PolishIteration(iter_folder, prefix)

    # 1) Meryl DB of repetitive k=15 from draft_fasta

    # Ensure the repetitive k‑mer output path is available as a directory.
    if os.path.exists(paths.repet_db) and not os.path.isdir(paths.repet_db):
        logger.warning("%s exists as a file; removing it so meryl can create a directory.", paths.repet_db)
        os.remove(paths.repet_db)

    with open(os.devnull, 'w') as devnull:
        with redirect_stdout(devnull), redirect_stderr(devnull):
            conditional_run(
                outfile=paths.repet_db,
                resume_flag=(resume and resume_from <= 1),
                step_desc=f"[Iter {iteration_id} Step 1] Generation of Repetitive 15-mers from Genome",
                command=[
                    MERYL, "k=15", "count", f"threads={num_threads}", f"memory={meryl_memory}",
                    "output", paths.repet_db, draft_fasta
                ],
                max_log_lines=0
            )

    if not (resume and os.path.exists(paths.repet_k15_txt)):
        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                result = conditional_run(
                    outfile=paths.repet_k15_txt,
                    resume_flag=(resume and resume_from <= 1),
                    step_desc=f"[Iter {iteration_id} Step 1B] Filter Repetitive 15-mers",
                    command=[
                        MERYL, "print", "greater-than", "distinct=0.9998", paths.repet_db],
                    max_log_lines=0
                )
                if result:
                    with open(paths.repet_k15_txt, "w") as fh:
                        fh.write(result.stdout)

    # 2) winnowmap -> samtools view -> samtools sort
    # define the sorted BAM inside the iteration folder
    winnowmap_sorted_bam = os.path.join(iter_folder, f"{prefix}.winnowmap.sorted.bam")

    # Step 2A: Winnowmap → SAM
    sam_file = os.path.join(iter_folder, f"{prefix}.winnowmap.sam")
    
    conditional_run(
        outfile=sam_file,
        resume_flag=(resume and resume_from <= 2),
        step_desc=f"[Iter {iteration_id} Step 2A] Winnowmap → SAM",
        command=[
            WINNOWMAP, "-k", "15",
            "-W", paths.repet_k15_txt,
            "-t", str(num_threads),
            "-ax", f"map-{winnowmap_seq_type}",
            "--MD", draft_fasta,
            reads
        ],
        stream=True
    )

    # Step 2B: SAM → unsorted BAM
    unsorted_bam = os.path.join(iter_folder, f"{prefix}.winnowmap.unsorted.bam")
    conditional_run(
        outfile=unsorted_bam,
        resume_flag=(resume and resume_from <= 2),
        step_desc=f"[Iter {iteration_id} Step 2B] Samtools view",
        command=[
            SAMTOOLS, "view",
            "--threads", str(num_threads),
            "-hb",
            "-T", draft_fasta,
            sam_file
        ],
        stream=True
    )

    # Step 2C: sort into final BAM
    conditional_run(
        outfile=winnowmap_sorted_bam,
        resume_flag=(resume and resume_from <= 2),
        step_desc=f"[Iter {iteration_id} Step 2C] Samtools sort",
        command=[
            SAMTOOLS, "sort",
            "--threads", str(num_threads),
            "-o", winnowmap_sorted_bam,
            unsorted_bam
        ],
        stream=True
    )

    # 3) falconc filter => SAM
    conditional_run(
        outfile=paths.sam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3] FalconC Filter of Alignment",
        command=[
            FALCONC, "bam-filter-clipped",
            "-t",
            "-F", "0x104",
            "--input-fn", winnowmap_sorted_bam,
            "--output-fn", paths.sam,
            "--output-count-fn", f"{paths.sam}.filtered_aln_count.txt"
        ],
        stream=True
    )

    # Temporary unsorted BAM produced directly from FalconC‑filtered SAM
    falconc_unsorted_bam = os.path.join(iter_folder, f"{prefix}.falconc.unsorted.bam")
    # Convert FalconC SAM → unsorted BAM
    conditional_run(
        outfile=falconc_unsorted_bam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Convert] Samtools view FalconC SAM → BAM (unsorted)",
        command=[
            SAMTOOLS, "view", "-@", str(num_threads), "-bh",
            "-o", falconc_unsorted_bam, paths.sam
        ],
        stream=True
    )

    # Sort the FalconC BAM
    conditional_run(
        outfile=paths.bam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Sort] Samtools sort FalconC BAM",
        command=[
            SAMTOOLS, "sort", "-@", str(num_threads),
            "-o", paths.bam, falconc_unsorted_bam
        ],
        stream=True
    )

    conditional_run(
        outfile=paths.bam + ".csi",
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Index] Samtools index BAM",
        command=[SAMTOOLS, "index", "-c", paths.bam]
    )

    # 4) DeepVariant-based variant calling
    # Corruption checks for BAM and VCF files
    if resume and resume_from <= 4 and os.path.exists(paths.bam):
        try:
            with pysam.AlignmentFile(paths.bam, "rb") as bamfile:
                bamfile.check_index()
        except Exception as e:
            logger.error("BAM file %s appears to be corrupted or unindexed: %s", paths.bam, e)
            sys.exit(1)
    
    if resume and resume_from <= 4 and os.path.exists(paths.dv_vcf):
        try:
            with pysam.VariantFile(paths.dv_vcf, "r"):
                pass
        except Exception as e:
            logger.error("VCF file %s is invalid or unreadable: %s", paths.dv_vcf, e)
            sys.exit(1)

    # choose and create a per-iteration tmp dir for DeepVariant
    dv_tmp = dv_tmp_dir or os.path.join(iter_folder, "dv_tmp")
    # make absolute so Singularity bind destination is valid
    dv_tmp = os.path.abspath(dv_tmp)
    os.makedirs(dv_tmp, exist_ok=True)

    os.environ["TMPDIR"]             = dv_tmp
    os.environ["TMP"]                = dv_tmp
    os.environ["DEEPVARIANT_TMPDIR"] = dv_tmp

    deepvariant_cmd = [
        "singularity", "exec", "--nv",
        "-B", f"{dv_tmp}:{dv_tmp}",
        singularity_sif,
        "run_deepvariant",
        f"--model_type={deepseq_type}",
        f"--ref={draft_fasta}",
        f"--reads={paths.bam}",
        f"--output_vcf={paths.dv_vcf}",
        f"--output_gvcf={paths.dv_gvcf}",
        f"--intermediate_results_dir={dv_tmp}",
        f"--num_shards={num_threads}"
    ]

    conditional_run(
        outfile=paths.dv_vcf,
        resume_flag=(resume and resume_from <= 4),
        step_desc=f"[Iter {iteration_id} Step 4] DeepVariant",
        command=deepvariant_cmd
    )

    if not os.path.exists(paths.dv_vcf):
        logger.error("DeepVariant did not produce the expected output VCF. Check logs and command.")
        sys.exit(1)

    # 5) Meryl DB from DeepVariant consensus
    if os.path.exists(paths.meryl_db) and not os.path.isdir(paths.meryl_db):
        logger.warning("%s exists as a file; removing it so meryl can create a directory.", paths.meryl_db)
        os.remove(paths.meryl_db)
    conditional_run(
        outfile=paths.meryl_db,
        resume_flag=(resume and resume_from <= 5),
        step_desc=f"[Iter {iteration_id} Step 5] MerylDB from Genome Draft",
        command=[
            MERYL,
            "-Q", f"k={k_mer_size}",
            f"threads={num_threads}",
            "memory=50",
            "count",
            "output", paths.meryl_db,
            draft_fasta
        ],
        max_log_lines=0
    )

    # 6) Merfin => produce .polish.vcf using DeepVariant Consensus
    # build merfin_cmd once, then optionally tweak or skip
    merfin_cmd = [
        MERFIN, "-polish",
        "-vcf", paths.dv_vcf,
        "-sequence", draft_fasta,
        "-seqmers", paths.meryl_db,
        "-readmers", readmers,
        "-output", paths.merfin_base,
        "-threads", str(num_threads)
    ]
    if optimized and ideal_kcov and fitted_hist_location:
        merfin_cmd += ["-peak", str(ideal_kcov), "-prob", fitted_hist_location]
        
    conditional_run(
        outfile=paths.merfin_vcf,
        resume_flag=(resume and resume_from <= 5),
        step_desc=f"[Iter {iteration_id} Step 5] Merfin Polish",
        command=merfin_cmd
    )


    # Step 6A: Compress the Merfin-polished VCF to .vcfgz
    conditional_run(
        outfile=paths.merfin_vcfgz,
        resume_flag=(resume and resume_from <= 6),
        step_desc=f"[Iter {iteration_id} Step 6A] BCFtools compress Merfin VCF",
        command=[
            BCFTOOLS, "view",
            "-Oz",
            paths.merfin_vcf,
            "-o", paths.merfin_vcfgz
        ]
    )

    # Step 6B: Index the compressed VCF
    conditional_run(
        outfile=paths.merfin_vcfgz + ".tbi",
        resume_flag=(resume and resume_from <= 6),
        step_desc=f"[Iter {iteration_id} Step 6B] BCFtools index Merfin VCF",
        command=[
            BCFTOOLS, "index",
            paths.merfin_vcfgz
        ]
    )

    # Step 6C: Final consensus FASTA
    conditional_run(
        outfile=paths.consensus,
        resume_flag=(resume and resume_from <= 6),
        step_desc=f"[Iter {iteration_id} Step 6 Final Consensus] BCFtools consensus Merfin-polished VCF → Final FASTA",
        command=[
            BCFTOOLS, "consensus",
            paths.merfin_vcfgz,
            "-f", draft_fasta,
            "-H", "1",
            "-o", paths.consensus
        ]
    )
    logger.info("Completed iteration %d: consensus at %s", iteration_id, paths.consensus)

    if resume_from <= 6 and resume and cleanup:
        logger.info("Cleanup: Removing intermediate files to conserve space")
        for file in [paths.sam, paths.merfin_vcf]:
            if os.path.exists(file):
                os.remove(file)
                logger.info("Cleanup: Removed %s", file)


    # 8) Parallel QV: Merfin + Merqury
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        logger.info("Running QV Module: QV and completeness estimation via Merfin and Merqury")
        fut_merfin = executor.submit(
            run_merfin_eval,
            paths.consensus,
            readmers,
            optimized,
            ideal_kcov,
            fitted_hist_location
        )
        fut_merqury = executor.submit(
            run_merqury_eval,
            paths.consensus,
            readmers,
            resume=resume
        )
        merfin_output = fut_merfin.result()
        merqury_output = fut_merqury.result()

    logger.info("Polishing pipeline end: %s", time.ctime())

    return paths.consensus, merfin_output, merqury_output

# ---------------------------------------------------------------------
# Sub-command: sub_polish (collection of parameters and execution of pipeline)
# ---------------------------------------------------------------------

def sub_polish(args):
    # Auto-reload kcov from JSON when resume is requested
    kcov_meta = f"{args.prefix}.kcov.json"
    if args.resume and os.path.exists(kcov_meta):
        with open(kcov_meta) as mf:
            meta = json.load(mf)
        args.fitted_hist = meta.get("fitted_hist")
        if args.ploidy == "haploid":
            args.ideal_dpeak = meta.get("kcov")
        else:
            args.ideal_hpeak = meta.get("kcov")
        logger.info(
            "Resume: reloaded kcov=%s, fitted_hist=%s from %s",
            meta.get("kcov"), args.fitted_hist, kcov_meta
        )
    
    if args.resume and args.resume_from > 0:
        logger.warning("Both --resume and --resume-from set. Using --resume-from as priority.")
    
    if args.config_file:
        with open(args.config_file) as cf:
            config_data = json.load(cf)
        if not args.reads and "reads" in config_data:
            args.reads = config_data["reads"]
        if not args.readmers and "readmers" in config_data:
            args.readmers = config_data["readmers"]
        if args.optimized:
            # carry over legacy keys from config, if present
            if not args.ideal_dpeak and "peak_haploid" in config_data:
                args.ideal_dpeak = config_data["peak_haploid"]
            if not args.ideal_hpeak and "peak_diploid" in config_data:
                args.ideal_hpeak = config_data["peak_diploid"]
        # Handle fitted histogram from config
        if not args.fitted_hist and "fitted_hist_location" in config_data:
            args.fitted_hist = config_data["fitted_hist_location"]

    if args.resume and args.fitted_hist and (
        (args.ploidy == "haploid" and args.ideal_dpeak) or
        (args.ploidy == "diploid" and args.ideal_hpeak)
    ):
        # Resume mode with precomputed coverage parameters
        discovered_kcov = args.ideal_dpeak if args.ploidy == "haploid" else args.ideal_hpeak
        fitted_hist_location = args.fitted_hist
        logger.info(
            "Resume: using existing coverage parameters with --ideal_peak=%s and --fitted_hist=%s",
            discovered_kcov,
            fitted_hist_location
        )
        # Restore chosen_peak for downstream polishing
        chosen_peak = discovered_kcov

    if args.reads.lower().endswith(".fasta") or args.reads.lower().endswith(".fa"):
        if not args.read_corrector:
            logger.warning("FASTA reads detected, but no corrector specified. Defaulting to 'unknown'. DO NOT use uncorrected FASTA reads!")
            args.read_corrector = "unknown"
        quality_score = DEFAULT_CORRECTOR_Q.get(args.read_corrector, 30)
        converted_fastq = args.prefix + ".converted_reads.fastq"
        logger.info(
            "Corrected Reads Detected: converting FASTA to FASTQ with quality Q=%d for corrector: %s",
            quality_score,
            args.read_corrector
        )
        fasta_to_fastq(args.reads, converted_fastq, quality_score=quality_score)
        args.reads = converted_fastq

    if not args.draft or not os.path.exists(args.draft):
        logger.error("Missing or invalid --draft")
        sys.exit(1)
    if not args.reads or not os.path.exists(args.reads):
        logger.error("Missing or invalid --reads")
        sys.exit(1)
    if not args.readmers or not os.path.exists(args.readmers):
        computed_readmers = args.prefix + ".readmers.meryl"
        logger.info("Running Module: Meryl Started: Readmers not provided")
        with open(os.devnull, 'w') as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                compute_readmers_db(args.reads, computed_readmers, k=31, threads=args.threads)
        args.readmers = computed_readmers

    if not args.seq_type:  # means user didn't specify --seq_type
        if args.deepseq_type == "PACBIO":
            args.seq_type = "pb"
            logger.info("Auto-Setting: setting winnowmap seq_type to 'pb' because deepseq_type = 'PACBIO'")
        elif args.deepseq_type == "ONT_R104":
            args.seq_type = "ont"
            logger.info("Auto-Setting: setting winnowmap seq_type to 'ont' because deepseq_type = 'ONT_R104'")
        else:
            # If deepseq_type is something else (WGS, WES, HYBRID_PACBIO_ILLUMINA, etc.),
            # we can't guess how winnowmap should run => user must specify.
            logger.error(
                "deepseq_type=%s isn't automatically mapped. Please specify --seq_type for winnowmap (e.g. pb or ont).",
                args.deepseq_type
            )
            sys.exit(1)

    if needs_kcov(args):
        # Optimized run: compute kcov
        logger.info("Optimized Run: running internal compute-kcov step to determine coverage parameters")
        cov_prefix = args.prefix + ".kcov"
        class CKArgs: pass
        ckargs = CKArgs()
        ckargs.resume = args.resume
        ckargs.jellyfish_hash_size = getattr(args, "jellyfish_hash_size", 1000000000)
        ckargs.reads = args.reads
        ckargs.prefix = cov_prefix
        ckargs.kmer_size = args.kmer_size
        ckargs.threads = args.threads
        ckargs.ploidy = args.ploidy
        discovered_kcov, fitted_hist_location = sub_computekcov(ckargs)
        args.fitted_hist = fitted_hist_location
        if args.ploidy == "haploid":
            args.ideal_dpeak = discovered_kcov
            chosen_peak = discovered_kcov
        else:
            args.ideal_hpeak = discovered_kcov
            chosen_peak = discovered_kcov
    elif args.optimized:
        # Resume mode or parameters already set
        if args.ploidy == "haploid":
            chosen_peak = args.ideal_dpeak
            if chosen_peak is None:
                logger.error("ploidy=haploid but no --ideal_dpeak set.")
                sys.exit(1)
        else:
            chosen_peak = args.ideal_hpeak
            if chosen_peak is None:
                logger.error("ploidy=diploid but no --ideal_hpeak set.")
                sys.exit(1)
    else:
        chosen_peak = None

    # Warn/error: --optimized with --resume requires kcov metadata file
    if args.optimized and args.resume:
        kcov_meta = f"{args.prefix}.kcov.json"
        if os.path.exists(kcov_meta):
            logger.warning(
                "Running with --optimized and --resume: loading kcov parameters from %s",
                kcov_meta
            )
        else:
            logger.error(
                "Missing kcov metadata file %s required for --optimized --resume",
                kcov_meta
            )
            sys.exit(1)

    logger.info(
        "Launching polishing with ideal_kcov=%s, fitted_hist=%s",
        (args.ideal_dpeak if args.ploidy=="haploid" else args.ideal_hpeak),
        args.fitted_hist
    )

    prefix = args.prefix
    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

    # iteration 0 => copy user draft
    iteration0 = f"{prefix}.iter_0.consensus.fasta"

    conditional_run(
    outfile=iteration0,
    resume_flag=(args.resume and args.resume_from <= 0),
    step_desc="[Iteration 0] Copy initial draft",
    command=["cp", args.draft, iteration0]
    )

    # QV data from each iteration
    all_qv_records = []

    old_consensus = os.path.abspath(iteration0)
    for i in range(args.iterations):
        iter_num = i + 1
        iter_folder = f"{prefix}_iteration_{iter_num}"
        os.makedirs(iter_folder, exist_ok=True)

        new_consensus, merfin_out, merqury_out = run_polish_iteration(
            iter_folder=iter_folder,
            iteration_id=iter_num,
            draft_fasta=old_consensus,
            reads=os.path.abspath(args.reads),
            readmers=os.path.abspath(args.readmers),
            k_mer_size=args.kmer_size,
            num_threads=args.threads,
            winnowmap_seq_type=args.seq_type,
            deepseq_type=args.deepseq_type,
            singularity_sif=args.singularity_sif,
            optimized=args.optimized,
            ideal_kcov=chosen_peak,
            fitted_hist_location=args.fitted_hist,
            resume=args.resume,
            resume_from=args.resume_from,
            meryl_memory=args.meryl_memory,
            dv_tmp_dir=args.dv_tmp_dir
        )

        # Copy final consensus to top-level
        final_cons_name = f"{prefix}.iter_{iter_num}.consensus.fasta"
        if not os.path.exists(final_cons_name):
            shutil.copyfile(new_consensus, final_cons_name)

        # Store iteration QV
        all_qv_records.append((iter_num, merfin_out, merqury_out))
        # Persist per‑iteration QC to disk
        iter_qc = os.path.join(iter_folder, "QC_summary.txt")
        with open(iter_qc, "w") as qcfh:
            qcfh.write("----- Merfin -----\n")
            qcfh.write(merfin_out + "\n")
            qcfh.write("----- Merqury -----\n")
            qcfh.write(merqury_out + "\n")
        logger.info("Saved QC summary to %s", iter_qc)
        old_consensus = os.path.abspath(final_cons_name)

    # Print final summary
    logger.info("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======")
    for (iter_num, merfin_text, merqury_text) in all_qv_records:
        logger.info("Iteration %d", iter_num)
        logger.info("----- Merfin -----")
        logger.info("%s", merfin_text)
        logger.info("----- Merqury -----")
        logger.info("%s", merqury_text)

    summary_file = f"{prefix}.QV_Completeness_summary.txt"
    with open(summary_file, "w") as fh:
        fh.write("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======\n\n")
        for (iter_num, merfin_text, merqury_text) in all_qv_records:
            fh.write(f"> Iteration {iter_num}\n")
            fh.write("----- Merfin -----\n")
            fh.write(merfin_text + "\n")
            fh.write("----- Merqury -----\n")
            fh.write(merqury_text + "\n")

    logger.info("Done: wrote final summary to %s", summary_file)

# ---------------------------------------------------------------------
# Helper code to pull up help for users
# ---------------------------------------------------------------------

class ShowHelpOnErrorArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        # Print the full help message
        self.print_help(sys.stderr)
        # Show the error message
        sys.stderr.write(f"\nerror: {message}\n")
        sys.exit(2)

# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    validate_dependencies()
    parser = ShowHelpOnErrorArgumentParser(
        description="Integrated T2T Polishing Pipeline with K-mer Coverage Optimization. Subcommands: polish, computekcov."
    )
    subparsers = parser.add_subparsers(dest="subcommand")

    # Sub-command: computekcov
    parser_ck = subparsers.add_parser("computekcov",
        help="Compute optimal K-coverage and a fitted K-mer histogram using Jellyfish + GenomeScope2.")
    parser_ck.add_argument("-r", "--reads", required=True,
        help="Path to reads (FASTA/FASTQ).")
    parser_ck.add_argument("-o", "--prefix", default="KCoverage",
        help="Prefix for Jellyfish/GenomeScope2 outputs (default: KCoverage)")
    parser_ck.add_argument("-k", "--kmer_size", type=int, default=21,
        help="K-mer size (default: 21, recommended by Jellyfish / GenomeScope)")
    parser_ck.add_argument("-t", "--threads", type=int, default=8,
        help="Number of threads (default: 8)")
    parser_ck.add_argument("--ploidy", choices=["haploid","diploid"], default="haploid",
        help="Genome ploidy for GenomeScope2. (default: haploid => ploidy = 1)")
    parser_ck.add_argument("--jellyfish-hash-size", type=int, default=1000000000,
        help="Hash size (number of k-mers) for Jellyfish counting (default: 1e9)")
    parser_ck.set_defaults(func=sub_computekcov)

    # Sub-command: required polish arguments
    parser_pol = subparsers.add_parser("polish",
        help="Run iterative polishing with DeepVariant-based variant calling, Merfin, etc. Use --optimized to automatically run computekcov and input the calculated values automatically.")
    parser_pol.add_argument("-d", "--draft",
        help="Draft assembly FASTA")
    parser_pol.add_argument("-r", "--reads",
        help="Reads in FASTQ or FASTA format, GZIP is currently not accepted")
    parser_pol.add_argument("--singularity_sif", required=True,
        help="Path to GPU-accelerated DeepVariant Singularity image.")
    parser_pol.add_argument("--deepseq_type",
        choices=["WGS", "WES", "PACBIO", "ONT_R104", "HYBRID_PACBIO_ILLUMINA"],
        default="PACBIO",
        help="DeepVariant model type for read input. (default: PACBIO)")
    
    # Sub-command: optional polish arguments 
    parser_pol.add_argument("-o", "--prefix", default="AutoPolisher",
        help="Output prefix (default: AutoPolisher)")
    parser_pol.add_argument("-t", "--threads", type=int, default=32,
        help="Number of threads. Equal to number of shards in DeepVariant. (default: 32)")
    parser_pol.add_argument("--optimized", action="store_true",
        help="Enable automatic coverage/hist calculation (Merfin -peak / -prob).")
    parser_pol.add_argument("--resume", action="store_true",
        help="Resume from a partially completed run. If set, computationally expensive steps can be skipped if their output already exists. Pipeline will automatically check for existing steps with the same prefix.")
    parser_pol.add_argument("--resume-from", type=int, default=0,
        help="Resume from a specific polishing step (0=all, 1=Meryl DB, 2=Winnowmap, 3=FalconC, 4=DeepVariant, 5=Merfin, 6=Consensus)")
    parser_pol.add_argument("-s", "--seq_type",
        help="Optionally set winnowmap read type (pb, ont). If not set, we may guess from --deepseq_type.")
    parser_pol.add_argument("--dv_tmp_dir", default=None,
        help="Directory to use for DeepVariant intermediate results (default: <iteration_folder>/dv_tmp)")
    parser_pol.add_argument("--read_corrector", choices=["hifiasm", "herro", "flye", "unknown"],
        help="Specify the read corrector used to assign default quality scores for FASTA input.")
    parser_pol.add_argument("-m", "--readmers", help="Meryl DB of read k-mers. Calculated automatically with --optimized.")
    parser_pol.add_argument("-k", "--kmer_size", type=int, default=DEFAULT_KMER_SIZE,
        help="K-mer size for Meryl. (default: 31, recommended by developers)")
    parser_pol.add_argument("-i", "--iterations", type=int, default=3,
        help="Number of polishing iterations (default: 3)")
    parser_pol.add_argument("--meryl-memory", type=int, default=50,
        help="Maximum memory (in GB) to hand off to Meryl (default: 50)")
    parser_pol.add_argument("--cleanup", action="store_true",
        help="If set, intermediate files (SAMs, uncompressed VCFs) are deleted after use to save disk space.")
    parser_pol.add_argument("--log-file",
        help="Path to write detailed log output (default: none—only console)")
    parser_pol.add_argument("--quiet", action="store_true",
        help="Suppress INFO‑level logs to console; only WARNING+ERROR will print")

    # Peak Settings
    parser_pol.add_argument("--ideal_dpeak", type=float,
                            help="Coverage peak to use for haploid assemblies (was --peak_haploid).")
    parser_pol.add_argument("--ideal_hpeak", type=float,
                            help="Coverage peak to use for diploid assemblies (was --peak_diploid).")
    parser_pol.add_argument("--fitted_hist", help="Fitted histogram location if not run as optimized")
    parser_pol.add_argument("--ploidy", choices=["haploid","diploid"], default="haploid",
                            help="Assembly ploidy (haploid or diploid) (default: haploid)")
    
    # Optional config file argument
    parser_pol.add_argument("--config_file",
                            help="Optional JSON config file with overrides to defaults, instead of command line entries.")
    
    parser_pol.set_defaults(func=sub_polish)

    parser_diag = subparsers.add_parser("diagnostics", help="Run a tool and environment diagnostic check")
    parser_diag.set_defaults(func=lambda args: diagnostic_mode())

    # Parse subcommand
    args = parser.parse_args()
    log_file = getattr(args, "log_file", None)
    quiet    = getattr(args, "quiet", False)
    setup_logging(log_file, quiet)
    logger = logging.getLogger(__name__)
    logger.info("Pipeline parameters: %s", vars(args))
    if hasattr(args, "threads"):
        if args.threads < 1:
            logger.error("Number of threads must be >= 1.")
            sys.exit(1)
        elif args.threads > 128:
            logger.warning("You requested more than 128 threads. This may exceed available system resources.")

    if not args.subcommand:
        parser.print_help()
        sys.exit(1)

    # Log run parameters
    try:
        param_log_path = args.prefix + ".run_parameters.json"
        params = vars(args).copy()
        params.pop("func", None)
        with open(param_log_path, "w") as f:
            json.dump(params, f, indent=4)
        logger.info("Saved run parameters to %s", param_log_path)
    except Exception as e:
        logger.warning("Could not save run parameters: %s", e)
    
    log_tool_versions(args.prefix)
    
    # Check input files
    critical_inputs = [args.draft, args.reads]
    for path in critical_inputs:
        if path and not os.path.exists(path):
            logger.error("Required input file not found: %s", path)
            sys.exit(1)
        elif path and os.path.getsize(path) == 0:
            logger.error("Input file %s exists but is empty.", path)
            sys.exit(1)
    
    # Execute subcommand and catch any unhandled exceptions
    try:
        args.func(args)
    except Exception as e:
        logger.exception("Pipeline failed with unhandled exception")
        sys.exit(1)

if __name__ == "__main__":
    main()

# T2T-Polish v4.0
# Any usage is subject to this software's license.
