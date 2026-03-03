#!/usr/bin/env python3
"""
T2T-Polish v4.1 — polish.py

Core polishing iteration and the ``sub_polish`` orchestrator.

Author: Patrick Grady
Anthropic Claude Opus 4.6 used for code formatting and cleanup assistance.
License: Public Domain (U.S. Government Work - NHGRI/NIH)
"""

from __future__ import annotations

import concurrent.futures
import json
import logging
import os
import shutil
import time
from dataclasses import dataclass

import pysam  # type: ignore

from t2t_polish.constants import (
    BCFTOOLS,
    DEFAULT_CORRECTOR_Q,
    DEFAULT_KMER_SIZE,
    FALCONC,
    MERFIN,
    MERYL,
    SAMTOOLS,
    WINNOWMAP,
)
from t2t_polish.evaluation import run_merfin_eval, run_merqury_eval
from t2t_polish.exceptions import InputValidationError, PipelineStepError
from t2t_polish.kcov import sub_computekcov
from t2t_polish.runner import conditional_run
from t2t_polish.utils import compute_readmers_db, fasta_to_fastq, needs_kcov

logger = logging.getLogger(__name__)


# ── Iteration paths ─────────────────────────────────────────────────


class PolishIteration:
    """Holds all per-iteration file paths for a single polishing round."""

    def __init__(self, folder: str, prefix: str) -> None:
        self.folder = folder
        self.prefix = prefix

        self.repet_db = os.path.join(folder, "merylDB")
        self.repet_k15_txt = os.path.join(folder, f"{prefix}.repetitive_k15.txt")
        self.sorted_bam = os.path.join(folder, f"{prefix}.winnowmap.sorted.bam")
        self.sam = os.path.join(folder, f"{prefix}.falconc.sam")
        self.bam = os.path.join(folder, f"{prefix}.falconc.sorted.bam")

        self.dv_outdir = os.path.join(folder, "deepvariant_results")
        os.makedirs(self.dv_outdir, exist_ok=True)
        self.dv_vcf = os.path.join(self.dv_outdir, "deepvariant_output.vcf.gz")
        self.dv_gvcf = os.path.join(self.dv_outdir, "deepvariant_output.g.vcf.gz")

        self.meryl_db = os.path.join(folder, f"{prefix}.dv.meryl")
        self.merfin_base = os.path.join(folder, f"{prefix}.dv.merfin")
        self.merfin_vcf = os.path.join(folder, f"{prefix}.dv.merfin.polish.vcf")
        self.merfin_vcfgz = os.path.join(folder, f"{prefix}.dv.merfin.polish.vcfgz")
        self.consensus = os.path.join(folder, f"{prefix}.consensus.fasta")


# ── Single polishing iteration ───────────────────────────────────────


@dataclass
class PolishConfig:
    """Typed configuration extracted from the CLI ``args`` namespace.

    Replaces the old pattern of mutating an argparse namespace with
    ``# type: ignore[attr-defined]`` annotations.  Built once at the
    top of :func:`sub_polish` and threaded through to helpers.
    """

    draft: str
    reads: str
    prefix: str = "AutoPolisher"
    threads: int = 32
    iterations: int = 3
    kmer_size: int = DEFAULT_KMER_SIZE
    meryl_memory: int = 50
    singularity_sif: str = ""
    deepseq_type: str = "PACBIO"
    seq_type: str = ""
    ploidy: str = "haploid"
    optimized: bool = False
    resume: bool = False
    resume_from: int = 0
    cleanup: bool = False
    dv_tmp_dir: str | None = None
    read_corrector: str | None = None
    readmers: str = ""
    ideal_dpeak: float | None = None
    ideal_hpeak: float | None = None
    fitted_hist: str | None = None
    config_file: str | None = None
    jellyfish_hash_size: int = 1_000_000_000
    json_summary: bool = False
    log_file: str | None = None
    quiet: bool = False

    @classmethod
    def from_args(cls, args: object) -> PolishConfig:
        """Build a :class:`PolishConfig` from an argparse-like namespace."""
        return cls(
            draft=getattr(args, "draft", "") or "",
            reads=getattr(args, "reads", "") or "",
            prefix=getattr(args, "prefix", "AutoPolisher"),
            threads=getattr(args, "threads", 32),
            iterations=getattr(args, "iterations", 3),
            kmer_size=getattr(args, "kmer_size", DEFAULT_KMER_SIZE),
            meryl_memory=getattr(args, "meryl_memory", 50),
            singularity_sif=getattr(args, "singularity_sif", ""),
            deepseq_type=getattr(args, "deepseq_type", "PACBIO"),
            seq_type=getattr(args, "seq_type", "") or "",
            ploidy=getattr(args, "ploidy", "haploid"),
            optimized=getattr(args, "optimized", False),
            resume=getattr(args, "resume", False),
            resume_from=getattr(args, "resume_from", 0),
            cleanup=getattr(args, "cleanup", False),
            dv_tmp_dir=getattr(args, "dv_tmp_dir", None),
            read_corrector=getattr(args, "read_corrector", None),
            readmers=getattr(args, "readmers", "") or "",
            ideal_dpeak=getattr(args, "ideal_dpeak", None),
            ideal_hpeak=getattr(args, "ideal_hpeak", None),
            fitted_hist=getattr(args, "fitted_hist", None),
            config_file=getattr(args, "config_file", None),
            jellyfish_hash_size=getattr(args, "jellyfish_hash_size", 1_000_000_000),
            json_summary=getattr(args, "json_summary", False),
            log_file=getattr(args, "log_file", None),
            quiet=getattr(args, "quiet", False),
        )


def run_polish_iteration(
    iter_folder: str,
    iteration_id: int,
    draft_fasta: str,
    reads: str,
    readmers: str,
    k_mer_size: int,
    num_threads: int,
    winnowmap_seq_type: str,
    deepseq_type: str,
    singularity_sif: str,
    cleanup: bool = False,
    optimized: bool = False,
    ideal_kcov: float | str | None = None,
    fitted_hist_location: str | None = None,
    resume: bool = False,
    resume_from: int = 0,
    meryl_memory: int = 50,
    dv_tmp_dir: str | None = None,
    dry_run: bool = False,
) -> tuple[str, str, str]:
    """Execute one full polish iteration and return
    ``(consensus_path, merfin_summary, merqury_summary)``.
    """
    os.makedirs(iter_folder, exist_ok=True)
    logger.debug(
        "run_polish_iteration: ideal_kcov=%s, fitted_hist_location=%s",
        ideal_kcov,
        fitted_hist_location,
    )
    logger.info("Polishing pipeline start: %s", time.ctime())

    prefix = f"iter_{iteration_id}"
    paths = PolishIteration(iter_folder, prefix)

    # ── Step 1 — Meryl repetitive k=15 from draft ───────────────────
    if os.path.exists(paths.repet_db) and not os.path.isdir(paths.repet_db):
        logger.warning(
            "%s exists as a file; removing it so meryl can create a directory.",
            paths.repet_db,
        )
        os.remove(paths.repet_db)

    conditional_run(
        outfile=paths.repet_db,
        resume_flag=(resume and resume_from <= 1),
        step_desc=f"[Iter {iteration_id} Step 1] Generation of Repetitive 15-mers from Genome",
        command=[
            MERYL, "k=15", "count", f"threads={num_threads}",
            f"memory={meryl_memory}", "output", paths.repet_db, draft_fasta,
        ],
        max_log_lines=0,
    )

    if not (resume and os.path.exists(paths.repet_k15_txt)):
        result = conditional_run(
            outfile=paths.repet_k15_txt,
            resume_flag=(resume and resume_from <= 1),
            step_desc=f"[Iter {iteration_id} Step 1B] Filter Repetitive 15-mers",
            command=[
                MERYL, "print", "greater-than", "distinct=0.9998", paths.repet_db,
            ],
            max_log_lines=0,
        )
        if result:
            with open(paths.repet_k15_txt, "w") as fh:
                fh.write(result.stdout)

    # ── Step 2 — Winnowmap → SAM → BAM ──────────────────────────────
    winnowmap_sorted_bam = os.path.join(iter_folder, f"{prefix}.winnowmap.sorted.bam")
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
            "--MD", draft_fasta, reads,
        ],
        stream=True,
    )

    unsorted_bam = os.path.join(iter_folder, f"{prefix}.winnowmap.unsorted.bam")
    conditional_run(
        outfile=unsorted_bam,
        resume_flag=(resume and resume_from <= 2),
        step_desc=f"[Iter {iteration_id} Step 2B] Samtools view",
        command=[
            SAMTOOLS, "view", "--threads", str(num_threads),
            "-hb", "-T", draft_fasta, sam_file,
        ],
        stream=True,
    )

    conditional_run(
        outfile=winnowmap_sorted_bam,
        resume_flag=(resume and resume_from <= 2),
        step_desc=f"[Iter {iteration_id} Step 2C] Samtools sort",
        command=[
            SAMTOOLS, "sort", "--threads", str(num_threads),
            "-o", winnowmap_sorted_bam, unsorted_bam,
        ],
        stream=True,
    )

    # ── Step 3 — FalconC filter ──────────────────────────────────────
    conditional_run(
        outfile=paths.sam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3] FalconC Filter of Alignment",
        command=[
            FALCONC, "bam-filter-clipped", "-t", "-F", "0x104",
            "--input-fn", winnowmap_sorted_bam,
            "--output-fn", paths.sam,
            "--output-count-fn", f"{paths.sam}.filtered_aln_count.txt",
        ],
        stream=True,
    )

    falconc_unsorted_bam = os.path.join(iter_folder, f"{prefix}.falconc.unsorted.bam")
    conditional_run(
        outfile=falconc_unsorted_bam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Convert] Samtools view FalconC SAM → BAM (unsorted)",
        command=[
            SAMTOOLS, "view", "-@", str(num_threads), "-bh",
            "-o", falconc_unsorted_bam, paths.sam,
        ],
        stream=True,
    )

    conditional_run(
        outfile=paths.bam,
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Sort] Samtools sort FalconC BAM",
        command=[
            SAMTOOLS, "sort", "-@", str(num_threads),
            "-o", paths.bam, falconc_unsorted_bam,
        ],
        stream=True,
    )

    conditional_run(
        outfile=paths.bam + ".csi",
        resume_flag=(resume and resume_from <= 3),
        step_desc=f"[Iter {iteration_id} Step 3 Index] Samtools index BAM",
        command=[SAMTOOLS, "index", "-c", paths.bam],
    )

    # ── Step 4 — DeepVariant variant calling ─────────────────────────
    if resume and resume_from <= 4 and os.path.exists(paths.bam):
        try:
            with pysam.AlignmentFile(paths.bam, "rb") as bamfile:
                bamfile.check_index()
        except Exception as e:
            raise InputValidationError(
                f"BAM file {paths.bam} appears to be corrupted or unindexed: {e}"
            ) from e

    if resume and resume_from <= 4 and os.path.exists(paths.dv_vcf):
        try:
            with pysam.VariantFile(paths.dv_vcf, "r"):
                pass
        except Exception as e:
            raise InputValidationError(
                f"VCF file {paths.dv_vcf} is invalid or unreadable: {e}"
            ) from e

    dv_tmp = dv_tmp_dir or os.path.join(iter_folder, "dv_tmp")
    dv_tmp = os.path.abspath(dv_tmp)
    os.makedirs(dv_tmp, exist_ok=True)

    os.environ["TMPDIR"] = dv_tmp
    os.environ["TMP"] = dv_tmp
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
        f"--num_shards={num_threads}",
    ]

    conditional_run(
        outfile=paths.dv_vcf,
        resume_flag=(resume and resume_from <= 4),
        step_desc=f"[Iter {iteration_id} Step 4] DeepVariant",
        command=deepvariant_cmd,
    )

    if not os.path.exists(paths.dv_vcf):
        raise PipelineStepError(
            step=f"[Iter {iteration_id}] DeepVariant",
            returncode=1,
            stderr="DeepVariant did not produce the expected output VCF. Check logs and command.",
        )

    # ── Step 5 — Meryl DB + Merfin polish VCF ───────────────────────
    if os.path.exists(paths.meryl_db) and not os.path.isdir(paths.meryl_db):
        logger.warning(
            "%s exists as a file; removing it so meryl can create a directory.",
            paths.meryl_db,
        )
        os.remove(paths.meryl_db)

    conditional_run(
        outfile=paths.meryl_db,
        resume_flag=(resume and resume_from <= 5),
        step_desc=f"[Iter {iteration_id} Step 5] MerylDB from Genome Draft",
        command=[
            MERYL, "-Q", f"k={k_mer_size}", f"threads={num_threads}",
            "memory=50", "count", "output", paths.meryl_db, draft_fasta,
        ],
        max_log_lines=0,
    )

    merfin_cmd = [
        MERFIN, "-polish",
        "-vcf", paths.dv_vcf,
        "-sequence", draft_fasta,
        "-seqmers", paths.meryl_db,
        "-readmers", readmers,
        "-output", paths.merfin_base,
        "-threads", str(num_threads),
    ]
    if optimized and ideal_kcov and fitted_hist_location:
        merfin_cmd += ["-peak", str(ideal_kcov), "-prob", fitted_hist_location]

    conditional_run(
        outfile=paths.merfin_vcf,
        resume_flag=(resume and resume_from <= 5),
        step_desc=f"[Iter {iteration_id} Step 5] Merfin Polish",
        command=merfin_cmd,
    )

    # ── Step 6 — BCFtools compress → index → consensus ───────────────
    conditional_run(
        outfile=paths.merfin_vcfgz,
        resume_flag=(resume and resume_from <= 6),
        step_desc=f"[Iter {iteration_id} Step 6A] BCFtools compress Merfin VCF",
        command=[BCFTOOLS, "view", "-Oz", paths.merfin_vcf, "-o", paths.merfin_vcfgz],
    )

    conditional_run(
        outfile=paths.merfin_vcfgz + ".tbi",
        resume_flag=(resume and resume_from <= 6),
        step_desc=f"[Iter {iteration_id} Step 6B] BCFtools index Merfin VCF",
        command=[BCFTOOLS, "index", paths.merfin_vcfgz],
    )

    conditional_run(
        outfile=paths.consensus,
        resume_flag=(resume and resume_from <= 6),
        step_desc=(
            f"[Iter {iteration_id} Step 6 Final Consensus] "
            "BCFtools consensus Merfin-polished VCF → Final FASTA"
        ),
        command=[
            BCFTOOLS, "consensus", paths.merfin_vcfgz,
            "-f", draft_fasta, "-H", "1", "-o", paths.consensus,
        ],
    )
    logger.info(
        "Completed iteration %d: consensus at %s", iteration_id, paths.consensus
    )

    if resume_from <= 6 and resume and cleanup:
        logger.info("Cleanup: Removing intermediate files to conserve space")
        for file in [paths.sam, paths.merfin_vcf]:
            if os.path.exists(file):
                os.remove(file)
                logger.info("Cleanup: Removed %s", file)

    # ── Step 7 — Parallel QV: Merfin + Merqury ──────────────────────
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as executor:
        logger.info(
            "Running QV Module: QV and completeness estimation via Merfin and Merqury"
        )
        fut_merfin = executor.submit(
            run_merfin_eval, paths.consensus, readmers, optimized,
            ideal_kcov, fitted_hist_location,
        )
        fut_merqury = executor.submit(
            run_merqury_eval, paths.consensus, readmers, resume=resume,
        )
        merfin_output = fut_merfin.result()
        merqury_output = fut_merqury.result()

    logger.info("Polishing pipeline end: %s", time.ctime())
    return paths.consensus, merfin_output, merqury_output


# ── sub_polish orchestrator ──────────────────────────────────────────


def sub_polish(args: object) -> None:
    """Entry-point for the ``polish`` subcommand — iterates through polishing rounds."""
    cfg = PolishConfig.from_args(args)
    cfg = _apply_resume_kcov(cfg)
    cfg = _apply_config_file(cfg)

    chosen_peak = _resolve_chosen_peak(cfg)
    cfg = _convert_fasta_reads(cfg)
    _validate_inputs(cfg)
    cfg = _ensure_readmers(cfg)
    cfg = _resolve_seq_type(cfg)
    chosen_peak, cfg = _compute_kcov_if_needed(cfg, chosen_peak)
    _check_optimized_resume(cfg)

    logger.info(
        "Launching polishing with ideal_kcov=%s, fitted_hist=%s",
        cfg.ideal_dpeak if cfg.ploidy == "haploid" else cfg.ideal_hpeak,
        cfg.fitted_hist,
    )

    _run_iterations(cfg, chosen_peak)


# ── sub_polish helpers ───────────────────────────────────────────────


def _apply_resume_kcov(cfg: PolishConfig) -> PolishConfig:
    """Reload kcov metadata from a previous run when resuming."""
    kcov_meta = f"{cfg.prefix}.kcov.json"
    if cfg.resume and os.path.exists(kcov_meta):
        with open(kcov_meta) as mf:
            meta = json.load(mf)
        cfg.fitted_hist = meta.get("fitted_hist")
        if cfg.ploidy == "haploid":
            cfg.ideal_dpeak = meta.get("kcov")
        else:
            cfg.ideal_hpeak = meta.get("kcov")
        logger.info(
            "Resume: reloaded kcov=%s, fitted_hist=%s from %s",
            meta.get("kcov"), cfg.fitted_hist, kcov_meta,
        )
    if cfg.resume and cfg.resume_from > 0:
        logger.warning(
            "Both --resume and --resume-from set. Using --resume-from as priority."
        )
    return cfg


def _apply_config_file(cfg: PolishConfig) -> PolishConfig:
    """Merge overrides from an optional JSON config file."""
    if not cfg.config_file:
        return cfg
    with open(cfg.config_file) as cf:
        data = json.load(cf)
    if not cfg.reads and "reads" in data:
        cfg.reads = data["reads"]
    if not cfg.readmers and "readmers" in data:
        cfg.readmers = data["readmers"]
    if cfg.optimized:
        if not cfg.ideal_dpeak and "peak_haploid" in data:
            cfg.ideal_dpeak = data["peak_haploid"]
        if not cfg.ideal_hpeak and "peak_diploid" in data:
            cfg.ideal_hpeak = data["peak_diploid"]
    if not cfg.fitted_hist and "fitted_hist_location" in data:
        cfg.fitted_hist = data["fitted_hist_location"]
    return cfg


def _resolve_chosen_peak(cfg: PolishConfig) -> float | str | None:
    """Determine the chosen coverage peak from existing config."""
    if (
        cfg.resume
        and cfg.fitted_hist
        and (
            (cfg.ploidy == "haploid" and cfg.ideal_dpeak)
            or (cfg.ploidy == "diploid" and cfg.ideal_hpeak)
        )
    ):
        peak = cfg.ideal_dpeak if cfg.ploidy == "haploid" else cfg.ideal_hpeak
        logger.info(
            "Resume: using existing coverage parameters with --ideal_peak=%s and --fitted_hist=%s",
            peak, cfg.fitted_hist,
        )
        return peak
    return None


def _convert_fasta_reads(cfg: PolishConfig) -> PolishConfig:
    """Convert FASTA reads to FASTQ with an appropriate quality score."""
    if not cfg.reads.lower().endswith((".fasta", ".fa")):
        return cfg
    corrector = cfg.read_corrector
    if not corrector:
        logger.warning(
            "FASTA reads detected, but no corrector specified. "
            "Defaulting to 'unknown'. DO NOT use uncorrected FASTA reads!"
        )
        corrector = "unknown"
        cfg.read_corrector = corrector
    quality_score = DEFAULT_CORRECTOR_Q.get(corrector, 30)
    converted_fastq = cfg.prefix + ".converted_reads.fastq"
    logger.info(
        "Corrected Reads Detected: converting FASTA to FASTQ with quality Q=%d for corrector: %s",
        quality_score, corrector,
    )
    fasta_to_fastq(cfg.reads, converted_fastq, quality_score=quality_score)
    cfg.reads = converted_fastq
    return cfg


def _validate_inputs(cfg: PolishConfig) -> None:
    """Check that draft and reads paths are valid."""
    if not cfg.draft or not os.path.exists(cfg.draft):
        raise InputValidationError("Missing or invalid --draft")
    if not cfg.reads or not os.path.exists(cfg.reads):
        raise InputValidationError("Missing or invalid --reads")


def _ensure_readmers(cfg: PolishConfig) -> PolishConfig:
    """Compute readmers database if not provided."""
    if cfg.readmers and os.path.exists(cfg.readmers):
        return cfg
    computed_readmers = cfg.prefix + ".readmers.meryl"
    logger.info("Running Module: Meryl Started: Readmers not provided")
    compute_readmers_db(
        cfg.reads, computed_readmers,
        k=31, threads=cfg.threads,
    )
    cfg.readmers = computed_readmers
    return cfg


def _resolve_seq_type(cfg: PolishConfig) -> PolishConfig:
    """Auto-detect winnowmap seq type from deepseq_type if not set."""
    if cfg.seq_type:
        return cfg
    if cfg.deepseq_type == "PACBIO":
        cfg.seq_type = "pb"
        logger.info(
            "Auto-Setting: setting winnowmap seq_type to 'pb' because deepseq_type = 'PACBIO'"
        )
    elif cfg.deepseq_type == "ONT_R104":
        cfg.seq_type = "ont"
        logger.info(
            "Auto-Setting: setting winnowmap seq_type to 'ont' because deepseq_type = 'ONT_R104'"
        )
    else:
        raise InputValidationError(
            f"deepseq_type={cfg.deepseq_type} isn't automatically mapped. "
            "Please specify --seq_type for winnowmap (e.g. pb or ont)."
        )
    return cfg


def _compute_kcov_if_needed(
    cfg: PolishConfig, chosen_peak: float | str | None
) -> tuple[float | str | None, PolishConfig]:
    """Run the compute-kcov sub-pipeline when --optimized requires it."""
    if needs_kcov(cfg):
        logger.info(
            "Optimized Run: running internal compute-kcov step to determine coverage parameters"
        )
        cov_prefix = cfg.prefix + ".kcov"
        ckargs = PolishConfig(
            draft=cfg.draft,
            reads=cfg.reads,
            prefix=cov_prefix,
            kmer_size=cfg.kmer_size,
            threads=cfg.threads,
            ploidy=cfg.ploidy,
            resume=cfg.resume,
            jellyfish_hash_size=cfg.jellyfish_hash_size,
        )
        discovered_kcov, fitted_hist_location = sub_computekcov(ckargs)
        cfg.fitted_hist = fitted_hist_location
        if cfg.ploidy == "haploid":
            cfg.ideal_dpeak = discovered_kcov
        else:
            cfg.ideal_hpeak = discovered_kcov
        chosen_peak = discovered_kcov
    elif cfg.optimized:
        if cfg.ploidy == "haploid":
            chosen_peak = cfg.ideal_dpeak
            if chosen_peak is None:
                raise InputValidationError("ploidy=haploid but no --ideal_dpeak set.")
        else:
            chosen_peak = cfg.ideal_hpeak
            if chosen_peak is None:
                raise InputValidationError("ploidy=diploid but no --ideal_hpeak set.")
    return chosen_peak, cfg


def _check_optimized_resume(cfg: PolishConfig) -> None:
    """Validate kcov metadata when running --optimized --resume."""
    if not (cfg.optimized and cfg.resume):
        return
    kcov_meta = f"{cfg.prefix}.kcov.json"
    if os.path.exists(kcov_meta):
        logger.warning(
            "Running with --optimized and --resume: loading kcov parameters from %s",
            kcov_meta,
        )
    else:
        raise InputValidationError(
            f"Missing kcov metadata file {kcov_meta} required for --optimized --resume"
        )


def _run_iterations(cfg: PolishConfig, chosen_peak: float | str | None) -> None:
    """Execute the polishing iteration loop and write summaries."""
    os.makedirs(os.path.dirname(cfg.prefix) or ".", exist_ok=True)

    # Iteration 0 — copy user draft
    iteration0 = f"{cfg.prefix}.iter_0.consensus.fasta"
    conditional_run(
        outfile=iteration0,
        resume_flag=(cfg.resume and cfg.resume_from <= 0),
        step_desc="[Iteration 0] Copy initial draft",
        command=["cp", cfg.draft, iteration0],
    )

    all_qv_records: list[tuple[int, str, str]] = []
    old_consensus = os.path.abspath(iteration0)

    for i in range(cfg.iterations):
        iter_num = i + 1
        iter_folder = f"{cfg.prefix}_iteration_{iter_num}"
        os.makedirs(iter_folder, exist_ok=True)

        new_consensus, merfin_out, merqury_out = run_polish_iteration(
            iter_folder=iter_folder,
            iteration_id=iter_num,
            draft_fasta=old_consensus,
            reads=os.path.abspath(cfg.reads),
            readmers=os.path.abspath(cfg.readmers),
            k_mer_size=cfg.kmer_size,
            num_threads=cfg.threads,
            winnowmap_seq_type=cfg.seq_type,
            deepseq_type=cfg.deepseq_type,
            singularity_sif=cfg.singularity_sif,
            cleanup=cfg.cleanup,
            optimized=cfg.optimized,
            ideal_kcov=chosen_peak,
            fitted_hist_location=cfg.fitted_hist,
            resume=cfg.resume,
            resume_from=cfg.resume_from,
            meryl_memory=cfg.meryl_memory,
            dv_tmp_dir=cfg.dv_tmp_dir,
        )

        final_cons_name = f"{cfg.prefix}.iter_{iter_num}.consensus.fasta"
        if not os.path.exists(final_cons_name):
            shutil.copyfile(new_consensus, final_cons_name)

        all_qv_records.append((iter_num, merfin_out, merqury_out))
        iter_qc = os.path.join(iter_folder, "QC_summary.txt")
        with open(iter_qc, "w") as qcfh:
            qcfh.write("----- Merfin -----\n")
            qcfh.write(merfin_out + "\n")
            qcfh.write("----- Merqury -----\n")
            qcfh.write(merqury_out + "\n")
        logger.info("Saved QC summary to %s", iter_qc)
        old_consensus = os.path.abspath(final_cons_name)

    _write_summaries(cfg, all_qv_records)


def _write_summaries(
    cfg: PolishConfig, all_qv_records: list[tuple[int, str, str]]
) -> None:
    """Write text and optional JSON summaries of all iterations."""
    logger.info("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======")
    for iter_num, merfin_text, merqury_text in all_qv_records:
        logger.info("Iteration %d", iter_num)
        logger.info("----- Merfin -----")
        logger.info("%s", merfin_text)
        logger.info("----- Merqury -----")
        logger.info("%s", merqury_text)

    summary_file = f"{cfg.prefix}.QV_Completeness_summary.txt"
    with open(summary_file, "w") as fh:
        fh.write("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======\n\n")
        for iter_num, merfin_text, merqury_text in all_qv_records:
            fh.write(f"> Iteration {iter_num}\n")
            fh.write("----- Merfin -----\n")
            fh.write(merfin_text + "\n")
            fh.write("----- Merqury -----\n")
            fh.write(merqury_text + "\n")
    logger.info("Done: wrote final summary to %s", summary_file)

    if cfg.json_summary:
        from t2t_polish import __version__

        json_records = []
        for iter_num, merfin_text, merqury_text in all_qv_records:
            json_records.append({
                "iteration": iter_num,
                "merfin": merfin_text,
                "merqury": merqury_text,
            })
        json_summary = {
            "pipeline": "T2T-Polish",
            "version": __version__,
            "prefix": cfg.prefix,
            "iterations": cfg.iterations,
            "optimized": cfg.optimized,
            "ploidy": cfg.ploidy,
            "records": json_records,
        }
        json_path = f"{cfg.prefix}.summary.json"
        with open(json_path, "w") as jf:
            json.dump(json_summary, jf, indent=2)
        logger.info("JSON summary written to %s", json_path)


# T2T-Polish v4.1
# Any usage is subject to this software's license.
