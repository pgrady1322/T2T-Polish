#!/usr/bin/env python3
"""
T2T-Polish v4.0 — polish.py

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
import sys
import time
from contextlib import redirect_stderr, redirect_stdout

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

    with open(os.devnull, "w") as devnull, redirect_stdout(devnull), redirect_stderr(devnull):
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
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
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
            logger.error(
                "BAM file %s appears to be corrupted or unindexed: %s", paths.bam, e
            )
            sys.exit(1)

    if resume and resume_from <= 4 and os.path.exists(paths.dv_vcf):
        try:
            with pysam.VariantFile(paths.dv_vcf, "r"):
                pass
        except Exception as e:
            logger.error(
                "VCF file %s is invalid or unreadable: %s", paths.dv_vcf, e
            )
            sys.exit(1)

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
        logger.error(
            "DeepVariant did not produce the expected output VCF. Check logs and command."
        )
        sys.exit(1)

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


def sub_polish(args: object) -> None:  # noqa: C901
    """Entry-point for the ``polish`` subcommand — iterates through polishing rounds."""
    # Auto-reload kcov from JSON on resume
    prefix: str = getattr(args, "prefix", "AutoPolisher")
    kcov_meta = f"{prefix}.kcov.json"
    if getattr(args, "resume", False) and os.path.exists(kcov_meta):
        with open(kcov_meta) as mf:
            meta = json.load(mf)
        args.fitted_hist = meta.get("fitted_hist")  # type: ignore[attr-defined]
        if getattr(args, "ploidy", "haploid") == "haploid":
            args.ideal_dpeak = meta.get("kcov")  # type: ignore[attr-defined]
        else:
            args.ideal_hpeak = meta.get("kcov")  # type: ignore[attr-defined]
        logger.info(
            "Resume: reloaded kcov=%s, fitted_hist=%s from %s",
            meta.get("kcov"), getattr(args, "fitted_hist", None), kcov_meta,
        )

    resume = getattr(args, "resume", False)
    resume_from = getattr(args, "resume_from", 0)
    if resume and resume_from > 0:
        logger.warning(
            "Both --resume and --resume-from set. Using --resume-from as priority."
        )

    # Config file overrides
    config_file = getattr(args, "config_file", None)
    if config_file:
        with open(config_file) as cf:
            config_data = json.load(cf)
        if not getattr(args, "reads", None) and "reads" in config_data:
            args.reads = config_data["reads"]  # type: ignore[attr-defined]
        if not getattr(args, "readmers", None) and "readmers" in config_data:
            args.readmers = config_data["readmers"]  # type: ignore[attr-defined]
        if getattr(args, "optimized", False):
            if not getattr(args, "ideal_dpeak", None) and "peak_haploid" in config_data:
                args.ideal_dpeak = config_data["peak_haploid"]  # type: ignore[attr-defined]
            if not getattr(args, "ideal_hpeak", None) and "peak_diploid" in config_data:
                args.ideal_hpeak = config_data["peak_diploid"]  # type: ignore[attr-defined]
        if not getattr(args, "fitted_hist", None) and "fitted_hist_location" in config_data:
            args.fitted_hist = config_data["fitted_hist_location"]  # type: ignore[attr-defined]

    chosen_peak: float | str | None = None

    if (
        resume
        and getattr(args, "fitted_hist", None)
        and (
            (getattr(args, "ploidy", "") == "haploid" and getattr(args, "ideal_dpeak", None))
            or (getattr(args, "ploidy", "") == "diploid" and getattr(args, "ideal_hpeak", None))
        )
    ):
        chosen_peak = (
            args.ideal_dpeak
            if getattr(args, "ploidy", "") == "haploid"
            else args.ideal_hpeak
        )
        logger.info(
            "Resume: using existing coverage parameters with --ideal_peak=%s and --fitted_hist=%s",
            chosen_peak,
            getattr(args, "fitted_hist", None),
        )

    # FASTA → FASTQ conversion
    reads: str = getattr(args, "reads", "")
    if reads.lower().endswith((".fasta", ".fa")):
        corrector = getattr(args, "read_corrector", None)
        if not corrector:
            logger.warning(
                "FASTA reads detected, but no corrector specified. "
                "Defaulting to 'unknown'. DO NOT use uncorrected FASTA reads!"
            )
            corrector = "unknown"
            args.read_corrector = corrector  # type: ignore[attr-defined]
        quality_score = DEFAULT_CORRECTOR_Q.get(corrector, 30)
        converted_fastq = prefix + ".converted_reads.fastq"
        logger.info(
            "Corrected Reads Detected: converting FASTA to FASTQ with quality Q=%d for corrector: %s",
            quality_score,
            corrector,
        )
        fasta_to_fastq(reads, converted_fastq, quality_score=quality_score)
        args.reads = converted_fastq  # type: ignore[attr-defined]
        reads = converted_fastq

    # Validate inputs
    draft: str = getattr(args, "draft", "")
    if not draft or not os.path.exists(draft):
        logger.error("Missing or invalid --draft")
        sys.exit(1)
    if not reads or not os.path.exists(reads):
        logger.error("Missing or invalid --reads")
        sys.exit(1)

    readmers_path: str = getattr(args, "readmers", "") or ""
    if not readmers_path or not os.path.exists(readmers_path):
        computed_readmers = prefix + ".readmers.meryl"
        logger.info("Running Module: Meryl Started: Readmers not provided")
        with open(os.devnull, "w") as devnull:
            with redirect_stdout(devnull), redirect_stderr(devnull):
                compute_readmers_db(
                    reads, computed_readmers,
                    k=31, threads=getattr(args, "threads", 8),
                )
        args.readmers = computed_readmers  # type: ignore[attr-defined]
        readmers_path = computed_readmers

    # Auto-detect winnowmap seq type
    seq_type: str = getattr(args, "seq_type", "") or ""
    deepseq_type: str = getattr(args, "deepseq_type", "PACBIO")
    if not seq_type:
        if deepseq_type == "PACBIO":
            seq_type = "pb"
            logger.info(
                "Auto-Setting: setting winnowmap seq_type to 'pb' because deepseq_type = 'PACBIO'"
            )
        elif deepseq_type == "ONT_R104":
            seq_type = "ont"
            logger.info(
                "Auto-Setting: setting winnowmap seq_type to 'ont' because deepseq_type = 'ONT_R104'"
            )
        else:
            logger.error(
                "deepseq_type=%s isn't automatically mapped. "
                "Please specify --seq_type for winnowmap (e.g. pb or ont).",
                deepseq_type,
            )
            sys.exit(1)
        args.seq_type = seq_type  # type: ignore[attr-defined]

    # K-mer coverage computation
    optimized: bool = getattr(args, "optimized", False)
    ploidy: str = getattr(args, "ploidy", "haploid")
    kmer_size: int = getattr(args, "kmer_size", DEFAULT_KMER_SIZE)
    threads: int = getattr(args, "threads", 8)

    if needs_kcov(args):
        logger.info(
            "Optimized Run: running internal compute-kcov step to determine coverage parameters"
        )
        cov_prefix = prefix + ".kcov"

        class CKArgs:
            pass

        ckargs = CKArgs()
        ckargs.resume = resume  # type: ignore[attr-defined]
        ckargs.jellyfish_hash_size = getattr(args, "jellyfish_hash_size", 1_000_000_000)  # type: ignore[attr-defined]
        ckargs.reads = reads  # type: ignore[attr-defined]
        ckargs.prefix = cov_prefix  # type: ignore[attr-defined]
        ckargs.kmer_size = kmer_size  # type: ignore[attr-defined]
        ckargs.threads = threads  # type: ignore[attr-defined]
        ckargs.ploidy = ploidy  # type: ignore[attr-defined]
        discovered_kcov, fitted_hist_location = sub_computekcov(ckargs)
        args.fitted_hist = fitted_hist_location  # type: ignore[attr-defined]
        if ploidy == "haploid":
            args.ideal_dpeak = discovered_kcov  # type: ignore[attr-defined]
        else:
            args.ideal_hpeak = discovered_kcov  # type: ignore[attr-defined]
        chosen_peak = discovered_kcov
    elif optimized:
        if ploidy == "haploid":
            chosen_peak = getattr(args, "ideal_dpeak", None)
            if chosen_peak is None:
                logger.error("ploidy=haploid but no --ideal_dpeak set.")
                sys.exit(1)
        else:
            chosen_peak = getattr(args, "ideal_hpeak", None)
            if chosen_peak is None:
                logger.error("ploidy=diploid but no --ideal_hpeak set.")
                sys.exit(1)

    # Warn on optimized + resume without metadata
    if optimized and resume:
        if os.path.exists(kcov_meta):
            logger.warning(
                "Running with --optimized and --resume: loading kcov parameters from %s",
                kcov_meta,
            )
        else:
            logger.error(
                "Missing kcov metadata file %s required for --optimized --resume",
                kcov_meta,
            )
            sys.exit(1)

    logger.info(
        "Launching polishing with ideal_kcov=%s, fitted_hist=%s",
        getattr(args, "ideal_dpeak" if ploidy == "haploid" else "ideal_hpeak", None),
        getattr(args, "fitted_hist", None),
    )

    os.makedirs(os.path.dirname(prefix) or ".", exist_ok=True)

    # Iteration 0 — copy user draft
    iteration0 = f"{prefix}.iter_0.consensus.fasta"
    conditional_run(
        outfile=iteration0,
        resume_flag=(resume and resume_from <= 0),
        step_desc="[Iteration 0] Copy initial draft",
        command=["cp", draft, iteration0],
    )

    all_qv_records: list[tuple[int, str, str]] = []
    old_consensus = os.path.abspath(iteration0)
    iterations: int = getattr(args, "iterations", 3)
    singularity_sif: str = getattr(args, "singularity_sif", "")
    meryl_memory: int = getattr(args, "meryl_memory", 50)
    dv_tmp_dir: str | None = getattr(args, "dv_tmp_dir", None)
    cleanup: bool = getattr(args, "cleanup", False)

    for i in range(iterations):
        iter_num = i + 1
        iter_folder = f"{prefix}_iteration_{iter_num}"
        os.makedirs(iter_folder, exist_ok=True)

        new_consensus, merfin_out, merqury_out = run_polish_iteration(
            iter_folder=iter_folder,
            iteration_id=iter_num,
            draft_fasta=old_consensus,
            reads=os.path.abspath(reads),
            readmers=os.path.abspath(readmers_path),
            k_mer_size=kmer_size,
            num_threads=threads,
            winnowmap_seq_type=seq_type,
            deepseq_type=deepseq_type,
            singularity_sif=singularity_sif,
            cleanup=cleanup,
            optimized=optimized,
            ideal_kcov=chosen_peak,
            fitted_hist_location=getattr(args, "fitted_hist", None),
            resume=resume,
            resume_from=resume_from,
            meryl_memory=meryl_memory,
            dv_tmp_dir=dv_tmp_dir,
        )

        final_cons_name = f"{prefix}.iter_{iter_num}.consensus.fasta"
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

    # Final summary
    logger.info("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======")
    for iter_num, merfin_text, merqury_text in all_qv_records:
        logger.info("Iteration %d", iter_num)
        logger.info("----- Merfin -----")
        logger.info("%s", merfin_text)
        logger.info("----- Merqury -----")
        logger.info("%s", merqury_text)

    summary_file = f"{prefix}.QV_Completeness_summary.txt"
    with open(summary_file, "w") as fh:
        fh.write("====== FINAL QV/COMPLETENESS SUMMARY (MERFIN + MERQURY) ======\n\n")
        for iter_num, merfin_text, merqury_text in all_qv_records:
            fh.write(f"> Iteration {iter_num}\n")
            fh.write("----- Merfin -----\n")
            fh.write(merfin_text + "\n")
            fh.write("----- Merqury -----\n")
            fh.write(merqury_text + "\n")
    logger.info("Done: wrote final summary to %s", summary_file)


# T2T-Polish v4.0
# Any usage is subject to this software's license.
