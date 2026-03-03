# T2T-Automated-Polishing

[![CI](https://github.com/arangrhie/T2T-Polish/actions/workflows/ci.yml/badge.svg)](https://github.com/arangrhie/T2T-Polish/actions/workflows/ci.yml)
[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Singularity](https://img.shields.io/badge/singularity-container-blue.svg)](APv4_Singularity.def)
[![Docker](https://img.shields.io/badge/docker-ready-blue.svg)](Dockerfile)
[![License: Public Domain](https://img.shields.io/badge/license-Public%20Domain-lightgrey.svg)](LICENSE)

Fully automatic K-mer based polishing of genome assemblies.

Current version is unpublished. Please cite this paper, Arang Rhie's T2T-Polish Git Repository, and this Git Repository if any of the code shared in this repo is used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) doi: https://doi.org/10.1038/s41592-022-01440-3

For further details on exact application to T2T-CHM13, read the corresponding section below.

## Repository Contents

This repository includes everything needed to run the APv4 polishing pipeline:

| File | Description |
|------|-------------|
| `t2t_polish/` | Modular Python package (CLI, constants, runner, polishing, evaluation, k-cov) |
| `APv4.py` | Backward-compatible shim that delegates to `t2t_polish.cli:main()` |
| `APv4.yml` | Conda environment specification with all dependencies |
| `Dockerfile` | Docker container definition for portable deployment |
| `APv4_Singularity.def` | Singularity container definition for HPC environments |
| `APv4_SLURM.sh` | Example SLURM submission script for cluster computing |
| `tests/` | pytest test suite |
| `legacy/` | Previous versions (v2, v3) for reference |

## What's New in Version 4

**APv4** represents a major reimplementation of the T2T polishing pipeline:

- **Python-based**: Complete rewrite in Python for improved speed, error handling, and maintainability
- **DeepVariant Integration**: Replaced Racon with GPU-accelerated DeepVariant for more accurate variant calling
- **Automatic QV Assessment**: Real-time quality evaluation at each iteration using both Merfin and Merqury
- **Optimized K-mer Coverage**: Automatic computation of optimal k-mer coverage using Jellyfish and GenomeScope2
- **Resume Capability**: Built-in checkpointing allows resuming from interrupted runs
- **Enhanced Logging**: Comprehensive logging with timestamped entries and detailed diagnostics
- **Parallel Evaluation**: Simultaneous Merfin and Merqury QV calculations for faster assessment

## Description and Best Practices

Auto-Polisher launches an iterative process that allows for more precise K-mer based polishing than typical consensus-only methods. Meryl and Winnowmap2 identify unique k-mers throughout the assembly, and map reads. These reads are extensively filtered in FalconC and Merfin to allow for the best base-level polishes. DeepVariant performs GPU-accelerated variant calling to identify corrections with high precision. Once corrections are made, this process repeats to now re-anchor on k-mers that are now present in the assembly from previous correction. Generally, base-level accuracy peaks at three iterations (the program default). 

Genome assembly accuracy is automatically assessed at each iteration using both Merfin and Merqury, providing real-time QV (Quality Value) and completeness metrics. For final assessment, it is highly recommended to use a hybrid k-mer database filtered for k-mers greater than one to obtain the most accurate Merqury QV. The steps for this via [Meryl and Merqury can be found here](https://github.com/arangrhie/T2T-Polish/tree/master/merqury#2-hybrid), as recommended by the developer, Arang Rhie. Using incomplete Meryl DBs to assess post auto-polisher can lead to inaccurate Merqury QV estimates.


## How to Run (Quick Start)

**Version 4** is implemented as a Python script (`APv4.py`) with improved performance and GPU acceleration support.

### Basic Usage

```bash
# Run diagnostics to check dependencies
python APv4.py diagnostics

# Compute optimal k-mer coverage (recommended first step)
python APv4.py computekcov \
    -r <reads.fastq> \
    -o <output_prefix> \
    -k 21 \
    -t 32 \
    --ploidy haploid

# Run polishing with optimized parameters (automatic k-cov calculation)
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <path/to/deepvariant.sif> \
    --deepseq_type PACBIO \
    -o AutoPolisher \
    -t 32 \
    -i 3 \
    --optimized

# Or run polishing with manual k-mer coverage parameters
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    -m <readmers.meryl> \
    --singularity_sif <path/to/deepvariant.sif> \
    --deepseq_type PACBIO \
    --fitted_hist <genomescope_output/lookup_table.txt> \
    --ideal_dpeak 106.7 \
    -o AutoPolisher \
    -t 32 \
    -i 3
```

### Subcommands

**APv4.py** has two main subcommands:

- **`computekcov`** - Calculate optimal k-mer coverage and fitted histogram using Jellyfish + GenomeScope2
- **`polish`** - Run iterative polishing with DeepVariant-based variant calling

For detailed help:
```bash
python APv4.py --help
python APv4.py computekcov --help
python APv4.py polish --help
```

### Quick Reference: Common Use Cases

**Local workstation (with GPU):**
```bash
python APv4.py polish -d draft.fasta -r reads.fq \
    --singularity_sif deepvariant.sif --deepseq_type PACBIO \
    --optimized -t 32 -i 3
```

**HPC cluster (SLURM):**
```bash
sbatch APv4_SLURM.sh
```

**Docker container:**
```bash
docker run -v $(pwd):/data t2t-polisher:v4 polish \
    -d /data/draft.fasta -r /data/reads.fq \
    --singularity_sif /data/deepvariant.sif --optimized
```

**Resume interrupted run:**
```bash
python APv4.py polish -d draft.fasta -r reads.fq \
    --singularity_sif deepvariant.sif --resume
```

### Key Features

- **GPU Acceleration**: DeepVariant utilizes GPU resources via Singularity for faster variant calling
- **Automatic Resume**: Use `--resume` flag to continue from interrupted runs
- **Real-time QV Assessment**: Both Merfin and Merqury evaluations run automatically after each iteration
- **Flexible Input**: Accepts FASTA or FASTQ reads (FASTQ recommended for HiFi; FASTA auto-converts with quality scores)
- **DeepVariant Models**: Supports WGS, WES, PACBIO, ONT_R104, and HYBRID_PACBIO_ILLUMINA model types

### System Requirements

- **GPU**: Required for DeepVariant acceleration (NVIDIA GPU with CUDA support). **Theoretically, this will work with CPU-enabled DeepVariant, but it has not been tested.**
- **RAM**: Allocate sufficient memory for read processing (e.g., ~400GB for a Revio flow cell on mammalian genomes)
- **Disk Space**: Ensure adequate space for intermediate files and output

## Dependencies 

### Core Dependencies
* [Winnowmap2](https://github.com/marbl/Winnowmap) - Read mapping with repeat-aware k-mer seeding
* [FalconC](https://github.com/bio-nim/pb-falconc/releases) - Alignment filtering (available in pbipa package)
* [DeepVariant](https://github.com/google/deepvariant) - GPU-accelerated variant calling (Singularity image required)
* [Meryl v1.3+](https://github.com/marbl/meryl) - K-mer database operations
* [Merfin v1.0+](https://github.com/arangrhie/merfin) - K-mer based variant filtering and QV evaluation
* [Samtools](https://github.com/samtools/samtools) - SAM/BAM manipulation
* [BCFtools](https://github.com/samtools/bcftools) - VCF processing and consensus generation

### Dependencies for K-mer Coverage Calculation
* [Jellyfish](https://github.com/gmarcais/Jellyfish) - K-mer counting
* [GenomeScope2](https://github.com/tbenavi1/genomescope2.0) - K-mer histogram analysis and coverage estimation

### QV Assessment
* [Merqury](https://github.com/marbl/merqury) - K-mer based assembly evaluation (merqury.sh must be on PATH)

### Python Dependencies
* Python 3.10+
* pysam
* tqdm (for progress bars)
* Standard library: argparse, subprocess, concurrent.futures, logging, pathlib


## Installation

### Conda Environment (Most tested, recommended)

A dedicated Conda environment is highly recommended. A YML file (`APv4.yml`) is available in this repo for easy setup:

```bash
conda env create -f APv4.yml -n t2t-auto-polisher-v4
conda activate t2t-auto-polisher-v4
```

The `APv4.yml` file includes all necessary dependencies:
- Core tools: Winnowmap, FalconC, Meryl, Merfin, Samtools, BCFtools
- K-mer analysis: Jellyfish, GenomeScope2
- Python packages: pysam, tqdm
- Other utilities: Merqury, Racon (for legacy support)

### Container-Based Installation

#### Docker

A Dockerfile is provided for containerized deployment:

```bash
# Build the Docker image
docker build -t t2t-polisher:v4 .

# Run with Docker
docker run -it --rm \
    -v $(pwd)/data:/data \
    t2t-polisher:v4 polish \
    -d /data/draft.fasta \
    -r /data/reads.fastq \
    --singularity_sif /data/deepvariant.sif \
    -o /data/AutoPolisher \
    -t 32
```

**Note**: The Docker container includes the APv4 pipeline and all dependencies except DeepVariant, which must be provided as a Singularity image.

#### Singularity

A Singularity definition file (`APv4_Singularity.def`) is provided for HPC environments:

```bash
# Build the Singularity container
singularity build APv4.sif APv4_Singularity.def

# Run with Singularity
singularity exec APv4.sif APv4.py polish \
    -d draft.fasta \
    -r reads.fastq \
    --singularity_sif deepvariant_gpu.sif \
    -o AutoPolisher \
    -t 32
```

The Singularity container includes:
- Complete APv4 environment with all dependencies
- Merfin v1.1 built from source
- Merqury for QV assessment
- Optimized for HPC cluster deployment

### DeepVariant Singularity Image

DeepVariant must be run via a Singularity container with GPU support. Download the GPU-enabled image:

```bash
# Download DeepVariant GPU Singularity image (recommended)
singularity pull docker://google/deepvariant:"${BIN_VERSION}-gpu"

# Or specify a specific version
singularity pull docker://google/deepvariant:1.6.1-gpu

# Or build from Docker
singularity build deepvariant_gpu.sif docker://google/deepvariant:1.6.1-gpu
```

**Note**: Ensure your system has NVIDIA GPU drivers and CUDA toolkit installed for GPU acceleration.

### HPC/SLURM Example

A SLURM submission script example (`APv4_SLURM.sh`) is provided for cluster environments:

```bash
#!/bin/bash
#SBATCH --job-name=APv4_Polish
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --partition=gpu
#SBATCH --qos=general
#SBATCH --mem=200g
#SBATCH -o %x_%j.stdout
#SBATCH -e %x_%j.stderr

python APv4.py polish \
    -d draft.fasta \
    -r reads.fastq \
    --optimized \
    --singularity_sif deepvariant_1.6.1-gpu.sif \
    -t 32 \
    -o AutoPolisher
```

Submit with:
```bash
sbatch APv4_SLURM.sh
```

Key SLURM considerations:
- Request GPU partition for DeepVariant acceleration
- Allocate sufficient memory (200GB+ for large genomes)
- Match thread count (`-c`) with APv4 threads (`-t`)
- Use `--resume` flag for long-running jobs that may timeout

### Manual Installation

If installing dependencies manually or using HPC modules, ensure all tools are available on PATH. You can verify dependencies by running:

```bash
python APv4.py --diagnostics
```

This will check for all required tools and display their versions.

### Configuration

The pipeline automatically detects tools on PATH. If tools are installed in non-standard locations, you can modify the tool names in `t2t_polish/constants.py`:

```python
# Tool names (modify if needed)
WINNOWMAP = "winnowmap"
FALCONC = "falconc"
MERYL = "meryl"
MERFIN = "merfin"
BCFTOOLS = "bcftools"
JELLYFISH = "jellyfish"
GENOMESCOPE = "genomescope2"
SAMTOOLS = "samtools"
MERQURY_SH = "merqury.sh"
```

These can be changed to full paths if necessary.  


## Advanced Features

### Resume Capability

APv4 includes robust checkpoint/resume functionality. If a run is interrupted, you can resume from the last completed step:

```bash
# Resume from last checkpoint
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --resume

# Resume from a specific step (0=all, 1=Meryl, 2=Winnowmap, 3=FalconC, 4=DeepVariant, 5=Merfin, 6=Consensus)
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --resume \
    --resume-from 4
```

### Optimized Mode with Automatic K-mer Coverage

When using `--optimized`, the pipeline automatically:
1. Computes optimal k-mer coverage using Jellyfish and GenomeScope2
2. Generates fitted histogram for Merfin probability calculations
3. Applies optimal parameters to Merfin polishing
4. Saves coverage parameters to JSON for reliable resume

```bash
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --optimized \
    --ploidy diploid \
    -t 64
```

### Quality Assessment

APv4 automatically runs both Merfin and Merqury evaluations in parallel after each iteration, providing:
- **QV (Quality Value)**: Phred-scaled base accuracy
- **Completeness**: Percentage of expected k-mers present
- **Per-iteration tracking**: Monitor improvement across iterations

Results are written to:
- `<prefix>.QV_Completeness_summary.txt` - Final summary of all iterations
- Per-iteration Merfin outputs: `<prefix>.iter_N.consensus.fasta.merfin_hist.txt`
- Per-iteration Merqury outputs: `<prefix>.iter_N.consensus_merqury.qv`

### Logging

Comprehensive logging is built-in:

```bash
# Enable detailed logging to file
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --log-file polishing.log

# Quiet mode (warnings only to console)
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --quiet
```

### Cleanup Options

Control intermediate file retention:

```bash
# Automatically clean up intermediate files after each iteration
python APv4.py polish \
    -d <draft.fasta> \
    -r <reads.fastq> \
    --singularity_sif <deepvariant.sif> \
    --cleanup
```

## Output Files

APv4 generates organized output with clear naming:

```
<prefix>.iter_0.consensus.fasta           # Initial draft (iteration 0)
<prefix>.iter_1/                          # Iteration 1 working directory
    ├── iter_1.repet_k15.meryl/          # Repetitive k-mer database
    ├── iter_1.winnowmap.sorted.bam      # Aligned reads
    ├── iter_1.falconc.sorted.bam        # Filtered alignments
    ├── iter_1.deepvariant.vcf.gz        # DeepVariant variants
    ├── iter_1.meryl_db/                 # Draft k-mer database
    ├── iter_1.merfin.polish.vcf         # Merfin-filtered variants
    └── iter_1.consensus.fasta           # Polished consensus
<prefix>.iter_2/                          # Iteration 2 working directory
    └── ...
<prefix>.QV_Completeness_summary.txt      # Final QV/completeness summary
<prefix>.tool_versions.txt                # Tool versions used
<prefix>.kcov.json                        # K-mer coverage parameters (if --optimized)
```

## Troubleshooting

### GPU Issues

If DeepVariant fails with GPU errors:
1. Verify NVIDIA drivers: `nvidia-smi`
2. Check CUDA compatibility with DeepVariant version
3. Ensure Singularity has `--nv` flag (automatically added by APv4)
4. Try CPU-only DeepVariant image (slower but more compatible)

### Memory Issues

If jobs fail due to memory:
1. Increase SLURM memory allocation (`--mem=400g` for large genomes)
2. Use `--meryl-memory` to limit Meryl memory usage
3. Clean up intermediate files with `--cleanup`

### Tool Not Found Errors

Run diagnostics to check dependencies:
```bash
python APv4.py diagnostics
```

For missing tools:
- **Conda**: Ensure environment is activated
- **Docker**: Tools are pre-installed in container
- **Singularity**: All tools except DeepVariant are included
- **Manual**: Add tool paths to system PATH or modify `APv4.py`

### Resume Issues

If resume fails:
1. Check for corrupted BAM/VCF files (APv4 validates automatically)
2. Verify `<prefix>.kcov.json` exists when using `--optimized`
3. Use `--resume-from N` to skip specific steps
4. Delete problematic iteration folder and restart

### Container-Specific Issues

**Docker:**
- Mount volumes correctly: `-v $(pwd):/data`
- Provide absolute paths in container: `/data/file.fasta`

**Singularity:**
- Bind mount paths: `singularity exec -B /data:/data`
- Check file permissions on HPC shared filesystems

## Deployment Recommendations

Choose the deployment method that best fits your environment:

| Environment | Recommended Method | Notes |
|-------------|-------------------|-------|
| **Local workstation with GPU** | Conda (`APv4.yml`) | Direct installation, easiest to customize |
| **HPC cluster (SLURM/PBS)** | Singularity (`APv4_Singularity.def`) | Best for shared environments, reproducible |
| **Cloud computing** | Docker (`Dockerfile`) | Portable across cloud providers |
| **Testing/Development** | Conda | Fast iteration, easy debugging |
| **Production pipelines** | Singularity or Docker | Reproducible, version-controlled |

### Choosing Your Setup

**Use Conda if:**
- You have admin/sudo access
- You want to customize tool versions
- You're developing or testing

**Use Singularity if:**
- Running on HPC without root access
- Need reproducible results across different systems
- Want to isolate from system libraries

**Use Docker if:**
- Running on cloud infrastructure
- Have root/Docker access
- Need maximum portability

**Use SLURM script if:**
- Submitting to job scheduler
- Need to queue long-running jobs
- Running on shared HPC resources

## Migrating from Version 3

If you're upgrading from the shell-based v3 pipeline:

### Key Differences

| Feature | Version 3 | Version 4 |
|---------|-----------|-----------|
| Implementation | Bash script | Python |
| Variant Caller | Racon | DeepVariant (GPU) |
| QV Assessment | Manual post-processing | Automatic (Merfin + Merqury) |
| Resume | Limited | Full checkpoint support |
| K-mer Coverage | Manual calculation | Automatic with `--optimized` |
| Input Format | GZIP required | GZIP not required |
| Logging | Basic stdout | Comprehensive with levels |

### Command Translation

**v3 fullauto mode:**
```bash
# Old (v3)
automated-polishing_v3.sh fullauto -d draft.fasta -r reads.fq.gz -s pb -t 32

# New (v4)
python APv4.py polish -d draft.fasta -r reads.fq \
    --singularity_sif deepvariant.sif --deepseq_type PACBIO \
    --optimized -t 32
```

**v3 optimizedpolish mode:**
```bash
# Old (v3)
automated-polishing_v3.sh optimizedpolish -d draft.fasta -r reads.fq.gz \
    -s pb --fitted_hist lookup_table.txt --peak 106.7

# New (v4)
python APv4.py polish -d draft.fasta -r reads.fq \
    --singularity_sif deepvariant.sif --deepseq_type PACBIO \
    --fitted_hist lookup_table.txt --ideal_dpeak 106.7
```

### Notable Changes

- **No GZIP requirement**: v4 accepts uncompressed FASTQ/FASTA
- **Automatic readmers**: When using `--optimized`, readmers DB is computed automatically
- **Read corrector info**: FASTA inputs can specify corrector type (hifiasm, herro, flye) for quality score assignment

## Legacy Version 3

The bash-based version 3 is still available in the `legacy/` directory for users who prefer the original implementation or don't have GPU access. See `legacy/README.md` for v3 documentation.

## Future Roadmap

- Multi-GPU support for DeepVariant
- Integration with cloud computing platforms
- Advanced QV visualization and reporting
- Support for additional variant callers

## T2T-CHM13 Original Resources

For exact command lines and workflows used to generate the T2T-CHM13v1.0 and T2T-CHM13v1.1 assemblies, please refer to the [Methods](https://github.com/marbl/CHM13-issues#methods) section in the [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo. Note that some of the tools have been updated since then, and are tracked on this repo.

This README contains details about applying the automated polishing on general genome assemblies using the latest version 4 implementation.

The [original script](https://github.com/arangrhie/T2T-Polish/blob/master/automated_polishing/automated-polishing-legacy.sh) used in [McCartney et al, 2021](https://doi.org/10.1101/2021.07.02.450803) has been substantially enhanced through version 3 (shell-based with Racon) and now version 4 (Python-based with DeepVariant). Each version represents significant improvements in speed, accuracy, and usability.
