# T2T-Automated-Polishing

**Forked from main branch for modification & independent use** Evaluation and polishing workflows for T2T genome assemblies. 
Please cite if any of the codes shared in this repo was used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) doi: https://doi.org/10.1038/s41592-022-01440-3



For exact command lines and workflows used to generate the T2T-CHM13v1.0 and T2T-CHM13v1.1 assemblies, please refer to the [Methods](https://github.com/marbl/CHM13-issues#methods) section in the [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo. Note that some of the tools have been updated since then, and are tracked on this repo.

This README contains details about applying the automated polishing on general genome assemblies. Step by step detail is available in the automated_polshing folder.

The [original script](https://github.com/arangrhie/T2T-Polish/blob/master/automated_polishing/automated-polishing-legacy.sh) used 
in [McCartney et al, 2021](https://doi.org/10.1101/2021.07.02.450803). In this version, the original script from McCartney et al. 2021 and subsequently Arang Rhie's manual edits has been further updated with a manually changed version of the Racon liftover chain.

## Dependencies 
* [Winnowmap2](https://github.com/marbl/Winnowmap)
* [Falconc available in pbipa package](https://github.com/bio-nim/pb-falconc/releases)
* [Racon (liftover branch)](https://github.com/pgrady1322/racon)
* [Meryl v1.3](https://github.com/marbl/meryl)
* [Merfin v1.0](https://github.com/arangrhie/merfin)
* [Samtools](https://github.com/samtools/samtools)
* [BCFtools](https://github.com/samtools/bcftools)


## How to run (Quick Start)

Allocate a fairly large amount of RAM relative to the size of your read set. The Racon step requires the loading of all reads into memory. For instance, a Revio flow cell (~100Gb) requires approximately 400Gb of RAM on a mammalian genome. This pipeline also accepts Herro-corrected ONT reads (use the ONT setting to account for ONT-specific base biases).

```
Automated polishing of draft genomes, version 2.2

Usage:
 automated-polishing_v2.2.sh (options) -d <draft fasta> -r <reads> -m <readmers> -s <sequencer type>

Required Arguments:

-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing
-r	Reads - path to the input reads file, in FASTA/FASTQ format (can be gzipped)
-m	Readmers - path to a Meryl database of read (in order of preference: Illumina - PacBio Hifi - ONT Duplex) k-mers
-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types

Optional Arguments:

-o	Out Prefix - prefix of the output files. Default: AutoPolisher
-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31
-t	Num Threads - number of threads to use. Default: 32
-i	Iterations - number of polishing iterations to perform. Default: 3
```

## Description

This script automatically launches Winnowmap2 to align a read set of choice, then uses falconc, a coordinate based version of Racon, and Merfin to call and filter polishing edits, and finalizes with bcftools to generate the polished consensus in one iteration. The number of iterations can be specified.

## Future Roadmap

1) Winnowmap2 will remap based on a new Meryl k-mer db in successive iterations.

## Installation

A dedicated conda environment is highly recommended. A YML file is available in this repo. Otherwise, install each package independently (or load them on a SLURM-like cluster environment, etc). Note that using package managers for the installation of Racon will lead to a pipeline error with an error code of 'invalid option -L'. Racon must be installed from the following Git repo: https://github.com/pgrady1322/racon, or any available liftover branch from the main Racon repository.

### Conda installation

Use these directions to build into a conda environment for the T2T Automated pipeline with the Conda YML in this Git repo.

```bash
conda env create -f t2tauto.yml -n t2t-auto-polisher
```

Confirm that the modified version of Racon is installed by calling `racon -h` and checking for the `-L` flag in the options. If it is not, install manually via the steps below.

### Self Installation / Installation of Specialized Racon

If you are performing a self install of all dependencies or already have them loaded on a HPC, the modified version of Racon must still be installed. All other dependencies can be from default repositories.

```bash
git clone --recursive https://github.com/pgrady1322/racon.git racon
cd racon
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
