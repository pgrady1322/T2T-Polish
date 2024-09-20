# T2T-Automated-Polishing

Fully automatic K-mer based polishing of genome assemblies.

Current version is unpublished. Please cite this paper, Arang Rhie's T2T-Polish Git Repository, and this Git Repository if any of the code shared in this repo is used:

Mc Cartney AM, Shafin K, Alonge M et al. Chasing perfection: validation and polishing strategies for telomere-to-telomere genome assemblies. Nat Methods (2022) doi: https://doi.org/10.1038/s41592-022-01440-3

For further details on exact application to T2T-CHM13, read the corresponding section below.

## Description and Best Practices

Auto-Polisher launches an iterative process that allows for more precise K-mer based polishing than typical consensus-only methods. Meryl and Winnowmap2 identify unique k-mers throughout the assembly, and map reads. These reads are extensively filtered in falconc and Merfin to allow for the best base-level polishes. Once corrections are made, this process repeats to now re-anchor on k-mers that are now present in the assembly from previous correction. Generally, base-level accuracy peaks at three iterations (the program default). Genome assembly accuracy can be assessed post-polish by Merqury, and it is highly recommended to use a hybrid k-mer database filtered for k-mers greater than one to obtain the most accurate Merqury QV. The steps for this via [Meryl and Merqury can be found here](https://github.com/arangrhie/T2T-Polish/tree/master/merqury#2-hybrid), as recommended by the developer, Arang Rhie. Using incomplete Meryl DBs to assess post auto-polisher can lead to inaccurate Merqury QV estimates.


## How to Run (Quick Start)

Allocate a fairly large amount of RAM relative to the size of your read set. The Racon step requires the loading of all reads into memory. For instance, a Revio flow cell (~100Gb) requires approximately 400Gb of RAM on a mammalian genome. This pipeline accepts (**and highly recommends**) Herro-corrected ONT reads (use the ONT setting to account for ONT-specific base biases) produced with either the original [Herro](https://github.com/lbcb-sci/herro) (R9 and R10 reads) or [Dorado correct](https://github.com/nanoporetech/dorado) (only R10 at time of v3 release).

```
Automated polishing of draft genomes, version 3

Subcommands:

automated-polishing_v3.sh fullauto - Runs the complete pipeline, including automated GenomeScope k-coverage analysis
automated-polishing_v3.sh polish - Runs a basic automated polishing run
automated-polishing_v3.sh optimizedpolish - Runs an optimized polishing run, requires manual k-cov peak and fitted lookup table from GenomeScope2
automated-polishing_v3.sh computekcov - Calculate kcov and fitted histogram for Merfin using GenomeScope

For help with each subcommand run:
automated-polishing_v3.sh <subcommand> -h|--help
```

```
Fully Automatic Polisher Usage:
 automated-polishing_v3.sh fullauto (options) -d <draft fasta> -r <reads.gz> -s <sequencer type>

Required Arguments:

-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing
-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped).
-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types

Optional Arguments:

-o	Out Prefix - prefix of the output files. Default: AutoPolisher
-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31
-t	Num Threads - number of threads to use. Default: 32
-i	Iterations - number of polishing iterations to perform. Default: 3
```

## Dependencies 
* [Winnowmap2](https://github.com/marbl/Winnowmap)
* [Falconc, available in pbipa package](https://github.com/bio-nim/pb-falconc/releases)
* [Racon (liftover branch)](https://github.com/pgrady1322/racon)
* [Meryl v1.3](https://github.com/marbl/meryl)
* [Merfin v1.0](https://github.com/arangrhie/merfin)
* [Samtools](https://github.com/samtools/samtools)
* [BCFtools](https://github.com/samtools/bcftools)

## Dependencies for FullAuto and ComputeKCov modes
* [Jellyfish](https://github.com/gmarcais/Jellyfish)
* [GenomeScope2](https://github.com/tbenavi1/genomescope2.0)


## Installation

A dedicated environment is highly recommended, but not strictly necessary. A YML file is available in this repo. Otherwise, each dependency must be installed independently (or loaded on an HPC with shared modules, etc).

**Important Note**: Using package managers for the installation of Racon will lead to a pipeline error with an error code of 'invalid option -L'. Racon must be installed from the following Git repo: https://github.com/pgrady1322/racon, or any available liftover branch from the main Racon repository.

### Conda Installation of Most Dependencies

Use these directions to build into a conda environment for the T2T Automated pipeline with the Conda YML in this Git repo.

```bash
conda env create -f auto_polisher_env_simple_v3.yml -n t2t-auto-polisher
```

### Installation of Specialized Racon

If you are performing a self install of all dependencies or already have them loaded on a HPC, the modified version of Racon must still be installed. All other dependencies can be from default repositories.

```bash
git clone --recursive https://github.com/pgrady1322/racon.git racon
cd racon
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```

### Installation of Dependencies (Non-Conda YML)

The remaining dependencies can be installed via pip, conda, or whatever your preferred installation method is (module loaded on an HPC, compiled, etc.). As long as they are on PATH, the auto-polisher will be able to use them.

If strictly necessary, find the dependencies block in the automated_polisher_v# script and modify the values as follows:

```bash
# Dependencies.
RACON=racon
WINNOWMAP=winnowmap
FALCONC=falconc
MERYL=meryl
MERFIN=merfin
BCFTOOLS=bcftools
```

Modify the **lowercase** values to whatever the call for that program is on your computing environment. These can be full paths. 


## Future Roadmap

1) GPU support.
2) Automatic assessment of QV score at each step (the pipeline performs this internally, but it is not yet output).


## T2T-CHM13 Original Resources

For exact command lines and workflows used to generate the T2T-CHM13v1.0 and T2T-CHM13v1.1 assemblies, please refer to the [Methods](https://github.com/marbl/CHM13-issues#methods) section in the [CHM13-Issues](https://github.com/marbl/CHM13-issues) repo. Note that some of the tools have been updated since then, and are tracked on this repo.

This README contains details about applying the automated polishing on general genome assemblies. Step by step detail is available in the automated_polshing folder.

The [original script](https://github.com/arangrhie/T2T-Polish/blob/master/automated_polishing/automated-polishing-legacy.sh) used 
in [McCartney et al, 2021](https://doi.org/10.1101/2021.07.02.450803). In this version, the original script from McCartney et al. 2021 and subsequently Arang Rhie's manual edits has been further updated with a manually changed version of the Racon liftover chain.