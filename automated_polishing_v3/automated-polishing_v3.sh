#!/bin/bash

#################################################
### Original Design: Ivan Sovic and Ann McCartney
### Date: September 2021
###
### Rewritten: Patrick Grady
### Date: Jan 2024
### Version 2: Jun 2024
### Version 3: Sept 2024
##################################################

### Basic Dependencies: meryl, minimap2, merfin, bcftools, raconL, winnowmap, falconc (pbipa)
### Extended Dependencies: jellyfish, genomescope2

whos=$(basename "$0")

# Dependencies.
RACON=racon
WINNOWMAP=winnowmap
FALCONC=falconc
MERYL=meryl
MERFIN=merfin
BCFTOOLS=bcftools
GENOMESCOPE=genomescope2
JELLYFISH=jellyfish

subcommandusage () {
	echo ""
	echo "Automated polishing of draft genomes, version 3"
	echo ""
	echo "Subcommands:"
	echo ""
	echo "$whos fullauto - Runs the complete pipeline, including automated GenomeScope k-coverage analysis"
	echo "$whos polish - Runs a basic automated polishing run"
	echo "$whos optimizedpolish - Runs an optimized polishing run, requires manual k-cov peak and fitted histogram"
	echo "$whos computekcov - Calculate kcov and fitted histogram for Merfin using GenomeScope"
	echo ""
	echo "For help with each subcommand run:"
	echo "$whos <subcommand> -h|--help"
	echo ""
	exit 0;
}

computekcov_help () {
    echo ""
    echo "Compute K-mer Coverage Usage"
    echo " $0 computekcov (options) -r <reads.gz>"
    echo ""
    echo "Required Arguments:"
    echo ""
    echo "-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped)"
    echo ""
    echo "Optional Arguments:"
    echo ""
    echo "-o	Out Prefix - prefix of the output files. Default: KCoverage"
    echo "-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31"
    echo "-t	Num Threads - number of threads to use. Default: 32"
    echo ""
    exit 0;
}

polish_help () {    
    echo "Polisher Usage:"
    echo " $0 polish (options) -d <draft fasta> -r <reads.gz> -m <readmers> -s <sequencer type>"
    echo ""
    echo "Required Arguments:"
    echo ""
    echo "-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing"
    echo "-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped)"
    echo "-m	Readmers - path to a Meryl database of read (in order of preference: Illumina - PacBio Hifi - ONT Duplex) k-mers"
    echo "-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types"
    echo ""
    echo "Optional Arguments:"
    echo ""
    echo "-o	Out Prefix - prefix of the output files. Default: AutoPolisher"
    echo "-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31"
    echo "-t	Num Threads - number of threads to use. Default: 32"
    echo "-i	Iterations - number of polishing iterations to perform. Default: 3"
	echo ""
	exit 0;
}

optimized_polish_help () {    
    echo "Optimized Polisher Usage:"
    echo " $0 optimizedpolish (options) -d <draft fasta> -r <reads.gz> -m <readmers> -s <sequencer type> -p <Merfin peak value> -f <fitted table>"
    echo ""
    echo "Required Arguments:"
    echo ""
    echo "-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing"
    echo "-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped)"
    echo "-m	Readmers - path to a Meryl database of read (in order of preference: Illumina - PacBio Hifi - ONT Duplex) k-mers"
    echo "-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types"
    echo "-p	Peak Value - The calculated haploid or pseudohaploid peak value (recommend computekcov subcommand)"
    echo "-f	Fitted Table - A lookup table of kmer multiplicities with fitted copy numbers and probabilities (recommend computekcov subcommand)"
    echo ""
    echo "Optional Arguments:"
    echo ""
    echo "-o	Out Prefix - prefix of the output files. Default: AutoPolisher"
    echo "-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31"
    echo "-t	Num Threads - number of threads to use. Default: 32"
    echo "-i	Iterations - number of polishing iterations to perform. Default: 3"
	echo ""
	exit 0;
}

fullauto_help () {    
    echo "Fully Automatic Polisher Usage:"
    echo " $0 fullauto (options) -d <draft fasta> -r <reads.gz> -s <sequencer type>"
    echo ""
    echo "Required Arguments:"
    echo ""
    echo "-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing"
    echo "-r	Reads - path to the input reads file, in FASTA/FASTQ format (MUST be gzipped)"
    echo "-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types"
    echo ""
    echo "Optional Arguments:"
    echo ""
    echo "-o	Out Prefix - prefix of the output files. Default: AutoPolisher"
    echo "-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31"
    echo "-t	Num Threads - number of threads to use. Default: 32"
    echo "-i	Iterations - number of polishing iterations to perform. Default: 3"
	echo ""
	exit 0;
}
 
run_one_standard_iteration () {

	set -e

    # Get the absolute paths.
    mkdir -p $(dirname ${out_prefix})
    out_prefix=$(basename ${out_prefix})
	
    # Generate repetitive 15-mers to downweight.
    local out_winnowmap_bam=${out_prefix}.winnowmap.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} k=15 threads=${num_threads} memory=50 count ${in_draft} output merylDB
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} threads=${num_threads} memory=50 print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

    # Map the reads using Winnowmap.
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${WINNOWMAP} -k 15 -W ${out_winnowmap_bam}.repetitive_k15.txt -t ${num_threads} -ax map-${read_types} --MD ${in_draft} ${in_reads} > ${out_winnowmap_bam}.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    samtools view --threads ${num_threads} -hb -T ${in_draft} ${out_winnowmap_bam}.sam > ${out_winnowmap_bam}

    # Sort the BAM file.
    local out_winnowmap_sorted_bam=${out_prefix}.winnowmap.sorted.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    samtools sort --threads ${num_threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}

    # Filtering the BAM file.
    local out_falconc_sam=${out_prefix}.falconc.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    falconc bam-filter-clipped -t -F 0x104 --input-fn ${out_winnowmap_sorted_bam} --output-fn ${out_falconc_sam} --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt 2>&1 | tee ${out_falconc_sam}.falconc.log

    # Polish using Racon.
    local out_racon_fasta=${out_prefix}.racon.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${RACON} -t ${num_threads} ${in_reads} ${out_falconc_sam} ${in_draft} -L ${out_racon_fasta} -S > ${out_racon_fasta}

    # Generate the Meryl database.
    local out_meryl=${out_prefix}.racon.meryl
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} k=${k_mer_size} threads=${num_threads} memory=50 count ${in_draft} output ${out_meryl}

    # Run Merfin.
    local out_merfin=${out_prefix}.racon.merfin
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERFIN} -polish -sequence ${in_draft} -seqmers ${out_meryl} -readmers ${in_readmers} -peak 106.7 -min 1 -vcf ${out_racon_fasta}.vcf -output ${out_merfin} -threads ${num_threads}

    # Call Consensus
    local out_consensus=${out_prefix}.consensus.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${BCFTOOLS} view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${BCFTOOLS} index ${out_merfin}.polish.vcf.gz
	/usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
	${BCFTOOLS} consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
}

run_one_optimized_iteration () {

	set -e

    # Get the absolute paths.
    mkdir -p $(dirname ${out_prefix})
    out_prefix=$(basename ${out_prefix})
	
    # Generate repetitive 15-mers to downweight.
    local out_winnowmap_bam=${out_prefix}.winnowmap.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} k=15 threads=${num_threads} memory=50 count ${in_draft} output merylDB
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} threads=${num_threads} memory=50 print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

    # Map the reads using Winnowmap.
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${WINNOWMAP} -k 15 -W ${out_winnowmap_bam}.repetitive_k15.txt -t ${num_threads} -ax map-${read_types} --MD ${in_draft} ${in_reads} > ${out_winnowmap_bam}.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    samtools view --threads ${num_threads} -hb -T ${in_draft} ${out_winnowmap_bam}.sam > ${out_winnowmap_bam}

    # Sort the BAM file.
    local out_winnowmap_sorted_bam=${out_prefix}.winnowmap.sorted.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    samtools sort --threads ${num_threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}

    # Filtering the BAM file.
    local out_falconc_sam=${out_prefix}.falconc.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    falconc bam-filter-clipped -t -F 0x104 --input-fn ${out_winnowmap_sorted_bam} --output-fn ${out_falconc_sam} --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt 2>&1 | tee ${out_falconc_sam}.falconc.log

    # Polish using Racon.
    local out_racon_fasta=${out_prefix}.racon.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${RACON} -t ${num_threads} ${in_reads} ${out_falconc_sam} ${in_draft} -L ${out_racon_fasta} -S > ${out_racon_fasta}

    # Generate the Meryl database.
    local out_meryl=${out_prefix}.racon.meryl
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERYL} k=${k_mer_size} threads=${num_threads} memory=50 count ${in_draft} output ${out_meryl}

    # Run Merfin.
    local out_merfin=${out_prefix}.racon.merfin
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${MERFIN} -polish -sequence ${in_draft} -seqmers ${out_meryl} -readmers ${in_readmers} -peak ${ideal_kcov} -min 1 -prob ${fitted_hist_location} -vcf ${out_racon_fasta}.vcf -output ${out_merfin} -threads ${num_threads}

    # Call Consensus
    local out_consensus=${out_prefix}.consensus.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${BCFTOOLS} view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
    ${BCFTOOLS} index ${out_merfin}.polish.vcf.gz
	/usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
	${BCFTOOLS} consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
}

sub_computekcov () {

	set -e

	num_threads='32'
	out_prefix='KCoverage'
	k_mer_size='31'

	OPTSTRING=":t:r:o:k:h"

	while getopts ${OPTSTRING} opt; do
	case ${opt} in
		t)
		num_threads=${OPTARG}
			echo "Threads set to ${threads}."
			;;
		r)
		in_reads=${OPTARG}
			echo "Input reads located at ${in_reads}."
			;;
		o)
		out_prefix=${OPTARG}
			echo "Prefix for all output files set to ${out_prefix}."
			;;
		k)
		k_mer_size=${OPTARG}
			echo "K-mer size set to ${OPTARG}."
			;;
		:)
		echo "Option -${OPTARG} requires an argument."
		exit 1
			;;
		h | * | ?) # Display help.
		echo "Invalid option(s) chosen"
		computekcov_help
		exit 0
		;;
	esac
	done

	shift "$(( OPTIND - 1 ))"

	if [ -z "$in_reads" ]; then
			echo ""
		echo 'Missing -r, input read set. Required arguments are -r.'
			echo ""
			computekcov_help
		exit 1
	fi

	local out_reads_jf=${out_prefix}.reads.jf
	local out_histo_jf=${out_prefix}.reads.histo

	echo ""
	echo "Counting Readmers with Jellyfish"
	echo ""

	${JELLYFISH} count -C -m ${k_mer_size} -s 10000000000 -t ${num_threads} -o ${out_reads_jf} <(zcat $in_reads)

	/usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
	${JELLYFISH} histo -t ${num_threads} ${out_reads_jf} > ${out_histo_jf}

	ideal_kcov=$(${GENOMESCOPE} -i ${out_histo_jf} -o genomescope_out -k ${k_mer_size} --fitted_hist | grep -oP 'kcov:\s*\K\d*\.?\d*')

	fitted_hist_location=$(realpath genomescope_out/lookup_table.txt)

	echo "K Coverage Value = "
	echo ${ideal_kcov}
	echo ""
	echo "Fitted Histogram Location"
	echo ${fitted_hist_location}

}

sub_polish () {

	set -e

	num_threads='32'
	iterations='3'
	out_prefix='AutoPolisher'
	k_mer_size='31'

	OPTSTRING=":t:i:d:r:m:o:s:k:h"

	while getopts ${OPTSTRING} opt; do
	case ${opt} in
		t)
		num_threads=${OPTARG}
			echo "Threads set to ${threads}."
			;;
		i)
		iterations=${OPTARG}
			echo "Iterations set to ${iterations}."
			;;
		d)
		in_draft=${OPTARG}
			echo "Draft fasta located at ${in_draft}."
			;;
		r)
		in_reads=${OPTARG}
			echo "Input reads located at ${in_reads}."
			;;
		m)
		in_readmers=${OPTARG}
			echo "Input readmers located at ${in_readmers}."
			;;
		o)
		out_prefix=${OPTARG}
			echo "Prefix for all output files set to ${out_prefix}."
			;;
		s)
		read_types=${OPTARG}
			((s == pb || s == ont)) || usage
			echo "Read type set to ${read_types}."
			;;
		k)
		k_mer_size=${OPTARG}
			echo "K-mer size set to ${OPTARG}."
			;;
		:)
		echo "Option -${OPTARG} requires an argument."
		exit 1
		;;
		h | * | ?) # Display help.
		echo "Invalid option(s) chosen"
		polish_help
		exit 0
		;;
	esac
	done

	shift "$(( OPTIND - 1 ))"

	if [ -z "$in_draft" ]; then
			echo ""
		echo 'Missing -d, input draft fasta. Required arguments are -d, -r, -m, and -s.'
			echo ""
			polish_help
		exit 1
	fi

	if [ -z "$in_reads" ]; then
			echo ""
		echo 'Missing -r, input read set. Required arguments are -d, -r, -m, and -s.'
			echo ""
			polish_help
		exit 1
	fi

	if [ -z "$in_readmers" ]; then
			echo ""
		echo 'Missing -m, input readmers set. Required arguments are -d, -r, -m, and -s.'
			echo ""
			polish_help
		exit 1
	fi

	if [ -z "$read_types" ]; then
			echo ""
		echo 'Missing -s, input read types. Required arguments are -d, -r, -m, and -s.'
			echo ""
			polish_help
		exit 1
	fi

	mkdir -p $(dirname ${out_prefix})
	cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta


	for (( i = 0 ; i < ${iterations} ; i++ ))
	do next_i=$((i + 1))
		run_one_standard_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${num_threads} ${in_reads} ${in_readmers} ${read_types} ${k_mer_size}
	done

}

sub_optimizedpolish () {

	set -e

	num_threads='32'
	iterations='3'
	out_prefix='AutoPolisher'
	k_mer_size='31'

	OPTSTRING=":t:i:d:r:m:o:s:k:p:h"

	while getopts ${OPTSTRING} opt; do
	case ${opt} in
		t)
		num_threads=${OPTARG}
			echo "Threads set to ${threads}."
			;;
		i)
		iterations=${OPTARG}
			echo "Iterations set to ${iterations}."
			;;
		d)
		in_draft=${OPTARG}
			echo "Draft fasta located at ${in_draft}."
			;;
		r)
		in_reads=${OPTARG}
			echo "Input reads located at ${in_reads}."
			;;
		m)
		in_readmers=${OPTARG}
			echo "Input readmers located at ${in_readmers}."
			;;
		o)
		out_prefix=${OPTARG}
			echo "Prefix for all output files set to ${out_prefix}."
			;;
		s)
		read_types=${OPTARG}
			((s == pb || s == ont)) || usage
			echo "Read type set to ${read_types}."
			;;
		k)
		k_mer_size=${OPTARG}
			echo "K-mer size set to ${OPTARG}."
			;;
		p)
		ideal_kcov=${OPTARG}
			echo "Merfin peak set to ${OPTARG}."
			;;
		f)
		fitted_hist_location=${OPTARG}
			echo "Fitted histogram location ${OPTARG}"
			;;
		:)
		echo "Option -${OPTARG} requires an argument."
		exit 1
			;;
		h | * | ?) # Display help.
		echo "Invalid option(s) chosen"
		optimized_polish_help
		exit 0
		;;
	esac
	done

	shift "$(( OPTIND - 1 ))"

	if [ -z "$in_draft" ]; then
			echo ""
		echo 'Missing -d, input draft fasta. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	if [ -z "$in_reads" ]; then
			echo ""
		echo 'Missing -r, input read set. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	if [ -z "$in_readmers" ]; then
			echo ""
		echo 'Missing -m, input readmers set. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	if [ -z "$read_types" ]; then
			echo ""
		echo 'Missing -s, input read types. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	if [ -z "$ideal_kcov" ]; then
			echo ""
		echo 'Missing -p, Merfin k-cov peak. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	if [ -z "$fitted_hist_location" ]; then
			echo ""
		echo 'Missing -f, fitted histogram location. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			optimized_polish_help
		exit 1
	fi

	mkdir -p $(dirname ${out_prefix})
	cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta

	for (( i = 0 ; i < ${iterations} ; i++ ))
	do next_i=$((i + 1))
		run_one_optimized_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${num_threads} ${in_reads} ${in_readmers} ${read_types} ${k_mer_size} ${ideal_kcov} ${fitted_hist_location}
	done
}

sub_fullauto () {

	set -e

	num_threads='32'
	iterations='3'
	out_prefix='AutoPolisher'
	k_mer_size='31'

	OPTSTRING=":t:i:d:r:m:o:s:k:h"

	while getopts ${OPTSTRING} opt; do
	case ${opt} in
		t)
		num_threads=${OPTARG}
			echo "Threads set to ${threads}."
			;;
		i)
		iterations=${OPTARG}
			echo "Iterations set to ${iterations}."
			;;
		d)
		in_draft=${OPTARG}
			echo "Draft fasta located at ${in_draft}."
			;;
		r)
		in_reads=${OPTARG}
			echo "Input reads located at ${in_reads}."
			;;
		m)
		in_readmers=${OPTARG}
			echo "Input readmers located at ${in_readmers}."
			;;
		o)
		out_prefix=${OPTARG}
			echo "Prefix for all output files set to ${out_prefix}."
			;;
		s)
		read_types=${OPTARG}
			((s == pb || s == ont)) || usage
			echo "Read type set to ${read_types}."
			;;
		k)
		k_mer_size=${OPTARG}
			echo "K-mer size set to ${OPTARG}."
			;;
		:)
		echo "Option -${OPTARG} requires an argument."
		exit 1
		;;
		h | * | ?) # Display help.
		echo "Invalid option(s) chosen"
		fullauto_help
		exit 0
		;;
	esac
	done

	shift "$(( OPTIND - 1 ))"

	if [ -z "$in_draft" ]; then
			echo ""
		echo 'Missing -d, input draft fasta. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			fullauto_help
		exit 1
	fi

	if [ -z "$in_reads" ]; then
			echo ""
		echo 'Missing -r, input read set. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			fullauto_help
		exit 1
	fi

	if [ -z "$read_types" ]; then
			echo ""
		echo 'Missing -s, input read types. Required arguments are -d, -r, -m, -p, -f, and -s. -p, -m, and -f will be autoset in full autopolish mode.'
			echo ""
			fullauto_help
		exit 1
	fi

	mkdir -p $(dirname ${out_prefix})
	cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta

	local out_reads_jf=${out_prefix}.reads.jf
	local out_histo_jf=${out_prefix}.reads.histo

	echo ""
	echo "Counting Readmers with Jellyfish"
	echo ""

	${JELLYFISH} count -C -m ${k_mer_size} -s 10000000000 -t ${num_threads} -o ${out_reads_jf} <(zcat $in_reads)

	/usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
	${JELLYFISH} histo -t ${num_threads} ${out_reads_jf} > ${out_histo_jf}

	ideal_kcov=$(${GENOMESCOPE} -i ${out_histo_jf} -o genomescope_out -k ${k_mer_size} --fitted_hist | grep -oP 'kcov:\s*\K\d*\.?\d*')

	fitted_hist_location=$(realpath genomescope_out/lookup_table.txt)

	echo "K Coverage Value = "
	echo ${ideal_kcov}
	echo ""
	echo "Fitted Histogram Location"
	echo ${fitted_hist_location}

	/usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x\n" >&2 \
	${MERYL} k=${k_mer_size} threads=${num_threads} memory=50 count ${in_reads} output ${out_prefix}.in_mers.meryl
	in_readmers=${out_prefix}.in_mers.meryl

	for (( i = 0 ; i < ${iterations} ; i++ ))
	do next_i=$((i + 1))
		run_one_optimized_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${num_threads} ${in_reads} ${in_readmers} ${read_types} ${k_mer_size} ${ideal_kcov} ${fitted_hist_location}
	done
}


subcommand=$1
case $subcommand in
"" | "-h" | "--help")
    subcommandusage
    ;;
    *)
        shift
        sub_${subcommand} $@
        if [ $? = 127 ]; then
        	echo ""
            echo "Error: '$subcommand' is not a known subcommand." >&2
            subcommandusage
            exit 1
        fi
        ;;
esac
