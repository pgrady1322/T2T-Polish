#!/bin/bash

##############################################
### Authors: Ivan Sovic and Ann Mc Cartney ###
### Date: September 2021 		   ###
###					   ###	
### Rewritten: Patrick Grady	           ###
### Date: Jan 2024                         ###
### Version 2: Jun 2024			   ###
##############################################

### Dependencies: meryl, minimap2, merfin, bcftools, raconL, winnowmap, falconc (pbipa)

usage() {
    echo "Automated polishing of draft genomes, version 2.2"
    echo ""
    echo "Usage:"
    echo " $0 (options) -d <draft fasta> -r <reads> -m <readmers> -s <sequencer type>"
    echo ""
    echo "Required Arguments:"
    echo ""
    echo "-d	Draft Fasta - path to the input FASTA/FASTQ containing draft sequences for polishing"
    echo "-r	Reads - path to the input reads file, in FASTA/FASTQ format (can be gzipped)"
    echo "-m	Readmers - path to a Meryl database of read (in order of preference: Illumina - PacBio Hifi - ONT Duplex) k-mers"
    echo "-s	Sequencing Type - pb or ont, use pb for HiFi and Illumina, ont for all ONT read types"
    echo ""
    echo "Optional Arguments:"
    echo ""
    echo "-o	Out Prefix - prefix of the output files. Default: AutoPolisher"
    echo "-k	K-mer Size - Meryl database k-mer size, should follow initial input. Default (recommended): k=31"
    echo "-t	Num Threads - number of threads to use. Default: 32"
    echo "-i	Iterations - number of polishing iterations to perform. Default: 3"
    echo "";
    exit 0;
}

threads='32'
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
      usage
      exit 0
      ;;
  esac
done

shift "$(( OPTIND - 1 ))"

if [ -z "$in_draft" ]; then
        echo ""
	echo 'Missing -d, input draft fasta. Required arguments are -d, -r, -m, and -s.'
        echo ""
	usage
	exit 1
fi

if [ -z "$in_reads" ]; then
        echo ""
	echo 'Missing -r, input read set. Required arguments are -d, -r, -m, and -s.'
        echo ""
        usage
	exit 1
fi

if [ -z "$in_readmers" ]; then
        echo ""
	echo 'Missing -m, input readmers set. Required arguments are -d, -r, -m, and -s.'
        echo ""
        usage
	exit 1
fi

if [ -z "$read_types" ]; then
        echo ""
	echo 'Missing -s, input read types. Required arguments are -d, -r, -m, and -s.'
        echo ""
        usage
	exit 1
fi

# Dependencies.
RACON=racon
WINNOWMAP=winnowmap
FALCONC=falconc
MERYL=meryl
MERFIN=merfin
BCFTOOLS=bcftools

run_one_iteration () {

    # Get the absolute paths.
    mkdir -p $(dirname ${out_prefix}_${next_i})
    out_prefix=$(realpath ${out_prefix})

    # Generate repeitive 15-mers to downweight.
    local out_winnowmap_bam=${out_prefix}.winnowmap.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${MERYL} count k=15 ${in_draft} output merylDB
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${MERYL} print greater-than distinct=0.9998 merylDB > ${out_winnowmap_bam}.repetitive_k15.txt

    # Map the reads using Winnowmap.
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${WINNOWMAP} -k 15 -W ${out_winnowmap_bam}.repetitive_k15.txt -t ${num_threads} -ax map-${read_type} --MD ${in_draft} ${in_dataset} | samtools view -hb -T ${in_draft} > ${out_winnowmap_bam}

    # Sort the BAM file.
    local out_winnowmap_sorted_bam=${out_prefix}.winnowmap.sorted.bam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    samtools sort --threads ${num_threads} -o ${out_winnowmap_sorted_bam} ${out_winnowmap_bam}

    # Filtering the BAM file.
    local out_falconc_sam=${out_prefix}.falconc.sam
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    falconc bam-filter-clipped -t -F 0x104 --input-fn ${out_winnowmap_sorted_bam} --output-fn ${out_falconc_sam} --output-count-fn ${out_falconc_sam}.filtered_aln_count.txt 2>&1 | tee ${out_falconc_sam}.falconc.log

    # Polish using Racon.
    local out_racon_fasta=${out_prefix}.racon.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${RACON} -t ${num_threads} ${in_dataset} ${out_falconc_sam} ${in_draft} -L ${out_racon_fasta} -S > ${out_racon_fasta}

    # Generate the Meryl database.
    local out_meryl=${out_prefix}.racon.meryl
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${MERYL} count k=${k_mer_size} ${in_draft} output ${out_meryl}

    # Run Merfin.
    local out_merfin=${out_prefix}.racon.merfin
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${MERFIN} -polish -sequence ${in_draft} -seqmers ${out_meryl} -readmers ${in_readmers} -peak 106.7 -vcf ${out_racon_fasta}.vcf -output ${out_merfin} -threads ${num_threads}

    # Call Consensus
    local out_consensus=${out_prefix}.consensus.fasta
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${BCFTOOLS} view -Oz ${out_merfin}.polish.vcf > ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${BCFTOOLS} index ${out_merfin}.polish.vcf.gz
    /usr/bin/time --format="cmd: %C\\nreal_time: %e s\\nuser_time: %U s\\nsys_time: %S s\\nmax_rss: %M kB\\nexit_status: %x" >&2 \
    ${BCFTOOLS} consensus ${out_merfin}.polish.vcf.gz -f ${in_draft} -H 1 > ${out_consensus}
}

run_all () {
    cp ${in_draft} ${out_prefix}.iter_0.consensus.fasta
    for (( i = 0 ; i < ${iterations} ; i++ ))
	do next_i=$((i + 1))
        run_one_iteration ${out_prefix}.iter_${next_i} ${out_prefix}.iter_${i}.consensus.fasta ${num_threads} ${in_reads} ${in_readmers} ${read_types} ${k_mer_size}
    done
}

run_all ${out_prefix} ${num_threads} ${iterations} ${in_draft} ${in_reads} ${in_readmers} ${read_types} ${k_mer_size}
