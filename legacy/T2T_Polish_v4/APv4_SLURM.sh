#!/bin/bash
#SBATCH --job-name=Example
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 32
#SBATCH --partition=gpu
#SBATCH --qos=general
#SBATCH --mem=200g
#SBATCH -o %x_%j.stdout
#SBATCH -e %x_%j.stderr

python APv4.py polish -d draft.fasta -r reads.fastq.gz --optimized --singularity_sif deepvariant_1.5.0-gpu.sif