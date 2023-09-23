#!/bin/bash

if [[ -z $1 ]] ; then
    echo "Usage: ./init.sh <asm.fasta>"
    echo "  <asm.fasta>: absolute path to asm.fasta,"
    echo "               expecting to have asm.fasta.fai in the same path"
    exit 0
fi

asm=$1
name=`echo $asm | sed 's/\.gz$//g' | sed 's/\.fa$//g' | sed 's/\.fasta$//g'`

if [[ "$asm" == "*\.gz" ]]; then
  pigz -dc $asm > $name.fa
else
  ln -sf $asm     $name.fa
  ln -sf $asm.fai $name.fa.fai
fi

module load seqtk

if [[ ! -s $asm.fai ]]; then
  echo "No $asm.fai found. Exit."
  exit -1
fi

echo "# Get sizes"
awk '{print $1"\t0\t"$2}' $asm.fai > $name.bed

echo "# Telomere"
seqtk telo $asm > $name.telo.bed

echo "# Exclude"
seqtk gap -l0 $asm > $name.exclude.bed

echo "# Pattern"
$tools/T2T-Polish/pattern/microsatellites.sh $name.fa
