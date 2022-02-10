#!/bin/bash

#run as bash ~/scripts/trycycler_scripts/quast_longstitch_improvement_assessment.sh FASTA LONGSTITCH_OUTPUT_DIRECTORY THREADS OUTPUT_DIRECTORY REFERENCE_GENOME
#run in loop or parallel for multiple assemblies
#reads must be gzipped

fasta=$1
ls_dir=$2
threads=$3
outdir=$4
ref=$5

mkdir -p $outdir

lfasta=${1%.fasta}
lfasta=${lfasta%.fa}
lfasta=$(basename $lfasta)

quast.py -r $ref -o $outdir/$lfasta -t $threads --min-identity 80.0 $fasta $ls_dir/$lfasta.tigmint-ntLink.longstitch-scaffolds.fa $ls_dir/$lfasta.tigmint-ntLink-arks.longstitch-scaffolds.fa
