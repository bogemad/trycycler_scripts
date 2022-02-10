#!/bin/bash

#run as bash run_longstitch_tigmint-ntLink-arks.sh FASTA READS GENOME_SIZE THREADS OUTPUT_DIRECTORY
#run in loop or parallel for multiple assemblies
#reads must be gzipped

fasta=$1
reads=$2
genome_size=$3
threads=$4
outdir=$5

mkdir -p $outdir

lfasta=${1%.fasta}
lfasta=${lfasta%.fa}
lfasta=$(basename $lfasta)
lreads=${2%.fastq.gz}
lreads=${lreads%.fq.gz}
lreads=$(basename $lreads)

if [ ! -e $outdir/$lfasta.tigmint-ntLink.longstitch-scaffolds.fa ]; then

mkdir $outdir/$lfasta
cp -av $fasta $outdir/$lfasta/$lfasta.fa
cp -av $reads $outdir/$lfasta/$lreads.fq.gz

cd $outdir/$lfasta

longstitch tigmint-ntLink-arks draft=$lfasta reads=$lreads G=9500000 t=$threads

out1=$lfasta.*.tigmint-ntLink.longstitch-scaffolds.fa
out2=$lfasta.*.tigmint-ntLink-arks.longstitch-scaffolds.fa

cp -av $(readlink $out1) ../$lfasta.tigmint-ntLink.longstitch-scaffolds.fa
cp -av $(readlink $out2) ../$lfasta.tigmint-ntLink-arks.longstitch-scaffolds.fa

cd -

fi
