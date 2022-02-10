#!/bin/bash

sample_name=$1
reads_dir=/scratch/xantho_minion/guppy/fastq
threads=16
freads=filtlong/$sample_name.minion.filtlong.fastq
genome_size=5000000
mkdir -p assemblies

filtlong --min_length 1000 --keep_percent 95 $reads_dir/$sample_name.minion.fastq.gz > $freads

#Assemble w/ Flye
flye --nano-raw $freads --genome-size "$genome_size" --threads "$threads" --plasmids --out-dir $sample_name.flye
cp $sample_name.flye/assembly.fasta assemblies/$sample_name.flye.fasta
rm -r $sample_name.flye

#Assemble w/ miniasm & minipolish
miniasm_and_minipolish.sh $freads "$threads" > $sample_name.miniasm.gfa
any2fasta $sample_name.miniasm.gfa > assemblies/$sample_name.miniasm.fasta
rm $sample_name.miniasm.gfa

#Assemble w/ raven
mkdir $sample_name.raven
cd $sample_name.raven
raven --threads "$threads" ../$freads > ../assemblies/$sample_name.raven.fasta
cd ..
rm -r $sample_name.raven

#Assemble w/ Redbean
wtdbg2.pl -o $sample_name.redbean -g "$genome_size" -t "$threads" -x ont $freads
mv $sample_name.redbean.cns.fa assemblies/$sample_name.redbean.fasta
rm $sample_name.redbean.*

trycycler cluster --assemblies assemblies/$sample_name.*.fasta --reads $freads --out_dir $sample_name.trycycler --threads $threads

