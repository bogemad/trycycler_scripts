#!/bin/bash

sample_name=$1
threads=32
freads=filtlong/$sample_name.minion.filtlong.fastq
genome_size=5000000

for c in $sample_name.trycycler/cluster_*; do
	trycycler msa --cluster_dir $c --threads $threads
done

trycycler partition --reads $freads --threads $threads --cluster_dirs $sample_name.trycycler/cluster_*

for c in $sample_name.trycycler/cluster_*; do
	trycycler consensus --threads $threads --cluster_dir $c
done

