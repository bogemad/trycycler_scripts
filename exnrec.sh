#!/bin/bash

#run as bash exnrec.sh ${path_to_cluster_directory} ${reads} ${threads} ${contig_letters_to_exclude} [${max_add_seq_percent}] [${max_length_diff}]

#${contig_letters_to_exclude} is single letter sequence case insensitive and no spaces, if you want to exclude contigs A, B and C from reconcile run then run as run as bash exnrec.sh ${path_to_cluster_directory} ${reads} abc

a=$1
reads=$2
threads=$3
e=${4^^}
maxadd=${5:-2}
max_len_diff=${6:-1.1}

contig_d=$a/1_contigs

while read -n1 d; do
  mkdir -p $contig_d/exclude
  mv $contig_d/${d}_* $contig_d/exclude
done < <(echo -n $e)

trycycler reconcile --reads $reads --threads $threads --cluster_dir $a --max_length_diff $max_len_diff --max_add_seq_percent $maxadd --max_add_seq 200000 --max_trim_seq 200000 --max_trim_seq_percent $maxadd

