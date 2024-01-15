#!/bin/bash

#run as bash sr_polish.sh <long_read_assembly.fasta> <short_reads.r1.fastq.gz> <short_reads.r2.fastq.gz> <polished_output.fasta> <threads>

infasta=$(readlink -f $1)
reads1=$(readlink -f $2)
reads2=$(readlink -f $3)
out=$(readlink -f $4)
tmp=$infasta.polish_tmp
t=$5

mkdir $tmp
cd $tmp
cp -av $infasta .
raw=$(basename $infasta)

bwa index $raw
bwa mem -t $t -a $raw $reads1 > alignments_1.sam
bwa mem -t $t -a $raw $reads2 > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
rm alignments_1.sam alignments_2.sam
polypolish $raw filtered_1.sam filtered_2.sam > 1_polypolish.fasta
polca.sh -a 1_polypolish.fasta -r "$reads1 $reads2" -t $t -m 1G
mv *.PolcaCorrected.fa 2_polypolish_polca.fasta
rm filtered_1.sam filtered_2.sam
bwa index 2_polypolish_polca.fasta
bwa mem -t $t -a 2_polypolish_polca.fasta $reads1 > alignments_1.sam
bwa mem -t $t -a 2_polypolish_polca.fasta $reads2 > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
rm alignments_1.sam alignments_2.sam
polypolish 2_polypolish_polca.fasta filtered_1.sam filtered_2.sam  > 3_polypolish_polca_polypolish.fasta
polca.sh -a 3_polypolish_polca_polypolish.fasta -r "$reads1 $reads2" -t $t -m 1G
mv *.PolcaCorrected.fa 4_polypolish_polca_polypolish_polca.fasta
rm filtered_1.sam filtered_2.sam

mv 4_polypolish_polca_polypolish_polca.fasta $out

cd ..
rm -r $tmp
