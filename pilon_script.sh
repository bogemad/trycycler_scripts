#!/bin/bash

sample_name=$1
reads_dir=/scratch/xantho_minion/illumina_reads/
threads=32
mkdir -p $sample_name.trycycler/pilon

for c in $sample_name.trycycler/cluster_*; do
    medaka_consensus -i "$c"/4_reads.fastq -d "$c"/7_final_consensus.fasta -o "$c"/medaka -m r941_min_high_g360
    mv "$c"/medaka/consensus.fasta "$c"/8_medaka.fasta
    rm -r "$c"/medaka "$c"/*.fai "$c"/*.mmi  # clean up
done

cat $sample_name.trycycler/cluster_*/8_medaka.fasta > $sample_name.trycycler/pilon/consensus.fasta
cd $sample_name.trycycler/pilon

fastp --in1 $reads_dir/${sample_name}_1.fastq.gz --in2 $reads_dir/${sample_name}_2.fastq.gz --out1 1.fastq.gz --out2 2.fastq.gz --unpaired1 u.fastq.gz --unpaired2 u.fastq.gz

bowtie2-build consensus.fasta consensus.fasta
bowtie2 -1 1.fastq.gz -2 2.fastq.gz -x consensus.fasta --fast --threads 16 -I 0 -X 1000 -S insert_size_test.sam
insert_min=$(python /scratch/xantho_minion/trycycler/insert_size_low.py)
insert_max=$(python /scratch/xantho_minion/trycycler/insert_size_hi.py)
rm *.bt2 insert_size_test.sam

before=consensus

for after in round_1 round_2 round_3 round_4 round_5; do
bowtie2-build "$before".fasta "$before".fasta
bowtie2 -1 1.fastq.gz -2 2.fastq.gz -x "$before".fasta --threads "$threads" -I "$insert_min" -X "$insert_max" --local --very-sensitive-local | samtools sort > illumina_alignments.bam
bowtie2 -U u.fastq.gz -x "$before".fasta --threads "$threads" --local --very-sensitive-local | samtools sort > illumina_alignments_u.bam
samtools index illumina_alignments.bam
samtools index illumina_alignments_u.bam
pilon --genome "$before".fasta --frags illumina_alignments.bam --unpaired illumina_alignments_u.bam --output "$after" --changes --threads $threads
before=$after
rm *.bam *.bam.bai *.bt2
sed -i 's/_pilon//' "$after".fasta
done

cp -av round_5.fasta ../$sample_name.trycycler.fasta
