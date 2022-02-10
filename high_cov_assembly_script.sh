#!/bin/bash

#run as bash high_cov_assembly_script.sh <sample_name> <reads> <genome_size> <threads> <out_dir>

sample_name=$1
target_depth=50
reads=$2
genome_size=$3
threads=$4
out_dir=$(readlink -f $5)
logs=$out_dir/logs

mkdir -p $out_dir/filtlong/
freads=$out_dir/filtlong/$sample_name.minion.fastq.gz

if [ ! -e $freads ]; then
filtlong --min_length 1000 --keep_percent 95 $reads > $freads || exit 1
fi

mean_length=$(seqtk comp $freads | awk '{count++; bases += $2} END{print bases/count}')
read_count=$(echo $target_depth"*"$genome_size"/"$mean_length | bc)
mkdir -p $out_dir/assemblies

mkdir -p $logs


for ass in flye miniasm raven redbean necat shasta; do
	if [[ $ass == "flye" ]]; then #Assemble w/ Flye
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				echo "Running $ass assembly $i"
				flye --nano-raw $sreads --genome-size "$genome_size" --threads "$threads" --plasmids --out-dir $workdir/assembly &> $logs/$ass.$i.log.txt && echo "$ass assembly $i successful" || echo "$ass assembly $i failed" 
				cp $workdir/assembly/assembly.fasta $out_dir/assemblies/$sample_name.$ass."$i".fasta
				rm -r $workdir
			fi
		done
	elif [[ $ass == "miniasm" ]]; then #Assemble w/ miniasm & minipolish
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				echo "Running $ass assembly $i"
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				miniasm_and_minipolish.sh $sreads "$threads" 1> $workdir/assembly.gfa 2> $logs/$ass.$i.log.txt && echo "$ass assembly $i successful" || echo "$ass assembly $i failed" 
				any2fasta $workdir/assembly.gfa > $out_dir/assemblies/$sample_name.$ass."$i".fasta
				rm -r $workdir
			fi
		done
	elif [[ $ass == "raven" ]]; then #Assemble w/ raven
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				echo "Running $ass assembly $i"
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				cd $workdir
				raven --threads "$threads" $(basename $sreads) 1> $out_dir/assemblies/$sample_name.$ass."$i".fasta 2> $logs/$ass.$i.log.txt && echo "$ass assembly $i successful" || echo "$ass assembly $i failed" 
				cd -
				rm -r $workdir
			fi
		done
	elif [[ $ass == "redbean" ]]; then #Assemble w/ Redbean
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				echo "Running $ass assembly $i"
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				wtdbg2 -o $workdir/assembly -g "$genome_size" -t "$threads" -x ont -i $sreads &> $logs/$ass.$i.log.txt && echo "$ass assembly $i successful" || echo "$ass assembly $i failed" 
				mv $workdir/assembly.cns.fa $out_dir/assemblies/$sample_name.$ass."$i".fasta
				rm -r $workdir
			fi
		done
	elif [[ $ass == "necat" ]]; then #Assemble w/ necat
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				echo "Running $ass assembly $i"
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				echo "$sreads" > $workdir/readlist.txt
				echo "PROJECT=assembly
ONT_READ_LIST=$workdir/readlist.txt
GENOME_SIZE=$genome_size
THREADS=$threads
MIN_READ_LENGTH=1000
PREP_OUTPUT_COVERAGE=40
OVLP_FAST_OPTIONS=-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000
OVLP_SENSITIVE_OPTIONS=-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000
CNS_FAST_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
CNS_SENSITIVE_OPTIONS=-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0
TRIM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400
ASM_OVLP_OPTIONS=-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400
NUM_ITER=2
CNS_OUTPUT_COVERAGE=30
CLEANUP=1
USE_GRID=false
GRID_NODE=0
GRID_OPTIONS=
SMALL_MEMORY=0
FSA_OL_FILTER_OPTIONS=
FSA_ASSEMBLE_OPTIONS=
FSA_CTG_BRIDGE_OPTIONS=
POLISH_CONTIGS=true
" > $workdir/assembly.config
				cd $workdir/
				necat correct $workdir/assembly.config &> $logs/$ass.$i.correct.log.txt && echo "$ass correct $i successful" || echo "$ass correct $i failed" 
				necat assemble $workdir/assembly.config &> $logs/$ass.$i.assemble.log.txt && echo "$ass assemble $i successful" || echo "$ass assemble $i failed" 
				necat bridge $workdir/assembly.config &> $logs/$ass.$i.bridge.log.txt && echo "$ass bridge $i successful" || echo "$ass bridge $i failed" 
				cd -
				mv $workdir/assembly/6-bridge_contigs/polished_contigs.fasta $out_dir/assemblies/$sample_name.$ass."$i".fasta
				rm -r $workdir
			fi
		done
	elif [[ $ass == "shasta" ]]; then #Assemble w/ necat
		for i in {1..3}; do
			if [ ! -e $out_dir/assemblies/$sample_name.$ass."$i".fasta ]; then
				workdir=$out_dir/$sample_name.$ass
				mkdir -p $workdir
				sreads=$workdir/read_sample."$i".fastq
				echo "Running $ass assembly $i"
				seqtk sample -s "$i" $freads "$read_count" | paste - - - - | shuf | tr '\t' '\n' > $sreads
				shasta --input $sreads --config Nanopore-Oct2021 --assemblyDirectory $workdir/assembly --threads $threads --Reads.minReadLength 1000 &> $logs/$ass.$i.log.txt && echo "$ass assembly $i successful" || echo "$ass assembly $i failed" 
				mv $workdir/assembly/Assembly.fasta $out_dir/assemblies/$sample_name.$ass."$i".fasta
				rm -r $workdir
			fi
		done
	else
		echo "Ass is wrong. BITCH!"
		exit 1
	fi
done


trycycler cluster --assemblies $out_dir/assemblies/$sample_name.*.fasta --reads $freads --out_dir $sample_name.trycycler --threads $threads

