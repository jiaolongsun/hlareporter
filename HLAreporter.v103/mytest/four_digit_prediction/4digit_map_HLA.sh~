#!/bin/bash

if [ $# -ne 3 ]; then
	echo Usage: command fq1 fq2 sample_name
	exit
fi

#bwa index -a is exon23_high_resolution_multi_ref.fa
bwa aln exon23_high_resolution_multi_ref.fa $1 > $3_1_exon23_high_resolution_multi_ref.sai
bwa aln exon23_high_resolution_multi_ref.fa $2 > $3_2_exon23_high_resolution_multi_ref.sai
bwa sampe exon23_high_resolution_multi_ref.fa $3_1_exon23_high_resolution_multi_ref.sai $3_2_exon23_high_resolution_multi_ref.sai $1 $2 > $3_exon23_high_resolution_multi_ref.sam 
samtools view -bS $3_exon23_high_resolution_multi_ref.sam > $3_exon23_high_resolution_multi_ref.bam
samtools view -b -F 4 $3_exon23_high_resolution_multi_ref.bam > $3_exon23_high_resolution_multi_ref_mappedreads.bam

samtools sort $3_exon23_high_resolution_multi_ref_mappedreads.bam $3_exon23_high_resolution_multi_ref_mappedreads_sorted
samtools index $3_exon23_high_resolution_multi_ref_mappedreads_sorted.bam
rm $3_exon23_high_resolution_multi_ref.sam 
rm $3_exon23_high_resolution_multi_ref.bam
rm $3_exon23_high_resolution_multi_ref_mappedreads.bam
