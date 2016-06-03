#!/bin/bash

#align to abc reference using all mapped reads that include all noise paired reads
../../bam2fastq-1.1.0/bam2fastq $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam -o $1#_exon23_high_resolution_multi_ref_mappedreads_sorted.fq
bwa aln -n 0 -k 0 $2_purification_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sai
bwa aln -n 0 -k 0 $2_purification_ref.fa $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sai
bwa sampe $2_purification_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sai $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sai $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted.fq $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sam
samtools view -bS $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.sam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.bam
samtools view -b -F 4 $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads.bam
samtools sort $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads.bam $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads_sorted
samtools index $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads_sorted.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads_sorted.sam

#generate ambiguous mapping candidates to abc genes by checking sam flag info
samtools view -b -F 12 $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_candidates.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_candidates.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_candidates.sam

#deal with abc candidates
linecount=`wc -l $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_candidates.sam | sed 's/ /\t/' | cut -f 1`
echo -e "start filtering out $linecount(whole gene count) ambiguous reads to abc..."

if [[ "$linecount" != "" && "$linecount" -gt "0" ]]; then

i=0
ambicount=0

while [ "$i" -lt "$linecount" ]; do
	i=`expr $i + 1`
	ambiguousread=`sed -n ${i}p $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri_candidates.sam | cut -f 1`
	matchreadinfo=`grep "$ambiguousread" $1_exon23_high_resolution_multi_ref_mappedreads_target.sam`
	matchread=`echo $matchreadinfo | cut -d ' ' -f 1`
	if [ "$matchread" = "$ambiguousread" ]; then
		echo "$matchread" >> $1_ambiguousread_abc.txt
		ambicount=`expr $ambicount + 1`
		echo "$ambicount : filtering out $ambiguousread ..."
		sed -i 's/.*'"$ambiguousread"'\t.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam
	fi

done

sed -i '/^$/d' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam

samtools view -bt exon23_high_resolution_multi_ref.fa.fai $1_exon23_high_resolution_multi_ref_mappedreads_target.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

fi

#clean up temporary files
rm $1_*_exon23_high_resolution_multi_ref_mappedreads_sorted.fq 
rm $1_exon23_high_resolution_multi_ref_mappedreads_sorted_puri.*am
mv *exon23_high_resolution_multi_ref_mappedreads_sorted_puri* ./$2/$1_exon23_high_resolution_multi_ref

#generate purified files
#../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq
