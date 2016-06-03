#!/bin/bash

#align to corresponding reference using all mapped reads that include all noise paired reads
../../bam2fastq-1.1.0/bam2fastq $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam -o $1#_exon23_high_resolution_multi_ref_mappedreads_sorted.fq
bwa aln -n 0 -k 0 drdqdp_etc_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sai
bwa aln -n 0 -k 0 drdqdp_etc_ref.fa $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sai
bwa sampe drdqdp_etc_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sai $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sai $1_1_exon23_high_resolution_multi_ref_mappedreads_sorted.fq $1_2_exon23_high_resolution_multi_ref_mappedreads_sorted.fq > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sam
samtools view -bS $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.sam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.bam
samtools view -b -F 4 $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads.bam
samtools sort $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads.bam $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads_sorted
samtools index $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads_sorted.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads_sorted.sam

#generate ambiguous mapping candidates to targeted genes by checking sam flag info
samtools view -b -F 12 $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_candidates.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_candidates.bam > $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_candidates.sam

#deal with corresponding candidates
linecount=`wc -l $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_candidates.sam | sed 's/ /\t/' | cut -f 1`
echo -e "start filtering out $linecount(whole gene count) ambiguous reads to targeted gene..."

if [[ "$linecount" != "" && "$linecount" -gt "0" ]]; then

i=0
ambicount=0

while [ "$i" -lt "$linecount" ]; do
	i=`expr $i + 1`
	ambiguousread=`sed -n ${i}p $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp_candidates.sam | cut -f 1`
	matchreadinfo=`grep "$ambiguousread" $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam`
	matchread=`echo $matchreadinfo | cut -d ' ' -f 1`
	if [ "$matchread" = "$ambiguousread" ]; then
		echo "$matchread" >> $1_ambiguousread_drdqdp.txt
		ambicount=`expr $ambicount + 1`
		echo "$ambicount : filtering out $ambiguousread ..."
		sed -i 's/.*'"$ambiguousread"'\t.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam
	fi

done

sed -i '/^$/d' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam

samtools view -bt exon23_high_resolution_multi_ref.fa.fai $1_exon23_high_resolution_multi_ref_mappedreads_target.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

fi

#clean up temporary files
rm $1_exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp.*am
rm $1_*_exon23_high_resolution_multi_ref_mappedreads_sorted.fq 
mv *exon23_high_resolution_multi_ref_mappedreads_sorted_drdqdp* ./$2/$1_exon23_high_resolution_multi_ref

#generate purified files
#../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq

