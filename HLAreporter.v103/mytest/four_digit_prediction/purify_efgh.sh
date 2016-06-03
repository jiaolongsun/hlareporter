#!/bin/bash

../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq

bwa aln -n 0 -k 0 efgh_etc_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_1_exon23_high_resolution_multi_ref_efgh.sai
bwa aln -n 0 -k 0 efgh_etc_ref.fa $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_2_exon23_high_resolution_multi_ref_efgh.sai

#handle each fq seperately
bwa samse efgh_etc_ref.fa $1_1_exon23_high_resolution_multi_ref_efgh.sai $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_1_exon23_high_resolution_multi_ref_efgh.sam
bwa samse efgh_etc_ref.fa $1_2_exon23_high_resolution_multi_ref_efgh.sai $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_2_exon23_high_resolution_multi_ref_efgh.sam
samtools view -bS $1_1_exon23_high_resolution_multi_ref_efgh.sam > $1_1_exon23_high_resolution_multi_ref_efgh.bam
samtools view -bS $1_2_exon23_high_resolution_multi_ref_efgh.sam > $1_2_exon23_high_resolution_multi_ref_efgh.bam
samtools view -b -F 4 $1_1_exon23_high_resolution_multi_ref_efgh.bam > $1_1_exon23_high_resolution_multi_ref_efgh_mappedreads.bam
samtools view -b -F 4 $1_2_exon23_high_resolution_multi_ref_efgh.bam > $1_2_exon23_high_resolution_multi_ref_efgh_mappedreads.bam

samtools view $1_1_exon23_high_resolution_multi_ref_efgh_mappedreads.bam > $1_1_exon23_high_resolution_multi_ref_efgh_mappedreads.sam
samtools view $1_2_exon23_high_resolution_multi_ref_efgh_mappedreads.bam > $1_2_exon23_high_resolution_multi_ref_efgh_mappedreads.sam

#deal with sam 1
linecount=`wc -l $1_1_exon23_high_resolution_multi_ref_efgh_mappedreads.sam | sed 's/ /\t/' | cut -f 1`
echo -e "start filtering $linecount ambiguous reads from sam 1..."

if [[ "$linecount" != "" && "$linecount" -gt "0" ]]; then

i=0
while [ "$i" -lt "$linecount" ]; do
	i=`expr $i + 1`
	ambiguousread=`sed -n ${i}p $1_1_exon23_high_resolution_multi_ref_efgh_mappedreads.sam | cut -f 1`
	sed -i 's/.*'"$ambiguousread"'\t.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam
	echo "$ambiguousread" >> $1_ambiguousread_efgh.txt
	echo "$i : filtering out $ambiguousread ..."
done

sed -i '/^$/d' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam

samtools view -bt exon23_high_resolution_multi_ref.fa.fai $1_exon23_high_resolution_multi_ref_mappedreads_target.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

fi

#deal with sam 2
linecount=`wc -l $1_2_exon23_high_resolution_multi_ref_efgh_mappedreads.sam | sed 's/ /\t/' | cut -f 1`
echo -e "start filtering $linecount ambiguous reads from sam 2..."

if [[ "$linecount" != "" && "$linecount" -gt "0" ]]; then

i=0
while [ "$i" -lt "$linecount" ]; do
	i=`expr $i + 1`
	ambiguousread=`sed -n ${i}p $1_2_exon23_high_resolution_multi_ref_efgh_mappedreads.sam | cut -f 1`
	sed -i 's/.*'"$ambiguousread"'\t.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam
	echo "$ambiguousread" >> $1_ambiguousread_efgh.txt
	echo "$i : filtering out $ambiguousread ..."
done

sed -i '/^$/d' $1_exon23_high_resolution_multi_ref_mappedreads_target.sam

samtools view -bt exon23_high_resolution_multi_ref.fa.fai $1_exon23_high_resolution_multi_ref_mappedreads_target.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

fi

#clean up temporary files
mv $1_*_exon23_high_resolution_multi_ref_efgh* ./$2/$1_exon23_high_resolution_multi_ref
rm $1_*_exon23_high_resolution_multi_ref_mappedreads_target.fq

#generate purified files
#../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq

