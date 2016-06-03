#!/bin/bash

#deal with abc candidates
linecount=`wc -l $1_ambiguousread_abc.txt | sed 's/ /\t/' | cut -f 1`

if [[ "$linecount" != "" && "$linecount" -gt "0" ]]; then

i=0
ambicount=0

while [ "$i" -lt "$linecount" ]; do
	i=`expr $i + 1`
	ambiguousread=`sed -n ${i}p $1_ambiguousread_abc.txt`
	echo "$i : filtering out $ambiguousread ..."
	sed -i 's/.*'"$ambiguousread"'\t.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.sam

done

sed -i '/^$/d' $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.sam

samtools view -bt exon23_high_resolution_multi_ref.fa.fai $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam

fi

#clean up temporary files
mv $1_ambiguousread_abc.txt ./$2/$1_exon23_high_resolution_multi_ref_flanking

#generate purified files
../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq
