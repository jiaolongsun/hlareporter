#!/bin/bash

echo "Formatting classII addition bed..."

#calculate fastq read length

LENGTH=`cat $2/$1_exon23_high_resolution_multi_ref/$1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq | awk '{if(NR%4==2) {print length($1); exit 0;}}'`

if [ ! $LENGTH ]; then
	echo "NULL read length! program teminated!"
	exit
fi

echo "Fastq read length: ${LENGTH}"

#calculate each addition bed range

RANGE1=1
RANGE2=`expr 151 - ${LENGTH} + 9`
RANGE3=`expr 151 + ${LENGTH} - 9`
RANGE4=302

#update classII addition bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DRB1_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DRB1_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DQB1_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DQB1_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DRB3_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DRB3_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DRB4_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DRB4_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DRB5_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DRB5_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DPB1_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DPB1_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_DQA1_addition.bed
sed -i 's/\t.*\t302/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_DQA1_addition.bed
