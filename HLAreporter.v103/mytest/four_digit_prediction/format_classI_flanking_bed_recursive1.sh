#!/bin/bash

echo "Formatting classI flanking bed..."

#calculate fastq read length

LENGTH=`cat $2/$1_exon23_high_resolution_multi_ref/$1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq | awk '{if(NR%4==2) {print length($1); exit 0;}}'`

if [ ! $LENGTH ]; then
	echo "NULL read length! program teminated!"
	exit
fi

echo "Fastq read length: ${LENGTH}"

#calculate each flanking bed range

RANGE1=1
RANGE2=`expr 185 - 0` #5
RANGE3=`expr 185 + 0` #5
RANGE4=`expr 370 + 188 - 0` #15
RANGE5=`expr 370 + 188 + 0`
RANGE6=746

#update classI flanking bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_A_flanking.bed
sed -i 's/\t.*\t746/\t'"$RANGE5"'\t'"$RANGE6"'/g' HLA_A_flanking.bed
sed -i 's/\t.*\t[3456]../\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_A_flanking.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_B_flanking.bed
sed -i 's/\t.*\t746/\t'"$RANGE5"'\t'"$RANGE6"'/g' HLA_B_flanking.bed
sed -i 's/\t.*\t[3456]../\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_B_flanking.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_C_flanking.bed
sed -i 's/\t.*\t746/\t'"$RANGE5"'\t'"$RANGE6"'/g' HLA_C_flanking.bed
sed -i 's/\t.*\t[3456]../\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_C_flanking.bed
