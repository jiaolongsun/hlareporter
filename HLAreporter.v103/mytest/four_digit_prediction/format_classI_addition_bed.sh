#!/bin/bash

echo "Formatting classI addition bed..."

#calculate fastq read length

LENGTH=`cat $2/$1_exon23_high_resolution_multi_ref/$1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq | awk '{if(NR%4==2) {print length($1); exit 0;}}'`

if [ ! $LENGTH ]; then
	echo "NULL read length! program teminated!"
	exit
fi

echo "Fastq read length: ${LENGTH}"

#calculate each addition bed range

RANGE1=1
RANGE2=`expr 148 - ${LENGTH} + 9`
RANGE3=`expr 148 + ${LENGTH} - 9`
RANGE4=296

#update classI addition bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_A_addition.bed
sed -i 's/\t.*\t296/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_A_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_B_addition.bed
sed -i 's/\t.*\t296/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_B_addition.bed

sed -i 's/\t1\t.*/\t'"$RANGE1"'\t'"$RANGE2"'/g' HLA_C_addition.bed
sed -i 's/\t.*\t296/\t'"$RANGE3"'\t'"$RANGE4"'/g' HLA_C_addition.bed
