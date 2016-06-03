#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_A HLA_B HLA_C
	exit
fi

mkdir $2/$1_exon23_high_resolution_multi_ref


#prepare an original exon file to be filtered out
samtools view -b -L $2.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

#~/software/Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq
~/HLAssake/bam2fastq-1.1.0/bam2fastq $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -o $1#_exon23_high_resolution_multi_ref_mappedreads_target.fq

#relocate the relevant temporary files
mv $1_exon23_high_resolution_multi_ref_mappedreads_target.bam ./$2/$1_exon23_high_resolution_multi_ref

cp $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq ../
cp $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq ../
mv $1_*_exon23_high_resolution_multi_ref_mappedreads_target.fq ./$2/$1_exon23_high_resolution_multi_ref

cd ../
vim patient.fof <<end > /dev/null 2>&1
:%s/.*_1_.*/$1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq/g
:%s/.*_2_.*/$1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq/g
:wq
end

./HPTASRwgs_classI-II_onestep.sh
mv TASRhla* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv HLAminer* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv *ncbi.coord ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv formatdb.log ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
cp patient.fof ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
rm $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq
rm $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq

cd ./four_digit_prediction
