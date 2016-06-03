#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_A HLA_B HLA_C
	exit
fi

mkdir $2
mkdir $2/$1_exon23_high_resolution_multi_ref


#prepare an original exon file to be filtered out
samtools view -b -L $2.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam

samtools view $1_exon23_high_resolution_multi_ref_mappedreads_target.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target.sam

#exclude reads with ambiguous mapping to other genes
./purify_efgh.sh $1 $2

#exclude reads with ambiguous mapping within abc genes
./purify_abc.sh $1 $2

#../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq
../../bam2fastq-1.1.0/bam2fastq $1_exon23_high_resolution_multi_ref_mappedreads_target.bam -o $1#_exon23_high_resolution_multi_ref_mappedreads_target.fq


#relocate the relevant temporary files
mv $1_exon23_high_resolution_multi_ref_mappedreads_target.sam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target.bam ./$2/$1_exon23_high_resolution_multi_ref

cp $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq ../
cp $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq ../
mv $1_*_exon23_high_resolution_multi_ref_mappedreads_target.fq ./$2/$1_exon23_high_resolution_multi_ref

cd ../

#replace patient.fof content
sed -i 's/.*_1_.*/'"$1"'_1_exon23_high_resolution_multi_ref_mappedreads_target.fq/g' patient.fof
sed -i 's/.*_2_.*/'"$1"'_2_exon23_high_resolution_multi_ref_mappedreads_target.fq/g' patient.fof

#creat corresponding database
cd ./four_digit_prediction
./establish_ssake_db.sh $1 $2
cd ../

./HPTASRwgs_classI-II_ssake.sh
mv paired.fa ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv paired.fa.ssake* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv TASRhla* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv HLAminer* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv *ncbi.coord ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv formatdb.log ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
cp patient.fof ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
rm $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq
rm $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq

cd ./four_digit_prediction
