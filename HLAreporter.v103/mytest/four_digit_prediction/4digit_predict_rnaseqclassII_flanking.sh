#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_DRB1 HLA_DQB1 HLA_DRB3 HLA_DRB4 HLA_DRB5
	exit
else
	if [[ $2 -eq "HLA_DRB1" || $2 -eq "HLA_DRB3" || $2 -eq "HLA_DRB4" || $2 -eq "HLA_DRB5" || $2 -eq "HLA_DQB1" ]]; then
		echo "Detecting $2"
	else
		exit
	fi
fi

mkdir $2/$1_exon23_high_resolution_multi_ref_flanking

./format_classII_flanking_bed.sh $1 $2

cd ./$2/$1_exon23_high_resolution_multi_ref

#remove the existing flanking miner database before generating a latest one
rm ../../../../database/HLA-I_II_GEN_FLANKING.fasta

#get the main step results
MAINSET=`grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | sed ':label;N;s/\n/\ /;b label'`

for PATTERN in ${MAINSET}; do

#there will be a special char star in the variable so that the star char needs to be converted to be matched
FIXEDPATTERN=`echo ${PATTERN} | sed 's/\*/\\\*/g'`

grep -A1 "$FIXEDPATTERN$" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_FLANKING.fasta

done


#fix the flanking database with NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
cd ../../../../database
cat HLA-I_II_GEN_FLANKING.fasta > TEMP.fasta
rm HLA-I_II_GEN_FLANKING.fasta

cat TEMP.fasta | while read STRING
do

	LENGTH=`expr length $STRING`

	if [ $LENGTH -gt 500 ]; then

		SUBSTR1=`expr substr "$STRING" 1 270`
		SUBSTR2=`expr substr "$STRING" 271 1000`
		SUBSTR3="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
		NEWSTR=${SUBSTR1}${SUBSTR3}${SUBSTR2}
		echo $NEWSTR >> HLA-I_II_GEN_FLANKING.fasta

	else

		echo $STRING >> HLA-I_II_GEN_FLANKING.fasta

	fi

done

rm TEMP.fasta


#compile the flanking database

../bin/formatdb -p F -i HLA-I_II_GEN_FLANKING.fasta
bwa index -a is HLA-I_II_GEN_FLANKING.fasta

if [[ ${PATTERN} != "" ]]; then
	echo "HLA-I_II_GEN_FLANKING.fasta established!"
else
	echo "HLA-I_II_GEN_FLANKING.fasta is EMPTY!"
fi

cd ../mytest/four_digit_prediction

#obtain the flanking reads based on bed file info
samtools view -b -L $2_flanking.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam

#start predicting
../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq

mv $1_*_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ./$2/$1_exon23_high_resolution_multi_ref_flanking
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam ./$2/$1_exon23_high_resolution_multi_ref_flanking


cd ./$2/$1_exon23_high_resolution_multi_ref_flanking

cp $1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ../../../
cp $1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ../../../


cd ../../../
vim patient.fof <<end > /dev/null 2>&1
:%s/.*_1_.*/$1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq/g
:%s/.*_2_.*/$1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq/g
:wq
end

./HPTASRwgs_classI-II_rnaseq_flankingII.sh
mv TASRhla* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking
mv HLAminer* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking
mv *ncbi.coord ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking
mv formatdb.log ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking
cp patient.fof ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking
rm $1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq
rm $1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq

#show flanking results
echo "Done examining the main exons"
grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log


cd ./four_digit_prediction/

./report_hla.sh $1 $2

./phase_analyser.sh $1 $2

