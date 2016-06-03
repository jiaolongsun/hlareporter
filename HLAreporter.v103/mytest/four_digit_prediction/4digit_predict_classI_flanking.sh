#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_A HLA_B HLA_C
	exit
fi

mkdir $2/$1_exon23_high_resolution_multi_ref_flanking

./format_classI_flanking_bed.sh $1 $2

cd ./$2/$1_exon23_high_resolution_multi_ref

#remove the existing flanking miner database before generating a latest one
rm ../../../../database/HLA-I_II_GEN_FLANKING.fasta

#get the main step results
MAINSET=`grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | sed ':label;N;s/\n/\ /;b label'`

for PATTERN in ${MAINSET}; do

#there will be a special char star in the variable so that the star char needs to be converted to be matched
FIXEDPATTERN=`echo ${PATTERN} | sed 's/\*/\\\*/g'`

grep -A1 ">$FIXEDPATTERN$" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_FLANKING.fasta

done


#fix the flanking database with NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
cd ../../../../database
cat HLA-I_II_GEN_FLANKING.fasta > TEMP.fasta
rm HLA-I_II_GEN_FLANKING.fasta

cat TEMP.fasta | while read STRING
do

	LENGTH=`expr length $STRING`

	if [ $LENGTH -gt 270 ]; then

		SUBSTR1=`expr substr "$STRING" 1 "$POS1"`
		SUBSTR2=`expr substr "$STRING" "$POS2" 1000`
		SUBSTR3="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
		NEWSTR=${SUBSTR1}${SUBSTR3}${SUBSTR2}
		echo $NEWSTR >> HLA-I_II_GEN_FLANKING.fasta

	else

		echo $STRING >> HLA-I_II_GEN_FLANKING.fasta

		#there will be a special char star in the variable so that the star char needs to be converted to be matched
		FIXEDSTRING=`echo ${STRING} | sed 's/\*/\\\*/g'`

		ALLELEINFO=`grep "$FIXEDSTRING" ../mytest/four_digit_prediction/special_allele_set.dat`
		if [[ ! $ALLELEINFO ]]; then

			POS1=270
			POS2=271

		else
			
			POS1=`echo $ALLELEINFO | cut -d ' ' -f 2`
			POS2=`expr $POS1 + 1`

		fi
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

samtools view $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.sam

#filter out ambiguous reads mapping to other genes
./purify_efgh_flanking.sh $1 $2

#filter out ambiguous reads mapping within abc genes
./purify_abc_flanking.sh $1 $2

#start predicting

mv $1_*_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ./$2/$1_exon23_high_resolution_multi_ref_flanking
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.bam ./$2/$1_exon23_high_resolution_multi_ref_flanking
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.sam ./$2/$1_exon23_high_resolution_multi_ref_flanking

cd ./$2/$1_exon23_high_resolution_multi_ref_flanking

cp $1_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ../../../
cp $1_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq ../../../


cd ../../../
sed -i 's/.*_1_.*/'"$1"'_1_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq/g' patient.fof
sed -i 's/.*_2_.*/'"$1"'_2_exon23_high_resolution_multi_ref_mappedreads_target_flanking.fq/g' patient.fof

./HPTASRwgs_classI-II_flanking.sh
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

#./report_hla.sh $1 $2

#./phase_analyser.sh $1 $2
cp ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR_round1.log
cp ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla.contigs ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla_round1.contigs
cp ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla90.contigs ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla90_round1.contigs
cp ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla.readposition ./$2/$1_exon23_high_resolution_multi_ref_flanking/TASRhla_round1.readposition
./4digit_predict_classI_flanking_recursive.sh $1 $2

./4digit_predict_classI_flanking_recursive1.sh $1 $2
