#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_A HLA_B HLA_C
	exit
fi

./format_classI_addition_bed.sh $1 $2

#remove the existing addition miner database before generating a latest one
rm ../../database/HLA-I_II_GEN_ADDITION.fasta

#get the g group results from the flanking results
GGROUPSET=`grep -E '\*.*G' ./report.out | sed ':label;N;s/\n/\ /;b label'`
PREGGROUP=""

for FLANKINGRESULT in ${GGROUPSET}; do

#do not handle if it has been handled before
THIS=`echo $FLANKINGRESULT | awk -F \* '{print $2}'`
CHECKFLANKING=`echo $PREGGROUP | grep "$THIS"`
if [[ $CHECKFLANKING != "" ]]; then
continue
fi
PREGGROUP=`echo -e "${PREGGROUP} ${FLANKINGRESULT}"`

#get the gene name
GENE=`echo $FLANKINGRESULT | awk -F \* '{print $1}'`

#get the group members, replace / with space

#first need to handle the special char star in the variable
FIXEDFLANKING=`echo ${FLANKINGRESULT} | sed 's/\*/\\\*/g'`

#now start grep process and replace char / with char space
SET=`grep "$FIXEDFLANKING" all_g_groups.txt | awk '{print $2}' | sed 's/\//\ /g'`

#start exploring allele sequence and save them to miner's ADDITION.fasta database
for TYPE in ${SET}; do

grep -A1 "$GENE\*$TYPE$" HLA-I_II_GEN_ADDITIONSET.fasta >> ../../database/HLA-I_II_GEN_ADDITION.fasta

done

done

#compile the addition database
cd ../../database
../bin/formatdb -p F -i HLA-I_II_GEN_ADDITION.fasta
bwa index -a is HLA-I_II_GEN_ADDITION.fasta

if [[ ${GGROUPSET} != "" ]]; then
	echo "HLA-I_II_GEN_ADDITION.fasta established!"
else
	echo "HLA-I_II_GEN_ADDITION.fasta is EMPTY!"
fi

cd ../mytest/four_digit_prediction

mkdir $2/$1_exon23_high_resolution_multi_ref_addition

samtools view -b -L $2_addition.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_addition_mappedreads_target.bam

../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_addition_mappedreads_target.bam -fq1 $1_1_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq -fq2 $1_2_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq

#relocate the relevant temporary files
mv $1_exon23_high_resolution_multi_ref_addition_mappedreads_target.bam ./$2/$1_exon23_high_resolution_multi_ref_addition

cp $1_1_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq ../
cp $1_2_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq ../
mv $1_*_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq ./$2/$1_exon23_high_resolution_multi_ref_addition

cd ../
#vim patient.fof <<end > /dev/null 2>&1
#:%s/.*_1_.*/$1_1_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq/g
#:%s/.*_2_.*/$1_2_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq/g
#:wq
#end

sed -i 's/.*_1_.*/'"$1"'_1_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq/g' patient.fof
sed -i 's/.*_2_.*/'"$1"'_2_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq/g' patient.fof

./HPTASRwgs_classI-II_addition.sh
mv TASRhla* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition
mv HLAminer* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition
mv *ncbi.coord ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition
mv formatdb.log ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition
cp patient.fof ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition
rm ./$1_1_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq
rm ./$1_2_exon23_high_resolution_multi_ref_addition_mappedreads_target.fq

#show additional results
echo "Done examining the additional exons"
grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref_addition/HLAminer_HPTASR.log

cd ./four_digit_prediction

./report_ggroup_hla.sh $1 $2
./report_ggroup_hla.sh $1 $2 >> report.out

