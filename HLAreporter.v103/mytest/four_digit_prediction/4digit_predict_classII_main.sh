#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_DRB1 HLA_DQB1 HLA_DRB3 HLA_DRB4 HLA_DRB5 
	exit
fi

mkdir $2
mkdir $2/$1_exon23_high_resolution_multi_ref

#prepare reads from flanking regions as reference information for exon read acquisition

samtools view -b -L $2_tempflanking.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.sam

#prepare an original exon file to be filtered out
samtools view -b -L $2.bed $1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target.bam
samtools view $1_exon23_high_resolution_multi_ref_mappedreads_target.bam > $1_exon23_high_resolution_multi_ref_mappedreads_target.sam
cp $1_exon23_high_resolution_multi_ref_mappedreads_target.sam $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam

#now start filtering out the exon file based on flanking read information

NUM=1
LINENUM=`wc -l $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.sam | sed 's/ /\t/' | cut -f 1`

while true; do
{
	STRING=`sed -n ${NUM}p $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.sam | cut -f 1,2`

	if [[ ${STRING} != "" ]]; then
		sed -i 's/'"${STRING}"'.*//g' $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam
	fi

	NUM=`expr ${NUM} + 1`

	if [ ${NUM} -gt ${LINENUM} ]; then
		break
	fi
}
done

#fix blank line
cp $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam temp.sam
sed '/^$/d' temp.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam
rm temp.sam

samtools view -bT exon23_high_resolution_multi_ref.fa $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam > $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.bam

#exclude reads with ambiguous mapping to other genes
if [[ $2 = "HLA_DRB1" || $2 = "HLA_DQB1" || $2 = "HLA_DPB1" || $2 = "HLA_DQA1" ]]; then
./purify_drb345.sh $1 $2
mv $1_ambiguousread_drb345.txt ./$2/$1_exon23_high_resolution_multi_ref
else
./purify_drdqdp.sh $1 $2
mv $1_ambiguousread_drdqdp.txt ./$2/$1_exon23_high_resolution_multi_ref
fi

../../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq

#relocate the relevant temporary files
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.sam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.bam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.bam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target_exon.sam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target.bam ./$2/$1_exon23_high_resolution_multi_ref
mv $1_exon23_high_resolution_multi_ref_mappedreads_target.sam ./$2/$1_exon23_high_resolution_multi_ref

cp $1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq ../
cp $1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq ../
mv $1_*_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq ./$2/$1_exon23_high_resolution_multi_ref

cd ../
#vim patient.fof <<end > /dev/null 2>&1
#:%s/.*_1_.*/$1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq/g
#:%s/.*_2_.*/$1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq/g
#:wq
#end

sed -i 's/.*_1_.*/'"$1"'_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq/g' patient.fof
sed -i 's/.*_2_.*/'"$1"'_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq/g' patient.fof

./HPTASRwgs_classI-II_mainII.sh
mv TASRhla* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv HLAminer* ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv *ncbi.coord ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
mv formatdb.log ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
cp patient.fof ./four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref
rm $1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq
rm $1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq

cd ./four_digit_prediction
