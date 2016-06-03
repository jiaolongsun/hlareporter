#!/bin/bash

if [ $# -ne 4 ]; then
	echo -e "HLAreporter version 1.03"
	echo "Usage: command sample_name gene_name fq1 fq2"
	echo "Gene: HLA_A HLA_B HLA_C HLA_DRB1 HLA_DRB3 HLA_DRB4 HLA_DRB5 HLA_DQB1 HLA_DPB1 HLA_DQA1"
	exit
fi

if [[ $2 != "HLA_A" && $2 != "HLA_B" && $2 != "HLA_C" && $2 != "HLA_DRB1" && $2 != "HLA_DRB3" && $2 != "HLA_DRB4" && $2 != "HLA_DRB5" && $2 != "HLA_DQB1" && $2 != "HLA_DPB1" && $2 != "HLA_DQA1" ]]; then
	echo "Gene entered out of scope"
	exit
fi

cd ./mytest/four_digit_prediction

#echo "Mapping required? Y/N (If you have done this step, please select N)"
#read option
#if [[ $option = "Y" || $option = "y" ]]; then
if [ ! -f "../../temp/$1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam" ]; then
	mv ../../$3 ./
	mv ../../$4 ./ 

	./4digit_map_HLA.sh $3 $4 $1
	
	mv ./$3 ../../
	mv ./$4 ../../
	mv $1_*_exon23_high_resolution_multi_ref.sai ../../temp/

else
	cp -r -p ../../temp/$1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam ./
fi

if [[ $2 = "HLA_A" || $2 = "HLA_B" || $2 = "HLA_C" ]]; then
	./4digit_predict_classI_ssake.sh $1 $2
	./4digit_predict_classI_flanking.sh $1 $2
	sed -i 's/>//g' report.out
 	cat report.out
	echo -e "\nHLA report created in report.out"

	echo -e "\nG group member detection required? Y/N (If you observe ambiguity, please select N)"
	read -t 15 option
	if [[ $option = "Y" || $option = "y"  ]]; then
        	./4digit_predict_classI_addition.sh $1 $2
		echo -e "\nHLA report created in report.out"
	fi

fi

if [[ $2 = "HLA_DRB1" || $2 = "HLA_DRB3" || $2 = "HLA_DRB4" || $2 = "HLA_DRB5" || $2 = "HLA_DQB1" || $2 = "HLA_DPB1" || $2 = "HLA_DQA1" ]]; then
	./4digit_predict_classII_main.sh $1 $2
	./4digit_predict_classII_flanking.sh $1 $2
        sed -i 's/>//g' report.out
        cat report.out
        echo -e "\nHLA report created in report.out"

        echo -e "\nG group member detection required? Y/N (If you observe ambiguity, please select N)"
        read -t 15 option
        if [[ $option = "Y" || $option = "y" ]]; then
                ./4digit_predict_classII_addition.sh $1 $2
		echo -e "\nHLA report created in report.out"
        fi

fi

cp report.out ../../$1_$2_report.out

cd ../../qualityprofile

echo -e "Generating HLA data quality profile..."
read -t 3 option
./bam_point_cov.sh $1 $2
echo -e "HLA data quality profile:" >> ../$1_$2_report.out
cat $2/$1/$1_$2_covsummary.out >> ../$1_$2_report.out

DEPTH=`grep "20xcov" $2/$1/$1_$2_covsummary.out | cut -f 4`
if [ $DEPTH -lt 98 ]; then
	echo -e "\nWARNING: Low data quality" >> ../$1_$2_report.out
fi

cd ../

./HLAfreq.sh $1 $2
sed -i 's/version 1.01/version 1.03/g' $1_$2_report.out

mv $1_$2_report.out results
mv ./mytest/four_digit_prediction/$1_exon23_high_resolution_multi_ref_mappedreads_sorted.bam* ./temp

