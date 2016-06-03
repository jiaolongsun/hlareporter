#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command argu1 argu2
	exit
fi

if [[ $2 = "HLA_A" || $2 = "HLA_B" || $2 = "HLA_C" || $2 = "HLA_DRB1" || $2 = "HLA_DQB1" || $2 = "HLA_DPB1" || $2 = "HLA_DQA1" ]]; then

	#output a format
	echo -e "------------------------------------------------------------------" >> $1_$2_report.out
	echo -e "HLA allele frequency" >> $1_$2_report.out
	echo -e "Four populations in Europe China Japan Africa are shown" >> $1_$2_report.out
	echo -e "By allele frequency net database (www.allelefrequencies.net)" >> $1_$2_report.out
	echo -e "------------------------------------------------------------------" >> $1_$2_report.out
	echo -e "[Allele]\t[EUR]\t[CHN]\t[JPN]\t[AFR]" >> $1_$2_report.out

	SET=`grep -E '\*' ./$1_$2_report.out | sed ':label;N;s/\n/\ /;b label'`
	PRESET=""

	for RESULT in ${SET}; do

	#normalize to 4 digit
	GENE=`echo $RESULT | awk -F \* '{print $1}'`
	DIGIT2=`echo $RESULT | awk -F \* '{print $2}' | awk -F \: '{print $1}'`
	DIGIT4=`echo $RESULT | awk -F \* '{print $2}' | awk -F \: '{print $2}'`
	THIS=`echo -e "${DIGIT2}:${DIGIT4}"`

	#do not handle if it has been handled before
	CHECK=`echo $PRESET | grep "$THIS"`
	if [[ $CHECK != "" ]]; then
		continue
	fi
	PRESET=`echo -e "${PRESET} ${GENE}*${THIS}"`

	#obtain the frequency
	ALLELE=`echo -e "${GENE}*${THIS}"`
	#first need to handle the special char star in the variable
	FIXEDALLELE=`echo ${ALLELE} | sed 's/\*/\\\*/g'`

	echo -n "$ALLELE	" >> $1_$2_report.out

	EURFREQ=`grep "$FIXEDALLELE	" ./freq/US_NMDP_European_Caucasian.txt | cut -f 4`
	if [[ $EURFREQ != "" ]]; then
	echo -n "$EURFREQ	" >> $1_$2_report.out
	else
	echo -n "0.0000	" >> $1_$2_report.out
	fi

	CHNFREQ=`grep "$FIXEDALLELE	" ./freq/China_Jiangsu_Han.txt | cut -f 4`
	if [[ $CHNFREQ != "" ]]; then
	echo -n "$CHNFREQ	" >> $1_$2_report.out
	else
	echo -n "0.0000	" >> $1_$2_report.out
	fi

	JPNFREQ=`grep "$FIXEDALLELE	" ./freq/Japan_pop_16.txt | cut -f 4`
	if [[ $JPNFREQ != "" ]]; then
	echo -n "$JPNFREQ	" >> $1_$2_report.out
	else
	echo -n "0.0000	" >> $1_$2_report.out
	fi

	#this time output char \n
	AFRFREQ=`grep "$FIXEDALLELE	" ./freq/South_Africa_Black.txt | cut -f 4`
	if [[ $AFRFREQ != "" ]]; then
	echo "$AFRFREQ" >> $1_$2_report.out
	else
	echo "0.0000" >> $1_$2_report.out
	fi

	done
	
else
	echo -e "------------------------------------------------------------------" >> $1_$2_report.out
	echo -e "HLA allele frequency not available" >> $1_$2_report.out

fi

