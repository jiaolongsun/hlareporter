#!/bin/bash

echo -e "------------------------------------------------------------------"
echo -e "HLAreporter version 1.01"
echo -e Report created 
date
echo -e "------------------------------------------------------------------"

cd ./$2/$1_exon23_high_resolution_multi_ref_flanking

#only report at most 2 alleles according to score value
ALLELECOUNT=0

PRESCORE=0

#get the main step results
MAINSET=`grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | sed ':label;N;s/\n/\ /;b label'`

for PATTERN in ${MAINSET}; do

	#there will be a special char star in the variable so that the star char needs to be converted to be matched
	FIXEDPATTERN=`echo ${PATTERN} | sed 's/\*/\\\*/g'`

	#get the flanking score of this pattern
	SCORE=`grep -E ".*,.*,$FIXEDPATTERN,.*,.*,.*,.*,.*," ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5`

	#convert to int
	FSCORE=`echo ${SCORE%.*}`

	#find the main score of this pattern
	SCORE=`grep -E ".*,.*,$FIXEDPATTERN,.*,.*,.*,.*,.*," ../$1_exon23_high_resolution_multi_ref/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5 `

	#convert to int
	MSCORE=`echo ${SCORE%.*}`


	#calculate the total score

	#first ignore the first line of discription words
	if [[ ! $FSCORE || ! $MSCORE ]]; then

		continue
	fi

	SCORE=`expr $FSCORE + $MSCORE`

	#output, [e] option is needed for TAB
	if [[ $SCORE -eq $PRESCORE ]]; then

		echo -e "\\t${PATTERN}\\n"
		
	else
		
		if [[ $ALLELECOUNT -eq 2 && $2 != "HLA_A" && $2 != "HLA_B" && $2 != "HLA_C" ]]; then
			exit
		fi
		echo -e "Allele\\n\\t${PATTERN}\\n"
		PRESCORE=$SCORE
		ALLELECOUNT=`expr ${ALLELECOUNT} + 1`

fi

done
