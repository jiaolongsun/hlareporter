#!/bin/bash

echo -e "------------------------------------------------------------------"
echo -e "G group report created"
echo -e "Alleles with the same minor exon will appear in the same group"
date
echo -e "------------------------------------------------------------------"

cd ./$2/$1_exon23_high_resolution_multi_ref_addition

PRESCORE=0

#get the main step results
MAINSET=`grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | sed ':label;N;s/\n/\ /;b label'`

for PATTERN in ${MAINSET}; do

	#there will be a special char star in the variable so that the star char needs to be converted to be matched
	FIXEDPATTERN=`echo ${PATTERN} | sed 's/\*/\\\*/g'`

	#get the flanking score of this pattern
	SCORE=`grep -E ".*,.*,$FIXEDPATTERN,.*,.*,.*,.*,.*," ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5`

	#convert to int
	INTSCORE=`echo ${SCORE%.*}`

	#first ignore the first line of discription words
	if [[ ! $INTSCORE  ]]; then

		continue
	fi

	#output, [e] option is needed for TAB
	if [[ $INTSCORE -eq $PRESCORE ]]; then

		echo -e "\\t${PATTERN}\\n"
		
	else
		
		echo -e "G group allele\\n\\t${PATTERN}\\n"
		PRESCORE=$INTSCORE

fi

done
