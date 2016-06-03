#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command input_filename gene
	echo Gene option: HLA_A HLA_B HLA_C
	exit
fi

cd ./$2/$1_exon23_high_resolution_multi_ref

#remove the existing flanking miner database before generating a latest one
rm ../../../../database/HLA-I_II_GEN_SSAKE.fasta

#get the main step results

MAINSET=`echo "$2" | sed 's/\_/\t/g' | cut -f 2`

grep -A1 ">$MAINSET\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta

#recruit E,F,G,H,J... when detecting A,B,C
if [[ ${MAINSET} = "A" || ${MAINSET} = "B" || ${MAINSET} = "C" ]]; then

	grep -A1 ">E\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta
	grep -A1 ">F\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta
	grep -A1 ">G\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta
	grep -A1 ">H\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta
	grep -A1 ">J\*" ../../../../database/HLA-I_II_GEN.fasta >> ../../../../database/HLA-I_II_GEN_SSAKE.fasta
fi

#fix the flanking database with NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
cd ../../../../database
cat HLA-I_II_GEN_SSAKE.fasta > TEMP.fasta
rm HLA-I_II_GEN_SSAKE.fasta

cat TEMP.fasta | while read STRING
do

	LENGTH=`expr length $STRING`

	if [ $LENGTH -gt 270 ]; then

		SUBSTR1=`expr substr "$STRING" 1 "$POS1"`
		SUBSTR2=`expr substr "$STRING" "$POS2" 1000`
		SUBSTR3="NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"
		NEWSTR=${SUBSTR1}${SUBSTR3}${SUBSTR2}
		echo $NEWSTR >> HLA-I_II_GEN_SSAKE.fasta

	else

		echo $STRING >> HLA-I_II_GEN_SSAKE.fasta

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


#compile the SSAKE database

../bin/formatdb -p F -i HLA-I_II_GEN_SSAKE.fasta
bwa index -a is HLA-I_II_GEN_SSAKE.fasta

if [[ ${MAINSET} != "" ]]; then
	echo "HLA-I_II_GEN_SSAKE.fasta established!"
else
	echo "HLA-I_II_GEN_SSAKE.fasta is EMPTY!"
fi

cd ../mytest/four_digit_prediction
