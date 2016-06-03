#!/bin/bash

declare string[2]
string[1]=$1
string[2]=$2

gaplocation1=`sed -n 1p gapinfo.dat | cut -d ' ' -f 1`
gaplocation2=`sed -n 1p gapinfo.dat | cut -d ' ' -f 2`
gaplocation3=`sed -n 2p gapinfo.dat | cut -d ' ' -f 1`
gaplocation4=`sed -n 2p gapinfo.dat | cut -d ' ' -f 2`

gaplocation1a=`expr ${gaplocation1} + 1`
gaplocation2a=`expr ${gaplocation2} + 1`
gaplocation3a=`expr ${gaplocation3} + 1`
gaplocation4a=`expr ${gaplocation4} + 1`

if [[ ${gaplocation2} != "" &&  ${gaplocation2} < 270 && ${gaplocation3} != "" &&  ${gaplocation3} > 270 ]]; then

	expr substr "${string[1]}" 1 ${gaplocation1} > temp
	read flagment1 < temp

	size2=`expr 271 - ${gaplocation2}`
	expr substr "${string[1]}" ${gaplocation2} ${size2} > temp
	read flagment2 < temp

	size3=`expr ${gaplocation3a} - 271`
	expr substr "${string[1]}" 271 ${size3} > temp
	read flagment3 < temp

	size4=`expr 547 - ${gaplocation4}`
	expr substr "${string[1]}" ${gaplocation4} ${size4} > temp
	read flagment4 < temp

        expr substr "${string[2]}" 1 ${gaplocation1} > temp
        read flagment5 < temp

        expr substr "${string[2]}" ${gaplocation2} ${size2} > temp
        read flagment6 < temp

        expr substr "${string[2]}" 271 ${size3} > temp
        read flagment7 < temp

	expr substr "${string[2]}" ${gaplocation4} ${size4} > temp
	read flagment8 < temp

	gapsize=`expr ${gaplocation2} - ${gaplocation1a}`
	expr substr "${string[1]}" ${gaplocation1a} ${gapsize} > temp
	read gapflagment < temp

	gap1size=`expr ${gaplocation4} - ${gaplocation3a}`
	expr substr "${string[1]}" ${gaplocation3a} ${gap1size} > temp
	read gap1flagment < temp

	allele1=`echo ${flagment1}${gapflagment}${flagment2}${flagment3}${gap1flagment}${flagment4}`
	allele2=`echo ${flagment1}${gapflagment}${flagment2}${flagment3}${gap1flagment}${flagment8}`
	allele3=`echo ${flagment1}${gapflagment}${flagment2}${flagment7}${gap1flagment}${flagment4}`
	allele4=`echo ${flagment1}${gapflagment}${flagment2}${flagment7}${gap1flagment}${flagment8}`
	allele5=`echo ${flagment1}${gapflagment}${flagment6}${flagment3}${gap1flagment}${flagment4}`
	allele6=`echo ${flagment1}${gapflagment}${flagment6}${flagment3}${gap1flagment}${flagment8}`
	allele7=`echo ${flagment1}${gapflagment}${flagment6}${flagment7}${gap1flagment}${flagment4}`
        allele8=`echo ${flagment1}${gapflagment}${flagment6}${flagment7}${gap1flagment}${flagment8}`
	allele9=`echo ${flagment5}${gapflagment}${flagment2}${flagment3}${gap1flagment}${flagment4}`
	allele10=`echo ${flagment5}${gapflagment}${flagment2}${flagment3}${gap1flagment}${flagment8}`
	allele11=`echo ${flagment5}${gapflagment}${flagment2}${flagment7}${gap1flagment}${flagment4}`
	allele12=`echo ${flagment5}${gapflagment}${flagment2}${flagment7}${gap1flagment}${flagment8}`
	allele13=`echo ${flagment5}${gapflagment}${flagment6}${flagment3}${gap1flagment}${flagment4}`
	allele14=`echo ${flagment5}${gapflagment}${flagment6}${flagment3}${gap1flagment}${flagment8}`
	allele15=`echo ${flagment5}${gapflagment}${flagment6}${flagment7}${gap1flagment}${flagment4}`
        allele16=`echo ${flagment5}${gapflagment}${flagment6}${flagment7}${gap1flagment}${flagment8}`

	a=`grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele1" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele2" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele3" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele4" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele5$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele5" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele6$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele6" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele7$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele7" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele8$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele8" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele9$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele9" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele10$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele10" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele11$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele11" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele12$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele12" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele13$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele13" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele14$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele14" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele15$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele15" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi

	a=`grep -B1 "$allele16$" ../../database/HLA-I_II_GEN.fasta`
	seq=`grep "$allele16" a.out`
	if [[ $seq = "" && $a != "" ]]; then
		echo -e "$a" >> a.out
	fi


fi


#report
rm report.out
echo -e "------------------------------------------------------------------" >> report.out
echo -e "HLAreporter version 1.01" >> report.out
echo -e "Report created" >> report.out
date >> report.out
echo -e "------------------------------------------------------------------" >> report.out
echo -e "Alleles detected" >> report.out
grep "\*" a.out >> report.out

#echo -e "\nAllele pair" >> report.out
#grep "\*" allele_pair.dat >> report.out

#gapsize=`expr \`sed -n 1p gapinfo.dat | cut -d ' ' -f 2\` - \`sed -n 1p gapinfo.dat | cut -d ' ' -f 1\`` 
#gap1size=`expr \`sed -n 2p gapinfo.dat | cut -d ' ' -f 2\` - \`sed -n 2p gapinfo.dat | cut -d ' ' -f 1\`` 
#echo -e "\nNon-polymorphic gap (potential intronic gap not shown):\n$gapsize bp\n$gap1size bp" >> report.out
#echo -e "Gap location:\n`sed -n 1p gapinfo.dat`\n`sed -n 2p gapinfo.dat`" >> report.out
