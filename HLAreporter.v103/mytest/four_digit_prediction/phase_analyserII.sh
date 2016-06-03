#!/bin/bash

#generate main report sequences

cd ./$2/$1_exon23_high_resolution_multi_ref_flanking

#remove the existing flanking miner database before generating a latest one
rm ../../allele_sequences.dat

#get the main step results
MAINSET=`grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | sed ':label;N;s/\n/\ /;b label'`

for PATTERN in ${MAINSET}; do

#there will be a special char star in the variable so that the star char needs to be converted to be matched
FIXEDPATTERN=`echo ${PATTERN} | sed 's/\*/\\\*/g'`

grep -A1 "$FIXEDPATTERN$" ../../../../database/HLA-I_II_GEN.fasta >> ../../allele_sequences.dat

done

cd ../../

#check if the number of line is OK

linecount=`wc -l allele_sequences.dat | sed 's/ /\t/' | cut -f 1`

if [[ linecount -ne 6 && linecount -ne 4 ]]; then

	#echo "The number of line is not expected"
	exit

else

#start reading each sequence and save it

num=1
declare string[3]
declare name[3]

#cat input.txt | while read str; this pipeline will create a subshell, where any changes in it cannot be accessed by outside
#now get all the info from a file, and assign an id 3 to describe this info, then read things from this info

exec 3<> allele_sequences.dat
while read str <&3
do

	length=`expr length $str`

	if [ $length -lt 20 ]; then

		name[$num]=${str}
		#echo ${name[$num]}
		
	else
		
		string[$num]=${str}
		#echo ${string[$num]}
		num=`expr ${num} + 1`

	fi

done
exec 3>&-
#above close the info discription space 3

#calculate fastq read length

LENGTH=`cat $2/$1_exon23_high_resolution_multi_ref/$1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq | awk '{if(NR%4==2) {print length($1); exit 0;}}'`

if [ ! $LENGTH ]; then
	echo "NULL read length! program teminated!"
	exit
fi

##############################################################
#echo "Analyse potential phase of 2-allele case"

if [ $linecount = 4 ]; then

Location=1
Prelocation=1
Gapsize=0
Gap=0
Gapstart=0
Gapend=0

while true; do
{
	#get one nucleotide at a time
	substr1=${string[1]:${Location}:1}
	substr2=${string[2]:${Location}:1}

	if [[ ${substr1} = "" ]]; then
	
		break
	
	fi

	if [[ ${substr1} = ${substr2} ]]; then
	
		Location=`expr ${Location} + 1`

	else

		Gapsize=`expr ${Location} - ${Prelocation}`
		if [[ $Gapsize -gt $LENGTH && $Prelocation != 1 ]]; then
			Gap=$Gapsize
			Gapstart=`expr ${Prelocation} + 1`
			Gapend=`expr ${Location} + 1`
		fi
		Prelocation=$Location
		Location=`expr ${Location} + 1`

	fi

}

done

#check if it is just about ambiguity with the 2 alleles

#obtain 2 Scores starting from row 2, then compare 2 Scores
grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5 > temprecord
Score1=`sed -n 2p temprecord`
Score2=`sed -n 3p temprecord`

#if the three candidates are irrelevant to phase, terminate the analysis
if [[ $Score1 = $Score2 || $Gap = 0 ]]; then
	
	exit

fi

if [[ $Gap -ne 0 ]]; then

	##############################################################
	expr substr "${string[1]}" 1 ${Gapstart} > temp
	read flagment1 < temp

	size2=`expr 271 - ${Gapend}`
	expr substr "${string[1]}" ${Gapend} ${size2} > temp
	read flagment2 < temp

        expr substr "${string[2]}" 1 ${Gapstart} > temp
        read flagment3 < temp

        expr substr "${string[2]}" ${Gapend} ${size2} > temp
        read flagment4 < temp

	Gapsize=`expr ${Gapend} - ${Gapstart} - 1`
	Gapstarta=`expr ${Gapstart} + 1`
	expr substr "${string[1]}" ${Gapstarta} ${Gapsize} > temp
	read Gapflagment < temp

	allele1=`echo ${flagment1}${Gapflagment}${flagment2}`
	allele2=`echo ${flagment1}${Gapflagment}${flagment4}`
	allele3=`echo ${flagment3}${Gapflagment}${flagment2}`
	allele4=`echo ${flagment3}${Gapflagment}${flagment4}`

        grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta > a.out
	
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
	
	#no need to output if allele count is still 2 after phase analysis
	phasecount=`grep "\*" a.out | wc -l`	
	if [ $phasecount -gt 2 ]; then
	echo -e "Alleles detected" >> report.out
	grep "\*" a.out | cut -d '>' -f 2 >> report.out

        echo -e "\nGap Location:\n$Gapstart\t$Gapend" >> report.out
        echo -e "\nNon-polymorphic Gap:\n$Gap bp" >> report.out
	fi
	##############################################################
	

fi

#after phase analysis of 2-allele case, exit here
exit

fi
##############################################################


#start comparing each nucleotide

location=1
prelocation=1
gapsize=0
gap=0
gapstart=0
gapend=0
flag1=0
flag2=0
flag3=0

while true; do
{
	#get one nucleotide at a time
	substr1=${string[1]:${location}:1}
	substr2=${string[2]:${location}:1}
	substr3=${string[3]:${location}:1}

	if [[ ${substr1} = "" ]]; then
	
		break
	
	fi

	if [[ ${substr1} = ${substr2} && ${substr1} = ${substr3} ]]; then

		location=`expr ${location} + 1`

	elif [[ ${substr1} = ${substr2} ]]; then 

		flag3=1
		gapsize=`expr ${location} - ${prelocation}`
		if [[ $gapsize -gt $LENGTH ]]; then
			gap=$gapsize
			gapstart=`expr ${prelocation} + 1`
			gapend=`expr ${location} + 1`
		fi
		prelocation=$location
		location=`expr ${location} + 1`

	elif [[ ${substr1} = ${substr3} ]]; then 

		flag2=1
		gapsize=`expr ${location} - ${prelocation}`
		if [[ $gapsize -gt $LENGTH ]]; then
			gap=$gapsize
			gapstart=`expr ${prelocation} + 1`
			gapend=`expr ${location} + 1`
		fi
		prelocation=$location
		location=`expr ${location} + 1`

	elif [[ ${substr2} = ${substr3} ]]; then 

		flag1=1
		gapsize=`expr ${location} - ${prelocation}`
		if [[ $gapsize -gt $LENGTH ]]; then
			gap=$gapsize
			gapstart=`expr ${prelocation} + 1`
			gapend=`expr ${location} + 1`
		fi
		prelocation=$location
		location=`expr ${location} + 1`

	else

		flag1=1
		flag2=1
		flag3=1
		location=`expr ${location} + 1`
		#no need to calculate gap size as we do not consider this situation

	fi
}

done


#start printing the results

#obtain 3 scores starting from row 2, then compare 3 scores
grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5 > temprecord
score1=`sed -n 2p temprecord`
score2=`sed -n 3p temprecord`
score3=`sed -n 4p temprecord`

#if the three candidates are irrelevant to phase, terminate the analysis
if [[ $flag1 = 1 && $flag2 = 1 && $flag3 = 1 || $score1 = $score2 || $score2 = $score3 || $score1 = $score3 || $gap = 0 ]]; then
	
	exit

fi

#echo "Seperating:"

echo -e "------------------------------------------------------------------" > report.out
echo -e "HLAreporter version 1.01" >> report.out
echo -e Report created >> report.out
date >> report.out
echo -e "------------------------------------------------------------------" >> report.out

if [ $flag1 = 0 ]; then

	echo -e "Allele\\n\\t${name[2]}\\n" >> report.out
	echo -e "Allele\\n\\t${name[3]}\\n" >> report.out
	echo -e "Phase\\n\\t${name[1]}\\n" >> report.out

	##############################################################
	expr substr "${string[2]}" 1 ${gapstart} > temp
	read flagment1 < temp

	size2=`expr 271 - ${gapend}`
	expr substr "${string[2]}" ${gapend} ${size2} > temp
	read flagment2 < temp

        expr substr "${string[3]}" 1 ${gapstart} > temp
        read flagment3 < temp

        expr substr "${string[3]}" ${gapend} ${size2} > temp
        read flagment4 < temp

	gapsize=`expr ${gapend} - ${gapstart} - 1`
	gapstarta=`expr ${gapstart} + 1`
	expr substr "${string[2]}" ${gapstarta} ${gapsize} > temp
	read gapflagment < temp

	allele1=`echo ${flagment1}${gapflagment}${flagment2}`
	allele2=`echo ${flagment1}${gapflagment}${flagment4}`
	allele3=`echo ${flagment3}${gapflagment}${flagment2}`
	allele4=`echo ${flagment3}${gapflagment}${flagment4}`

        grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta > a.out
        grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out

	echo -e "Alleles detected" >> report.out
	grep "\*" a.out | cut -d '>' -f 2 >> report.out
	##############################################################

	if [[ $gap -ne 0 ]]; then
		echo -e "\nGap location:\n$gapstart\t$gapend" >> report.out
		echo -e "\nNon-polymorphic gap:\n$gap bp" >> report.out
	fi

fi

if [ $flag2 = 0 ]; then

	echo -e "Allele\\n\\t${name[1]}\\n" >> report.out
	echo -e "Allele\\n\\t${name[3]}\\n" >> report.out
	echo -e "Phase\\n\\t${name[2]}\\n" >> report.out

	##############################################################
	expr substr "${string[1]}" 1 ${gapstart} > temp
	read flagment1 < temp

	size2=`expr 271 - ${gapend}`
	expr substr "${string[1]}" ${gapend} ${size2} > temp
	read flagment2 < temp

        expr substr "${string[3]}" 1 ${gapstart} > temp
        read flagment3 < temp

        expr substr "${string[3]}" ${gapend} ${size2} > temp
        read flagment4 < temp

	gapsize=`expr ${gapend} - ${gapstart} - 1`
	gapstarta=`expr ${gapstart} + 1`
	expr substr "${string[1]}" ${gapstarta} ${gapsize} > temp
	read gapflagment < temp

	allele1=`echo ${flagment1}${gapflagment}${flagment2}`
	allele2=`echo ${flagment1}${gapflagment}${flagment4}`
	allele3=`echo ${flagment3}${gapflagment}${flagment2}`
	allele4=`echo ${flagment3}${gapflagment}${flagment4}`

        grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta > a.out
        grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out

	echo -e "Alleles detected" >> report.out
	grep "\*" a.out | cut -d '>' -f 2 >> report.out
	##############################################################

	if [[ $gap -ne 0 ]]; then
		echo -e "\nGap location:\n$gapstart\t$gapend" >> report.out
		echo -e "\nNon-polymorphic gap:\n$gap bp" >> report.out
	fi

fi

if [ $flag3 = 0 ]; then

	echo -e "Allele\\n\\t${name[1]}\\n" >> report.out
	echo -e "Allele\\n\\t${name[2]}\\n" >> report.out
	echo -e "Phase\\n\\t${name[3]}\\n" >> report.out

	##############################################################
	expr substr "${string[1]}" 1 ${gapstart} > temp
	read flagment1 < temp

	size2=`expr 271 - ${gapend}`
	expr substr "${string[1]}" ${gapend} ${size2} > temp
	read flagment2 < temp

        expr substr "${string[2]}" 1 ${gapstart} > temp
        read flagment3 < temp

        expr substr "${string[2]}" ${gapend} ${size2} > temp
        read flagment4 < temp

	gapsize=`expr ${gapend} - ${gapstart} - 1`
	gapstarta=`expr ${gapstart} + 1`
	expr substr "${string[1]}" ${gapstarta} ${gapsize} > temp
	read gapflagment < temp

	allele1=`echo ${flagment1}${gapflagment}${flagment2}`
	allele2=`echo ${flagment1}${gapflagment}${flagment4}`
	allele3=`echo ${flagment3}${gapflagment}${flagment2}`
	allele4=`echo ${flagment3}${gapflagment}${flagment4}`

        grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta > a.out
        grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
        grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out

	echo -e "Alleles detected" >> report.out
	grep "\*" a.out | cut -d '>' -f 2 >> report.out
	##############################################################

	if [[ $gap -ne 0 ]]; then
		echo -e "\nGap location:\n$gapstart\t$gapend" >> report.out
		echo -e "\nNon-polymorphic gap:\n$gap bp" >> report.out
	fi

fi

fi
