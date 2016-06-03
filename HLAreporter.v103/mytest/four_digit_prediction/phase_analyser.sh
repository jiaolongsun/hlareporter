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

if [[ $linecount -lt 6 || $linecount -gt 20 ]]; then

	#echo ${linecount}
	#echo "The number of line is not expected"
	#exit

	#directly output all sequence info for phase combination analysis
	
	#reset
	rm a.out
	rm allele_pair.dat
	cat allele_sequences.dat > allele_pair.dat
	./phase_miner.sh $1 $2

elif [ $linecount -gt 6 ]; then
	
	#output their names to temp
	grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | grep -E ".*\*.*" > temprecord
	
	
	ROW=`expr $linecount / 2` 
	echo -e "row $ROW"
	i=1
	j=2
	
	rm a.out
	while [ "$i" -lt "$ROW" ]; do
		while [ "$j" -lt `expr $ROW + 1` ]; do
			echo "recursive phase analysis, please wait..."
			namei=`sed -n ${i}p temprecord`
			namej=`sed -n ${j}p temprecord`
						
			fixednamei=`echo ${namei} | sed 's/\*/\\\*/g'`
			fixednamej=`echo ${namej} | sed 's/\*/\\\*/g'`
			grep -A1 "$fixednamei$" ../../database/HLA-I_II_GEN.fasta >> allele_pair.dat
			grep -A1 "$fixednamej$" ../../database/HLA-I_II_GEN.fasta >> allele_pair.dat
			./recursive_phase_miner.sh $1 $2
			
			#cat allele_pair.dat
			rm allele_pair.dat #echo > allele_pair.dat
			j=`expr $j + 1`
		done
		i=`expr $i + 1`
		j=`expr $i + 1`
		
	done
	./recursive_phase_subminer.sh
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

#start comparing each nucleotide

location=1
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
		location=`expr ${location} + 1`

	elif [[ ${substr1} = ${substr3} ]]; then 

		flag2=1
		location=`expr ${location} + 1`

	elif [[ ${substr2} = ${substr3} ]]; then 

		flag1=1
		location=`expr ${location} + 1`

	else

		flag1=1
		flag2=1
		flag3=1
		location=`expr ${location} + 1`

	fi
}

done


#start printing the results

#reset
rm a.out
rm allele_pair.dat

#obtain 3 scores starting from row 2, then compare 3 scores
grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 5 > temprecord
score1=`sed -n 2p temprecord`
score2=`sed -n 3p temprecord`
score3=`sed -n 4p temprecord`

#if the three candidates are irrelevant to phase, output all and terminate
if [[ $flag1 = 1 && $flag2 = 1 && $flag3 = 1 || $score1 = $score2 || $score2 = $score3 || $score1 = $score3 ]]; then
	
	#output their names to temp
	grep -E '.*,.*,.*,.*,.*,.*,.*,.*,' ./$2/$1_exon23_high_resolution_multi_ref_flanking/HLAminer_HPTASR.log | sed 's/,/\t/g' | cut -f 4 | grep -E ".*\*.*" > temprecord
	
	
	ROW=`expr $linecount / 2` 
	echo -e "row $ROW"
	i=1
	j=2
	
	rm a.out
	while [ "$i" -lt "$ROW" ]; do
		while [ "$j" -lt `expr $ROW + 1` ]; do
			echo "recursive phase analysis, please wait..."
			namei=`sed -n ${i}p temprecord`
			namej=`sed -n ${j}p temprecord`
						
			fixednamei=`echo ${namei} | sed 's/\*/\\\*/g'`
			fixednamej=`echo ${namej} | sed 's/\*/\\\*/g'`
			grep -A1 "$fixednamei$" ../../database/HLA-I_II_GEN.fasta >> allele_pair.dat
			grep -A1 "$fixednamej$" ../../database/HLA-I_II_GEN.fasta >> allele_pair.dat
			./recursive_phase_miner.sh $1 $2
			
			#cat allele_pair.dat
			rm allele_pair.dat #echo > allele_pair.dat
			j=`expr $j + 1`
		done
		i=`expr $i + 1`
		j=`expr $i + 1`
		
	done

	#obtain contigs that support each allele
	ALLELE1=`echo ${name[1]} | cut -d '>' -f 2`
	FIXEDALLELE1=`echo ${ALLELE1} | sed 's/\*/\\\*/g'`
	CONFID1=`grep -E "cov.*${FIXEDALLELE1}," ./$2/$1_exon23_high_resolution_multi_ref/HLAminer_HPTASR.log | cut -d ',' -f 1 | sed ':label;N;s/\n/\ /;b label'`

	ALLELE2=`echo ${name[2]} | cut -d '>' -f 2`
	FIXEDALLELE2=`echo ${ALLELE2} | sed 's/\*/\\\*/g'`
	CONFID2=`grep -E "cov.*${FIXEDALLELE2}," ./$2/$1_exon23_high_resolution_multi_ref/HLAminer_HPTASR.log | cut -d ',' -f 1 | sed ':label;N;s/\n/\ /;b label'`

	ALLELE3=`echo ${name[3]} | cut -d '>' -f 2`
	FIXEDALLELE3=`echo ${ALLELE3} | sed 's/\*/\\\*/g'`
	CONFID3=`grep -E "cov.*${FIXEDALLELE3}," ./$2/$1_exon23_high_resolution_multi_ref/HLAminer_HPTASR.log | cut -d ',' -f 1 | sed ':label;N;s/\n/\ /;b label'`

	#print
	#echo -e "Alleles detected info" >> report.out
	#echo -e "${name[1]}\tSupported by ${CONFID1}" >> report.out
	#echo -e "${name[2]}\tSupported by ${CONFID2}" >> report.out
	#echo -e "${name[3]}\tSupported by ${CONFID3}" >> report.out

	#echo -e "\nAllele pair" >> report.out
	./recursive_phase_subminer.sh

	exit

fi

echo "Seperating:"

if [ $flag1 = 0 ]; then

	echo ${name[1]}
	#echo ${string[1]}

else

	echo ${name[1]} >> allele_pair.dat
	echo ${string[1]} >> allele_pair.dat	
	
fi

if [ $flag2 = 0 ]; then

	echo ${name[2]}
	#echo ${string[2]}

else

	echo ${name[2]} >> allele_pair.dat
	echo ${string[2]} >> allele_pair.dat		
fi

if [ $flag3 = 0 ]; then

	echo ${name[3]}	
	#echo ${string[3]}

else

	echo ${name[3]} >> allele_pair.dat
	echo ${string[3]} >> allele_pair.dat	
	
fi

./phase_miner.sh $1 $2

fi
