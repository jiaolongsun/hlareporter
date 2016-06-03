#!/bin/bash

#check if the number of line is OK

linecount=`wc -l allele_pair.dat | sed 's/ /\t/' | cut -f 1`

if [[ $linecount -ne 4 ]]; then
	#echo ${linecount}
	#echo "The number of line is not expected"
	exit
fi

#calculate fastq read length

LENGTH=`cat $2/$1_exon23_high_resolution_multi_ref/$1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq | awk '{if(NR%4==2) {print length($1); exit 0;}}'`

if [ ! $LENGTH ]; then
	echo "NULL read length! program teminated!"
	exit
fi

#echo "Fastq read length: ${LENGTH}"

#start reading each sequence and save it

num=1
declare string[2]
declare name[2]

#cat input.txt | while read str; this pipeline will create a subshell, where any changes in it cannot be accessed by outside
#now get all the info from a file, and assign an id 3 to describe this info, then read things from this info

exec 3<> allele_pair.dat
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

rm gapinfo.dat >/dev/null 2>&1
#rm a.out
location=0
oldlocation=1
flag1=0
flag2=0
string1=${string[1]}
string2=${string[2]}

while true; do
{
	#get one nucleotide at a time
	substr1=${string[1]:${location}:1}
	substr2=${string[2]:${location}:1}

	if [[ ${substr1} = "" ]]; then
	
		break
	
	fi

	if [[ ${substr1} = ${substr2} ]]; then

		location=`expr ${location} + 1`

	else

		flag1=1
		gap=`expr ${location} - ${oldlocation}`
		#echo $gap
		if [[ $gap -ge $LENGTH && $oldlocation != 1 && $oldlocation != 271 ]]; then

			echo -e "`expr $oldlocation + 1` `expr $location + 1`" >> gapinfo.dat

		fi
		oldlocation=${location}
		location=`expr ${location} + 1`

	fi

	if [[ $location = 271 ]]; then
		
		#due to meaning of substr statement, real location of currently readed nt is 1 greater than $location
		oldlocation=271

	fi
}
done


#start printing the results

gapnumber=`wc -l gapinfo.dat | cut -d ' ' -f 1`
if [[ $gapnumber > 1 ]]; then
	./recursive_multigap_miner.sh ${string[1]} ${string[2]}
	exit
fi

gaplocation1=`sed -n 1p gapinfo.dat | cut -d ' ' -f 1`
gaplocation2=`sed -n 1p gapinfo.dat | cut -d ' ' -f 2`

gaplocation1a=`expr ${gaplocation1} + 1`
gaplocation2a=`expr ${gaplocation2} + 1`


if [[ ${gaplocation1} = "" ]]; then

	expr substr "${string[1]}" 1 270 > temp
        read flagment1 < temp

        expr substr "${string[1]}" 271 276 > temp
        read flagment2 < temp

        expr substr "${string[2]}" 1 270 > temp
        read flagment3 < temp

        expr substr "${string[2]}" 271 276 > temp
        read flagment4 < temp

        allele1=`echo ${flagment1}${flagment2}`
        allele2=`echo ${flagment1}${flagment4}`
        allele3=`echo ${flagment3}${flagment2}`
        allele4=`echo ${flagment3}${flagment4}`

	#touch a.out
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

        #grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta >> a.out
        #grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
        #grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
        #grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out

fi

if [[ ${gaplocation2} != "" &&  ${gaplocation2} -le 270 ]]; then

	expr substr "${string[1]}" 1 ${gaplocation1} > temp
	read flagment1 < temp

	size2=`expr 271 - ${gaplocation2}`
	expr substr "${string[1]}" ${gaplocation2} ${size2} > temp
	read flagment2 < temp

	size3=`expr 547 - 271`
	expr substr "${string[1]}" 271 ${size3} > temp
	read flagment3 < temp

        expr substr "${string[2]}" 1 ${gaplocation1} > temp
        read flagment4 < temp

        expr substr "${string[2]}" ${gaplocation2} ${size2} > temp
        read flagment5 < temp

        expr substr "${string[2]}" 271 ${size3} > temp
        read flagment6 < temp

	gapsize=`expr ${gaplocation2} - ${gaplocation1a}`
	expr substr "${string[1]}" ${gaplocation1a} ${gapsize} > temp
	read gapflagment < temp

	allele1=`echo ${flagment1}${gapflagment}${flagment2}${flagment3}`
	allele2=`echo ${flagment1}${gapflagment}${flagment2}${flagment6}`
	allele3=`echo ${flagment1}${gapflagment}${flagment5}${flagment3}`
	allele4=`echo ${flagment1}${gapflagment}${flagment5}${flagment6}`
	allele5=`echo ${flagment4}${gapflagment}${flagment2}${flagment3}`
	allele6=`echo ${flagment4}${gapflagment}${flagment2}${flagment6}`
	allele7=`echo ${flagment4}${gapflagment}${flagment5}${flagment3}`
        allele8=`echo ${flagment4}${gapflagment}${flagment5}${flagment6}`

	#touch a.out
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
	#grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele5$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele6$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele7$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele8$" ../../database/HLA-I_II_GEN.fasta >> a.out

fi

if [[ ${gaplocation1} != "" && ${gaplocation1} > 270 ]]; then

        expr substr "${string[1]}" 1 270 > temp
        read flagment1 < temp

	size2=`expr ${gaplocation1a} - 271`
        expr substr "${string[1]}" 271 ${size2} > temp
        read flagment2 < temp

        size3=`expr 547 - ${gaplocation2}`
        expr substr "${string[1]}" ${gaplocation2} ${size3} > temp
        read flagment3 < temp

        expr substr "${string[2]}" 1 270 > temp
        read flagment4 < temp

        expr substr "${string[2]}" 271 ${size2} > temp
        read flagment5 < temp

        expr substr "${string[2]}" ${gaplocation2} ${size3} > temp
        read flagment6 < temp	

        gapsize=`expr ${gaplocation2} - ${gaplocation1a}`
        expr substr "${string[1]}" ${gaplocation1a} ${gapsize} > temp
        read gapflagment < temp

        allele1=`echo ${flagment1}${flagment2}${gapflagment}${flagment3}`
        allele2=`echo ${flagment1}${flagment2}${gapflagment}${flagment6}`
        allele3=`echo ${flagment1}${flagment5}${gapflagment}${flagment3}`
        allele4=`echo ${flagment1}${flagment5}${gapflagment}${flagment6}`
        allele5=`echo ${flagment4}${flagment2}${gapflagment}${flagment3}`
        allele6=`echo ${flagment4}${flagment2}${gapflagment}${flagment6}`
        allele7=`echo ${flagment4}${flagment5}${gapflagment}${flagment3}`
        allele8=`echo ${flagment4}${flagment5}${gapflagment}${flagment6}`

	#touch a.out
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


	#grep -B1 "$allele1$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele2$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele3$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele4$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele5$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele6$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele7$" ../../database/HLA-I_II_GEN.fasta >> a.out
	#grep -B1 "$allele8$" ../../database/HLA-I_II_GEN.fasta >> a.out

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

