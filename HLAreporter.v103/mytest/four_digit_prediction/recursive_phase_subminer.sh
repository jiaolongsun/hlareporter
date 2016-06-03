#!/bin/bash

#extract name
rm temp_seq.dat #reset temp files
grep "\*" report.out > temprecord

#extract sequence
exec 3<> temprecord
while read str <&3
do

	fixedname=`echo ${str} | sed 's/\*/\\\*/g'`
	grep -A1 "$fixedname$" ../../database/HLA-I_II_GEN.fasta >> temp_seq.dat

done
exec 3>&-
#above close the info discription space 3

#start analyzing
num=0
exec 3<> temp_seq.dat
while read str <&3
do

	length=`expr length $str`

	if [ $length -lt 20 ]; then

		name[$num]=${str}
		
	else
		
		string[$num]=${str}
		num=`expr ${num} + 1`

	fi

done
exec 3>&-
#above close the info discription space 3

i=0
while [ "$i" -lt "546" ]; do
		j=0
		NTa[$i]=0
		NTc[$i]=0
		NTg[$i]=0
		NTt[$i]=0
	while [ "$j" -lt "$num" ]; do
		
		#get one nucleotide at a time
		substr=${string[$j]:${i}:1}
		
		if [[ ${substr} = "A" ]]; then
			NTa[$i]=1
		fi
		if [[ ${substr} = "C" ]]; then
			NTc[$i]=1
		fi
		if [[ ${substr} = "G" ]]; then
			NTg[$i]=1
		fi
		if [[ ${substr} = "T" ]]; then
			NTt[$i]=1
		fi

		j=`expr $j + 1`
	done
	
	NT[$i]=`expr ${NTa[$i]} + ${NTc[$i]} + ${NTg[$i]} + ${NTt[$i]}`

	if [[ ${NT[$i]} -gt "2" ]]; then
		exit
	fi

	i=`expr $i + 1`
	
done

#check each allele pair
i=0
while [ "$i" -lt "$num" ]; do
	j=`expr $i + 1`

	while [ "$j" -lt "$num" ]; do
		
		k=0
		while [ "$k" -lt "546" ]; do
			#get one nucleotide at a time
			substri=${string[$i]:${k}:1}	
			substrj=${string[$j]:${k}:1}	
			
			if [[ ${substri} = "${substrj}" ]]; then
				nt[$k]=1
			else
				nt[$k]=2
			fi

			#if there is a location not the same, check next
			if [[ ${nt[$k]} != ${NT[$k]} ]]; then
				break
			fi
		
			k=`expr $k + 1`	
			#echo -e "i is $i j is $j k is $k"
			#if all locations are the same, output and exit
			if [[ $k -eq "546" ]]; then
				echo -e "\nAllele pair" >> report.out
				echo -e "${name[$i]}" >> report.out
				echo -e "${name[$j]}" >> report.out
				exit
			fi

		done
		
		j=`expr $j + 1`
	done

	i=`expr $i + 1`
	
done

