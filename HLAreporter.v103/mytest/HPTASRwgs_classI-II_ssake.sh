file1=`sed -n 1p patient.fof`
file2=`sed -n 2p patient.fof`
#merge two fq
paste -d ":" $file1 $file2 > paired.fa
#add one more \n in the end so that the following replacement command can handle
echo -e "\n" >> paired.fa
#normalize format
sed -i ':label;N;s/\+.*\n.*\n//;b label' paired.fa
sed -i 's/\/1.*/:200/g' paired.fa
sed -i 's/@/>/g' paired.fa
#delete the last two \n after normalization such that ssake can handle
sed -i ':label;N;s/\n\n//;b label' paired.fa
#delete all illegal rows with N nt
sed -i ':label;N;s/\n/;/;b label' paired.fa
sed -i 's/;>/\n>/g' paired.fa
sed -i 's/.*;.*N.*//g' paired.fa
sed -i '/^$/d' paired.fa
sed -i 's/;/\n/g' paired.fa

###Run SSAKE
echo "Running SSAKE..."
../ssake_v3-8-tar/ssake_v3-8/SSAKE -f paired.fa -m 50 -o 2 -r 0.7 -p 1 -c 1 -e 0.75 -k 2 -w 2
mv paired.fa.ssake*.contigs TASRhla.contigs

###Restrict 90nt+ contigs
cat TASRhla.contigs |perl -ne 'if(/size(\d+)/){if($1>=90){if(/cov(\d+)/){if($1>=5) {$flag=1;print;}else{$flag=0;}}}else{$flag=0;}}else{print if($flag);}' > TASRhla90.contigs
###Create a [NCBI] blastable database
echo "Formatting blastable database..."
../bin/formatdb -p F -i TASRhla90.contigs
###Align contigs against database
echo "Aligning TASR contigs to HLA references..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -d ../database/HLA-I_II_GEN_SSAKE.fasta -i TASRhla90.contigs -o 0 > tig_vs_hla-ncbi.coord
###Align HLA references to contigs
echo "Aligning HLA references to TASR contigs..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -i ../database/HLA-I_II_GEN_SSAKE.fasta -d TASRhla90.contigs -o 0 > hla_vs_tig-ncbi.coord
###Predict HLA alleles
echo "Predicting HLA alleles..."
../bin/HLAminer_ssake.pl -z 90 -i 100 -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla90.contigs -h ../database/HLA-I_II_GEN_SSAKE.fasta
