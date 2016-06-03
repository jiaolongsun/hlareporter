###Run TASR
echo "Running TASR..."
../bin/TASR -f patient.fof -m 20 -s ../database/HLA-I_II_GEN.fasta -i 1 -b TASRhla
#../bin/TASR -f patient.fof -m 20 -i 1 -b TASRhla
###Restrict 40nt+ contigs

#first calculate the threshold coverage for noise detection
mymax=`sed -n 1p TASRhla.contigs | awk -F \cov '{print $2}'`
mymax=`echo "$mymax / 5"|bc`
echo threshold$mymax > $mymax

#then extract candidate contigs based on threshold coverage
cat $mymax TASRhla.contigs |perl -ne 'if(/threshold(\d+)/) {$myMax=$1;} if(/cov(\d+)/){if($1>=$myMax){$flag=1;print;}else{$flag=0;}}else{print if($flag);}' > TASRhla40.contigs
rm $mymax

###Create a [NCBI] blastable database
echo "Formatting blastable database..."
../bin/formatdb -p F -i TASRhla40.contigs
###Align contigs against database
echo "Aligning TASR contigs to HLA references..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -d ../database/HLA-I_II_GEN_FLANKING.fasta -i TASRhla40.contigs -o 0 > tig_vs_hla-ncbi.coord
###Align HLA references to contigs
echo "Aligning HLA references to TASR contigs..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -i ../database/HLA-I_II_GEN_FLANKING.fasta -d TASRhla40.contigs -o 0 > hla_vs_tig-ncbi.coord
###Predict HLA alleles
echo "Predicting HLA alleles..."
../bin/HLAminer_rnaseq_classII.pl -z 40 -i 100 -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla40.contigs -h ../database/HLA-I_II_GEN_FLANKING.fasta
