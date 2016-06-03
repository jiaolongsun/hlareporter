###Run TASR
echo "Running TASR..."
../bin/TASR -f patient.fof -m 50 -o 5 -s ../database/HLA-I_II_GEN.fasta -i 1 -b TASRhla
#../bin/TASR -f patient.fof -m 20 -i 1 -b TASRhla
###Restrict 90nt+ contigs
cat TASRhla.contigs |perl -ne 'if(/size(\d+)/){if($1>=90){if(/cov(\d+)/){if($1>=5) {$flag=1;print;}else{$flag=0;}}}else{$flag=0;}}else{print if($flag);}' > TASRhla90.contigs
###Create a [NCBI] blastable database
echo "Formatting blastable database..."
../bin/formatdb -p F -i TASRhla90.contigs
###Align contigs against database
echo "Aligning TASR contigs to HLA references..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -d ../database/HLA-I_II_GEN.fasta -i TASRhla90.contigs -o 0 > tig_vs_hla-ncbi.coord
###Align HLA references to contigs
echo "Aligning HLA references to TASR contigs..."
../bin/parseXMLblast.pl -c ncbiBlastConfig.txt -i ../database/HLA-I_II_GEN.fasta -d TASRhla90.contigs -o 0 > hla_vs_tig-ncbi.coord
###Predict HLA alleles
echo "Predicting HLA alleles..."
../bin/HLAminer_classII.pl -z 90 -i 100 -b tig_vs_hla-ncbi.coord -r hla_vs_tig-ncbi.coord -c TASRhla90.contigs -h ../database/HLA-I_II_GEN.fasta
