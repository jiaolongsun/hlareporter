rm -rf *_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/*_gen.fasta
cat A_gen.fasta B_gen.fasta C_gen.fasta F_gen.fasta G_gen.fasta H_gen.fasta DP*_gen.fasta DQ*_gen.fasta DR*_gen.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA-I_II_GEN.fasta
../bin/formatdb -p F -i HLA-I_II_GEN.fasta
bwa index -a is HLA-I_II_GEN.fasta
