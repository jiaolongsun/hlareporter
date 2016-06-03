rm -rf *_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/A_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/B_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/C_gen.fasta
cat A_gen.fasta B_gen.fasta C_gen.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA_ABC_GEN.fasta
formatdb -p F -i HLA_ABC_GEN.fasta
bwa index -a is HLA_ABC_GEN.fasta
