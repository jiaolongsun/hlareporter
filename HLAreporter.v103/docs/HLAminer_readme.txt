Manual Reference Pages – HLAminer - Derivation of HLA class I and class II predictions from shotgun sequence datasets
* This manual assumes that you have a working knowledge of unix, and some shell and perl scripting experience

NAME
  HLAminer - Derivation of HLA class I and class II predictions from shotgun sequence datasets

CONTENTS
  SYNOPSIS
  LICENSE
  OVERVIEW
  DESCRIPTION
  INSTALL
  COMMANDS AND OPTIONS
  DATABASES
  AUTHORS
  SEE ALSO

--------
SYNOPSIS
========

  For RNAseq:
  1. Copy ./test-demo/    eg. cp -rf test-demo foo
  2. In folder "foo", edit the patient.fof file to point to your NGS RNAseq data.  Ensure all paths are ok.
  3. For HLA Predictions by Targeted Assembly of Shotgun Reads: execute ./HLAminer/foo/HPTASRrnaseq.sh 
     For HLA Predictions by Read Alignment: execute ./HLAminer/foo/HPRArnaseq.sh

LICENSE
=======

  HLAminer Copyright  (c) 2012 Canada's Michael Smith Genome Science Centre.  All rights reserved.
  TASR Copyright  (c) 2010-2011 Canada's Michael Smith Genome Science Centre.  All rights reserved.
  SSAKE Copyright (c) 2006-2011 Canada's Michael Smith Genome Science Centre.  All rights reserved.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.


OVERVIEW
========

Derivation of HLA class I and class II predictions from shotgun sequence datasets (HLAminer) by:
1) Targeted Assembly of Shotgun Reads (HPTASR)
2) Read Alignment (HPRA)


DESCRIPTION
===========

The HLA prediction by targeted assembly of short sequence reads (HPTASR), performs targeted de novo assembly of HLA NGS reads and align them to reference HLA alleles from the IMGT/HLA sequence repository using commodity hardware with standard specifications (<2GB RAM, 2GHz).  Putative HLA types are inferred by mining and scoring the contig alignments and an expect value is determined for each.  The method is accurate, simple and fast to execute and, for transcriptome data, requires low depth of coverage. Known HLA class I/class II reference sequences available from the IMGT/HLA public repository are read by TASR using default options (Warren and Holt 2011) to create a hash table of all possible 15 nt words (kmers) from these reference sequences.   Subsequently, NGS data sets are interrogated for the presence of one of these kmers (on either strand) at the 5’ or 3’ start. Whenever an HLA word is identified, the read is recruited as a candidate for de novo assembly. Upon de novo assembly of all recruited reads, a set of contigs is generated.  Only sequence contigs equal or larger than 200nt in length are considered for further analysis, as longer contigs better resolve HLA allelic variants.  Reciprocal BLASTN alignments are performed between the contigs and all HLA allelic reference sequences. HPTASR mines the alignments, scoring each possible HLA allele identified, computing and reporting an expect value (E-value) based on the chance of contigs characterizing given HLA alleles and, reciprocally, the chance of reference HLA alleles aligning best to certain assembled contig sequences

The HLA prediction from direct read alignment (HPRA) method is conceptually simpler and faster to execute, since paired reads are aligned up-front to reference HLA alleles.  Alignments from the HPTASR and HPRA methods are processed by the same software (HLAminer.pl) to derive HLA-I predictions by scoring and evaluating the probability of each candidate bearing alignments.


INSTALL
=======

1. Download and decompress the tar ball
gunzip HLAminer.tar.gz
tar -xvf HLAminer.tar
2. Make sure you see the following directories:
./bin
./databases
./docs
./test-demo

3. Read the docs in the ./docs/ folder
4. Change/Add/Adjust the perl shebang line of each .pl and .sh script in the ./bin/ folder as needed

From direct Read Alignment (HPRA, faster but less accurate):
HPRArnaseq_classI.sh
HPRArnaseq_classI-II.sh
HPRAwgs_classI.sh
HPRAwgs_classI-II.sh

From Targeted Assembly (HPTASR, longer but more accurate):
HPTASRrnaseq_classI.sh
HPTASRrnaseq_classI-II.sh
HPTASRwgs_classI.sh
HPTASRwgs_classI-II.sh

*Running HPTASRwgs(rnaseq)_classI-II.sh will take longer than HPTASRwgs(rnaseq)_classI.sh, due to the reciprocal BLAST step.  You may remove this step from the former (and HLAminer.pl command) to speed things up.  However, this step is helpful in weeding out spurious alignments to HLA references.  That said, if you're solely interested in HLA-I, you have the option to run the latter set of scripts [HPTASRwgs(rnaseq)_classI.sh].

HLAminer.pl
parseXMLblast.pl
TASR

5. You must install perl module Bio::SearchIO to use HPTASR
6. Edit the fullpath location of bwa and other software dependencies in the shell scripts in the ./bin/ folder, as needed
7. For your convenience, ncbi blastall and formatdb have been placed in the ./bin/ folder and executed from the following shell scripts:

NAME,PROCESS,NGS DATA TYPE,PREDICTIONS
HPRArnaseq_classI.sh,Paired read alignment,RNAseq (transcriptome),HLA-I A,B,C genes
HPRArnaseq_classI-II.sh,Paired read alignment,RNAseq (transcriptome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPRAwgs_classI.sh,Paired read alignment,Exon capture (exome) and WGS (genome),HLA-I A,B,C genes
HPRAwgs_classI-II.sh,Paired read alignment,Exon capture (exome) and WGS (genome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPTASRrnaseq_classI.sh,Targeted assembly of sequence reads,RNAseq (transcriptome),HLA-I A,B,C genes
HPTASRrnaseq_classI-II.sh,Targeted assembly of sequence reads,RNAseq (transcriptome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

HPTASRwgs_classI.sh,Targeted assembly of sequence reads,Exon capture (exome) and WGS (genome),HLA-I A,B,C genes
HPTASRwgs_classI-II.sh,Targeted assembly of sequence reads,Exon capture (exome) and WGS (genome),HLA-I A,B,C and HLA-II DP,DQ,DR genes

Run those scripts by specifying the relative path ../bin/blastall or ../bin/formatdb in the shell scripts AND config file "ncbiBlastConfig.txt".  Make sure HPRA and HPTASR are running in same-level directories to ../bin/ (eg. ../test-demo/)

8. Before running on your data, inspect the ./test-demo/ folder and familiarize yourself with the files and the execution.  
9. When you are ready and the demo works well, place the fullpath location of your short read fastq or fasta files in the "patient.fof" file.
10. Make sure that the following files are in your working directory:

patient.fof
ncbiBlastConfig.txt


COMMANDS AND OPTIONS
====================

The shell scripts are set to filter out short (<200) contigs that would blur HPTAR predictions.  Feel free to adjust as you see fit.

Likewise, HLAminer.pl runs with the set defaults:
-z minimum contig size.......................<200> (HPTASR)
-i minimum % sequence identity...............<99>  (HPTASR / HPRA)
-q minimum log10 (phred-like) expect value...<30>  (HPTASR / HPRA)
-s minimum score.............................<1000> (HPTASR / HPRA)
-n consider null alleles (1=yes/0=no)........<0> (HPTASR / HPRA)

The minimum sequence identity applies to the short read paired alignment or blast alignment, depending on the choice made.  HLA predictions with a phred-like expect value lower than -q or a score lower than -s will not be diplayed.  Because IMGT/HLA reports numerous null alleles, an option exist to consider or not these unexpressed alleles. 


DATABASES
=========

Follow these instructions to download updated HLA sequences from ebi/imgt (shell scripts to automatically download and format the databases exist in ./database/) and refer to README.txt in the ./database directory:


1) Coding HLA sequences

HLA CDS sequences from:

wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/A_nuc.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/B_nuc.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/C_nuc.fasta
cat A_nuc.fasta B_nuc.fasta C_nuc.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA_ABC_CDS.fasta
../bin/formatdb -p F -i HLA_ABC_CDS.fasta
/home/pubseq/BioSw/bwa/bwa-0.5.9/bwa index -a is HLA_ABC_CDS.fasta

2) HLA genomic sequences 

To make the HLA genomic sequence database, execute these unix commands:

wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/A_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/B_gen.fasta
wget ftp://ftp.ebi.ac.uk/pub/databases/imgt/mhc/hla/C_gen.fasta
cat A_gen.fasta B_gen.fasta C_gen.fasta | perl -ne 'chomp;if(/\>\S+\s+(\S+)/){print ">$1\n";}else{print "$_\n";}' > HLA_ABC_GEN.fasta
../bin/formatdb -p F -i HLA_ABC_GEN.fasta
/home/pubseq/BioSw/bwa/bwa-0.5.9/bwa index -a is HLA_ABC_GEN.fasta


3) P designation files

Upgrade the P designation info from:

info:
http://hla.alleles.org/wmda/index.html

file:
http://hla.alleles.org/wmda/hla_nom_p.txt


OUTPUT FILES
============

HLA predictions from read pair alignments:

HLAminer_HPRA.log
HLAminer_HPRA.csv

HLA predictions from targeted assemblies:

HLAminer_HPTASR.log
HLAminer_HPTASR.csv

The .log file tracks the process of HLA mining. It contains the following information:
-HLAminer command and parameters utilized
-Contig/read pair alignment output and best HLA hit for each
-Initial gene summary, score and expect value
-Final summary, listing all predictions by highest score (more likely).

The .csv file contains HLAminer predictions.  Predictions are listed by HLA
gene and ranked by highest score.  Predictions 1) and 2) are expected to
represent the two alleles for each.

eg.
----------------------------------------------------------------------
SUMMARY
MOST LIKELY HLA-I ALLELES (Confidence (-10 * log10(Eval)) >= 30, Score >= 500)
Allele,Score,Expect (Eval) value,Confidence (-10 * log10(Eval))
----------------------------------------------------------------------

HLA-A
Prediction #1 - A*26
        A*26:33,4179,5.22e-124,1232.8

Prediction #2 - A*33
        A*33:24,1791,2.41e-75,746.2

Prediction #3 - A*68
        A*68:05,597,1.85e-10,97.3


From these predictions, the individual is expected to be heterozygous with HLA-I A alleles A*26 and
A*33. The chance, expect value and confidence are different representations of the same
metric. The Confidence represents the Expect value - Eval - as a score in a
manner analoguous to the phred score employed in sequencing to quickly assess
the likelihood of a base being correct. 

Predictions/read pair are ambiguous when there are multiple predicted allele groups and/or protein coding alleles with the same score.


AUTHORS
=======

Rene Warren
rwarren at bcgsc.ca


SEE ALSO
========

Warren RL, Choe G, Castellarin M, Munro S, Moore R, Holt RA. 2012 
HLAminer: Derivation of HLA types from shotgun sequence datasets (submitted)

