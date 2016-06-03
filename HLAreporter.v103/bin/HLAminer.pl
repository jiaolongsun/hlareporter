#!/usr/bin/perl

#AUTHOR
#   Rene Warren (c) 2006-2012
#   rwarren at bcgsc.ca

#NAME
#   HLAminer - Derivation of HLA class I and class II predictions from shotgun sequence datasets 

#SYNOPSIS
#   HLA class I and II Prediction by Targeted Assembly of Short Sequence Reads or Short Read Sequence Alignment

#DOCUMENTATION
#   HLAminer.readme

#LICENSE
#   HLAminer Copyright  (c) 2011-2012 Canada's Michael Smith Genome Science Centre.  All rights reserved.
#   TASR Copyright   (c) 2010-2011 Canada's Michael Smith Genome Science Centre.  All rights reserved.
#   SSAKE Copyright  (c) 2006-2011 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 2
#   of the License, or (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

use strict;
use Data::Dumper;
require "getopts.pl";
use vars qw($opt_b $opt_r $opt_c $opt_h $opt_z $opt_i $opt_q $opt_s $opt_l $opt_n $opt_a $opt_p);
&Getopts('b:r:c:h:z:i:p:q:s:l:n:a:');
my $version = "[v1.0.5]\nDerivation of HLA class I and II predictions from Shotgun sequence datasets";

###CHANGE IF NEEDED - DEFAULTS
my($sam_file,$blast_file,$reci_file,$tig_file,$hla_file,$min_tig_size, $min_seq_id, $min_phred, $min_score, $label, $nullalleles,$p_file) = ("NA","NA","NA","NA","NA",200,99,30,1000,"",0,"../database/hla_nom_p.txt");
my @hla1=('A','B','C','E','F','G');
my @hla2=('DPA1','DPB1','DQA1','DQB1','DRA','DRB1','DRB2','DRB3','DRB4','DRB5','DRB6','DRB7','DRB8','DRB9');
my @hlagenes = (@hla1,@hla2);


if(! $opt_h){
   print "Usage: $0 $version\n";
   print "--------------------------------------------------------------------\n";
   print "HPTASR (HLA Predictions by Targeted Assembly of Shotgun Reads):\n";
   print "-b blastn alignments.........................<tig_vs_hla-ncbi.coord>\n";
   print "-r reciprocal blastn.........................<hla_vs_tig-ncbi.coord>\n";
   print "-c contig fasta file.........................<TASRhla200.contigs>\n";
   print "-z minimum contig size.......................<200>\n";
   print "------------------------------- OR ---------------------------------\n";
   print "HPRA (HLA Predictions by Read Alignment):\n";
   print "-a sam alignments............................<ngs_vs_hla.sam>\n";
   print "--------------------------------------------------------------------\n";
   print "-h hla fasta file............................<HLA_ABC_CDS.fasta>\n";
   print "-p P-designation file........................<hla_nom_p.txt>\n";
   print "-i minimum \% sequence identity...............<99>\n";
   print "-q minimum log10 (phred-like) expect value...<30>\n";
   print "-s minimum score.............................<1000>\n";
   print "-n consider null alleles (1=yes/0=no)........<0>\n";
   die "-l label (run name) -optional-\n";
}

$blast_file = $opt_b if($opt_b);
$reci_file = $opt_r if($opt_r);
$tig_file = $opt_c if($opt_c);
$hla_file = $opt_h if($opt_h);
$sam_file = $opt_a if($opt_a);
$p_file = $opt_p if($opt_p);
$min_tig_size = $opt_z if($opt_z);
$min_seq_id = $opt_i if($opt_i);
$min_phred = $opt_q if($opt_q);
$min_score = $opt_s if($opt_s);
$nullalleles = $opt_n if($opt_n);
$label = $opt_l;
my $column = 13;
my $parameters = "$0 $version\n-b $blast_file\n-r $reci_file\n-c $tig_file\n-a $sam_file\n-h $hla_file\n-p $p_file\n-z $min_tig_size\n-i $min_seq_id\n-q $min_phred\n-s $min_score\n-n $nullalleles";

print "\nRunning $parameters\n";


if(! -f $p_file){
   my $message = "P-designation file $p_file doesn't not exist - will not use the P designation to summarize ambiguous HLA alleles.\n";
   print $message;
   print LOG $message;
}

if(! -f $blast_file && $opt_b ne ""){
   my $message = "File -b $blast_file doesn't not exist - FATAL.\n";
   print LOG $message;
   die $message;
}

if(! -f $hla_file && $opt_h ne ""){
   my $message = "File -h $hla_file doesn't not exist - FATAL.\n";
   print LOG $message;
   die $message;
}

if(! -f $reci_file && $opt_r ne ""){
   my $message = "File -r $reci_file doesn't not exist - FATAL.\n";
   print LOG $message;
   die $message;
}

if(! -f $tig_file && $opt_c ne ""){
   my $message = "File -c $tig_file doesn't not exist - FATAL.\n";
   print LOG $message;
   die $message;
}

if(! -f $sam_file && $opt_a ne ""){
   my $message = "File -a $sam_file doesn't not exist - FATAL.\n";
   print LOG $message;
   die $message;
}

my $pdesign = readPfile($p_file);
my $hlaprob = readFasta($hla_file,1);
my($tigprob,$top_hits,$reci_hits,$rec,$enum,$recI,$enumI);
my $mode = "";
my $labelmod ="";

if($label ne ""){
   $labelmod = "_" . $label;
}

if($opt_a ne ""){

   $mode = "HLAminer_HPRA" . $labelmod;
   $top_hits = readSam($sam_file,$min_tig_size,$column,0);
   ($rec,$enum) = mineAlign($top_hits,0,$min_seq_id);

}else{

   $mode = "HLAminer_HPTASR" . $labelmod;
   $tigprob = readFasta($tig_file,0);
   $hlaprob = readFasta($hla_file,1);

   $top_hits = readBlast($blast_file,$min_tig_size,$column,0);
   $reci_hits = readBlast($reci_file,$min_tig_size,$column,1);

   ($rec,$enum) = mineAlign($top_hits,0,$min_seq_id);
   ($recI,$enumI) = mineAlign($reci_hits,1,$min_seq_id);

}

my $log = $mode . ".log";
open(LOG,">$log") || die "Can't open $log for writing - Fatal.\n";
print "$mode log written to $log\n";

my $out = $mode . ".csv";
open(OUT,">$out") || die "Can't open $out for writing - Fatal.\n";

print LOG "$parameters\n";


&reportAlignments($rec,$enum,$recI,$tigprob,$hlaprob,$pdesign,$min_phred,$min_score,$nullalleles,$label);

print "$mode predictions written to $out\n";
close LOG;
close OUT;
exit;



#------------------------------------
sub readSam{

   my ($file,$min,$column,$reci) = @_;

   print "Reading read alignments file $file, please be patient...\n";

   open (IN, $file) || die "Can't open $file for reading -- fatal.\n";
   my $top_hits;

   my ($match,$counter,$qualct,$qualtot)=(0,0,0,0);

   while(<IN>){
      chomp;
      my @a = split(/\s+/);

      if($a[6] eq "=" && $a[7] && $a[10] ne "" && ($a[8]>0 || $a[8]<0) ){### means there's pairing
         my $template = $a[0];
         my $hit= $a[2];
         $counter++;

         if($counter==1 || $counter==2){

            ### sequence match
            if($_=~/MD\:Z\:(\S+)/){
               my @mat = split(/\D/,$1);
               foreach my $nummat(@mat){
                  $match+=$nummat;
               }
            }
            ### qualities
            my @b=split(//,$a[10]);
            foreach my $ch(@b){
               $qualct++;
               #$qualtot+= ord($ch) - $encoded;
            }
         }

         if($counter==2){

            my $tigsz = length($a[9])*2;
            my $seqid = $match / $tigsz *100;
            $template .= "size$tigsz|cov1.00";
            $top_hits->{$template}{$seqid}{$hit}="$seqid^^^^^^^^^^$seqid^";
            #-----------reset pair counter and other variables
            $counter=0;
            $match=0;
            $qualct=0;
            $qualtot=0;
         }###counter=2
      }##pairing
    }##

    close IN;

    return $top_hits;
}


#------------------------------------
sub readPfile{

   my $file = shift;
   my $pdesign;

   print "Reading p-nomenclature file $file...\n";

   open(IN,$file) || die "Can't open $file for reading -- fatal.\n";
   while(<IN>){
      chomp;
      #A*;23:03:01/23:03:02;23:03P
      my @elements=split(/\;/);
      my @alleles=split(/\//,$elements[1]);
      foreach my $allele(@alleles){
         my $hla = $elements[0] . $allele;
         if($elements[2] eq ""){
            $pdesign->{$hla} = $hla;
         }else{
            my $pd = $elements[0] . $elements[2];
            $pdesign->{$hla} = $pd;
         }
      }
   }
   close IN;
   return $pdesign;
}

#------------------------------------
sub mineAlign{
 
   my ($top_hits,$reciprocal,$minid) = @_;

   my $rec;
   my $enum;

   QRY:
   foreach my $qry (keys %$top_hits){
      my $si_list = $top_hits->{$qry};

      foreach my $si (sort {$b<=>$a} keys %$si_list){
         my $ct=0;
         my $hitlist = $si_list->{$si};
         foreach my $hit (keys %$hitlist){
             $ct++;
             my @d = split(/\^/,$hitlist->{$hit});
             #print "$qry,$hit\n" if ($reciprocal);
             if($reciprocal || (! $reciprocal && $d[10]>=$minid)){
                $rec->{$qry}{$hit} = "$d[10] %H/$d[0] %C";
                $enum->{$hit}{$qry} = "$d[10] %H/$d[0] %C";
             }
         }
         next QRY; ### allows top id to show
      }
   }
   return ($rec,$enum);
}

#------------------------------------
sub readFasta{

   my ($file,$type) = @_;
   my $spec={};
   my $sum={};
   my $prob={};
   my $entries = 0;

   print "Reading HLA reference sequence file $file...\n";

   open(IN,$file) || die "Can't open $file for reading -- fatal.\n";
   while(<IN>){
      chomp;
   
      if(/\>(\S+)/){
         $entries++;
         my $head = $1;
         #print "$head,$type\n";
         if($type){
            my $codvar = "NA";
            my $ag = "NA";
            my $supertype = "NA";
            ($codvar,$supertype,$ag) = ($1,$2,$3) if($head=~/^(((\w+)\*\d+)\:[^\:]*)/);
            $codvar=$head; ### added 13oct2011 to see effect of higher res.
            $spec->{$codvar}++;     #Counts number of sequence for coding variant
            $spec->{$supertype}++;
            $sum->{$ag}++;          #Counts total number of IMGT/HLA major genes sequences
            #print "$ag,$supertype,$codvar\n";
         }
      }
   }
   close IN;

   if($type){
      foreach my $cv(keys %$spec){
         my $ag = "NA";
         $ag = $1 if($cv =~ /^(\w+)\*/);
         $prob->{$cv} = $spec->{$cv} / $sum->{$ag};
      }
   }else{
      $prob->{'contig'} = 1/$entries;
   }

   return $prob;
}


#------------------------------------
sub reportAlignments{

   my ($rec,$enum,$proc,$tigprob,$hlaprob,$pdesign,$min_phred,$min_score,$nullalleles,$label) = @_;
 
   my $cnt;
   foreach my $tig(keys %$rec){### where tig is a contig
      my $tlist = $rec->{$tig};
      $cnt->{$tig}{'support'} = keys (%$tlist);
   }
   my $tnc;
   foreach my $hit(keys %$enum){
      my $tlist = $enum->{$hit};
      $tnc->{$hit}{'support'} = keys (%$tlist);
   }

   my $been;
   my $unknown;
   my $score;
   my $probability;
   my $probst;
   my $tracktig;
   #print Dumper($proc);
   my $globaltigprob;

   print LOG "\nAlignments to HLA\n-----------------\nlegend\n+:Coding variant seen by another contig previously\nRBH:Reciprocal Best Hit (if applicable)\n";
   print LOG "\nContig/Read pair,HLA [% HSP/% CONTIG/READ PAIR]\n";
   #print "Most likely HLA types from the supplied sequences (in order):\n";
   foreach my $tig (sort {$cnt->{$a}{'support'}<=>$cnt->{$b}{'support'}  } keys %$cnt){### CYCLE THROUGH CONTIGS
      my $hitlist=$rec->{$tig};
      my ($sizepoint,$covpoint)=(0,0);

      if($tig=~/size(\d+).*cov(\d+)/){  #contig21|size338|read163|cov48.71
         $sizepoint = $1;
         $covpoint = $2;
      }

      my $tmpscore = {};
      my $tmpscorest = {};
      my $seen = 0;
      foreach my $hit (sort {$tnc->{$b}{'support'}<=>$tnc->{$a}{'support'}} keys %$tnc){### CYCLE THROUGH HLA
         if(defined $hitlist->{$hit}){
            my $enumhit=$enum->{$hit};
            my $sn = keys(%$enumhit);
            #print "\t$hit [$hitlist->{$hit}] || Top hit in $sn total tigs:\n";
            my $indicator = " ";
            my $reci = " ";

            if($been->{$hit}){
               $indicator = "+";
               $seen = 1;
            }

            my $rbhpoint = 1;
            my $contigprob = 1;
            if(defined $proc->{$hit}{$tig}){####if true, indicates reciprocal best hit
               $reci = "RBH";
               $rbhpoint = 1;
               $contigprob = $tigprob->{'contig'};### the probability of hitting contig by chance
            }
            my $sidpoint = $1/100 if($hitlist->{$hit}=~/\/(\d+\.\d+) \%C/);
            my $ag = "NA";
            my $st = "NA";
            ($st,$ag) = ($1,$2) if($hit=~/^((\w+)\*\d+)/);
            my $notnullallele=1;
            if($hit =~ /N$/ && ! $nullalleles){##added 20dec2011 :: no Null alleles will be reported as final predictions
               $notnullallele=0;
            }

            if($notnullallele){##added 20dec2011
               $tmpscore->{$ag}{$hit}{'score'} = ($sizepoint * $covpoint * $sidpoint * $rbhpoint);
               $tmpscore->{$ag}{$hit}{'seen'} = $indicator;
               $globaltigprob->{$tig} += $hlaprob->{$hit}; 
               $tracktig->{$ag}{$hit}{$tig} = $contigprob;
            }

            ### PRINTS CONTIG-CENTRIC TABLE TO STDOUT 
            print LOG "$tig,$hit,$indicator [$hitlist->{$hit}] $reci\n";

            $been->{$hit}=1;
         }
      }
      ####compute score and probability (expect-E-value)
      foreach my $ag(keys %$tmpscore){
         my $hlalist = $tmpscore->{$ag};
         foreach my $codvar(keys %$hlalist){
            $score->{$ag}{$codvar} += $hlalist->{$codvar}{'score'};
         }
      }
      $seen = 0;
      ####

   }###End of Contigs
   
   ###Calculate P (E-value), once finished cycling through contigs

   foreach my $ag(keys %$tracktig){ 
      my $hlalist = $tracktig->{$ag};
      foreach my $codvar(keys %$hlalist){
         my $contiglist = $tracktig->{$ag}{$codvar};
         $probability->{$ag}{$codvar} = 1;
         foreach my $tigin (keys %$contiglist){
            my $localprob = $globaltigprob->{$tigin} * $tracktig->{$ag}{$codvar}{$tigin};### second variable factors in RBH
            $probability->{$ag}{$codvar} *= $localprob;### independent events 
         }
      }
   }

   ### TABLE SUMMARY
   my $bt;
   my $pdesignprob;
   my $trackgroupprob;
   print LOG "\n";
   print LOG "-" x 70, "\n";
   print LOG "Gene,HLA Supertype/Group,Protein coding allele,Score,Validated by N Contigs,x=contig seen,Chance (%),Expect (Eval) value,Confidence Score (-10 * log10(P))\n*Score = sum(Sz * cov * seq.id. * 2 [RBH])\n*P HLA = (P tig1 is HLA * P HLA is tig 1) * (P tig2 is HLA * P HLA is tig 2) * .. * (P tig(n) is HLA * P HLA is tig n) where P tig x is HLA = sum(P HLA for tig x)\n";
   print LOG "-" x 70, "\n";

   foreach my $ag(@hlagenes){
      print LOG "\nHLA-$ag\n";
      #my $list = $probability->{$ag};
      my $list = $score->{$ag};

      my $ct=0;
      my $predct=1;
      my $prevscore = "-10";
      my $prevsupertype = "NA";
      my $prevmess = "";

      foreach my $cv (sort {$list->{$b}<=>$list->{$a}} keys %$list){###sorted by decreasing score
         my $supertype = $1 if($cv=~/^([^\:]*)/);#### 14may2012
         my $numtig = $tracktig->{$ag}{$cv};
         my $tt = keys (%$numtig);
         my $mess= "";
         foreach my $tig(keys %$numtig){if(defined $bt->{$tig}){$mess.="x"; }  $bt->{$tig}++;}
         if($prevscore == $list->{$cv}){### means same contigs contribute to call (score)
            $mess = $prevmess;
         }         

         if($mess eq "" || length($mess) < $tt ){

            my $chance = (1-$probability->{$ag}{$cv}) * 100;
            $probability->{$ag}{$cv} = 1 if(! $probability->{$ag}{$cv});
            my $chancescore = (log($probability->{$ag}{$cv})/log(10)) * -10;
            printf LOG "\t$ag,$supertype,$cv,$score->{$ag}{$cv},$tt,$mess,%.4f,%.2e,%.1f\n", ($chance,$probability->{$ag}{$cv},$chancescore );

            ### LAST SUMMARY added on 13oct2011
            my $designate = $cv;
            $designate = $pdesign->{$cv} if(defined $pdesign->{$cv});
            my $supertypefromp = $1 if($designate=~/^([^\:]*)\:/);
            #### BY SCORE
            my $phred_like = (log($probability->{$ag}{$cv})/log(10)) * -10;
            if($phred_like >= $min_phred && $score->{$ag}{$cv} >= $min_score){### conditions for score/E val met

               $pdesignprob->{$ag}{$supertypefromp}{$score->{$ag}{$cv}}{$designate}{'p'}=$probability->{$ag}{$cv}; #e val
               $pdesignprob->{$ag}{$supertypefromp}{$score->{$ag}{$cv}}{$designate}{'s'}=$score->{$ag}{$cv};       #score
               $trackgroupprob->{$ag}{$supertypefromp} = $score->{$ag}{$cv} if($score->{$ag}{$cv} > $trackgroupprob->{$ag}{$supertypefromp});
            }
         }
         #### $mess="x" indicates that this allele has been characterized equally well by another, higher-scoring, allele 
         $prevscore = $list->{$cv};
         $prevsupertype = $supertype;
         $prevmess = $mess;
      }
   }

   reportHLA("I",\@hla1,$pdesignprob,$trackgroupprob);
   reportHLA("II",\@hla2,$pdesignprob,$trackgroupprob);

}

#------------------------------------
sub reportHLA{

   my ($type,$hlaarr,$pdesignprob,$trackgroupprob) = @_;
   my @hlaarr = @$hlaarr;

   print LOG "\n";
   print LOG "-" x 90, "\n";
   print LOG "SUMMARY $label\n"; 
   print LOG "MOST LIKELY HLA-$type ALLELES (Confidence (-10 * log10(Eval)) >= $min_phred Score >= $min_score)\n";
   print LOG "Gene,Prediction rank,Group,Allele,Score,Expect (Eval) value,Confidence (-10 * log10(Eval))\n";
   print LOG "-" x 90, "\n";
   print OUT "\n";
   print OUT "-" x 90, "\n";
   print OUT "SUMMARY $label\n";
   print OUT "MOST LIKELY HLA class $type ALLELES (Confidence (-10 * log10(Eval)) >= $min_phred Score >= $min_score)\n-Arranged by decreasing score, from most to less likely.\n-The prediction rank factors in the maximum score for each predicted allele.\n";
   print OUT "Allele,Score,Expect (Eval) value,Confidence (-10 * log10(Eval))\n";
   print OUT "-" x 90, "\n";

#   my $mem;#tracks HLA alleles reported, prevent multiple reporting (esp. for those with P designation)
   foreach my $ag(@hlaarr){#sort {$a cmp $b} keys %$pdesignprob){
      my $subgroup = $trackgroupprob->{$ag};
      print LOG "\nHLA-$ag\n";
      print OUT "\nHLA-$ag\n";
      my $ctpred=0;
      my $lastscore=0;
      foreach my $group (sort {$subgroup->{$b}<=>$subgroup->{$a}} keys %$subgroup){
         my $listscore=$pdesignprob->{$ag}{$group};
         my $flag=0;
         my $mem;#tracks HLA alleles reported, prevent multiple reporting (esp. for those with P designation)
         foreach my $score (sort {$b<=>$a} keys %$listscore){
            if(! $flag){
               if($score == $lastscore){
                  print LOG "\tPrediction #$ctpred - $group (same score as above)\n";
                  print OUT "\tPrediction #$ctpred - $group (same score as above)\n";
               }else{
                  $ctpred++;
                  print LOG "\tPrediction #$ctpred - $group\n";
                  print OUT "\tPrediction #$ctpred - $group\n";
               }
               $flag=1;
            }
            my $hlalist=$listscore->{$score};
            foreach my $hptasr(sort {$hlalist->{$b}{'s'} <=> $hlalist->{$a}{'s'}} keys %$hlalist){
               if(! defined $mem->{$hptasr}){
                  my $chance = (1-$hlalist->{$hptasr}{'p'}) * 100;
                  my $chancescore = (log($hlalist->{$hptasr}{'p'})/log(10)) * -10;
                  printf LOG "\t\t$ag,$ctpred,$group,$hptasr,$score,%.2e,%.1f\n", ($hlalist->{$hptasr}{'p'},$chancescore );
                  printf OUT "\t\t$hptasr,$score,%.2e,%.1f\n", ($hlalist->{$hptasr}{'p'},$chancescore );
                  $mem->{$hptasr}=1;
               }
            }
            #print OUT "\n";
            #print LOG "\n";
            $lastscore=$score;
         }   
         #print OUT "\n";
         #print LOG "\n";
      }
   }

   print LOG "\n";
}

#------------------------------------
sub readBlast{
   my ($file,$min,$column,$reci) = @_;

   print "Reading BLAST alignments file $file, please be patient...\n";
 
   open (IN, $file) || die "Can't open $file for reading -- fatal.\n";
   my $top_hits;

   while(<IN>){
      chomp;
      my @s = split(/\^/);
      my $tigsz=$1 if(/size(\d+)/);
      #print "no screen..\b";
      $s[2]=$1 if($s[2]=~/(\S+)/);
      $s[6]=$1 if($s[6]=~/(\S+)/);

      my ($x,$y) = ("","");
      my $z = 0;
      if($reci){
         $x = $s[2];###added 13oct2011
         $y = $s[6];
         $z = ($s[12] / $s[7]) * 100 if($s[7]>0); #% seq id over tig length
      }else{
         $y = $s[6];###added 13oct2011
         $x = $s[2];
         $z = $s[$column];  # if column is 0, then % seq id over tig length
      }

      $top_hits->{$x}{$z}{$y}=$_ if($tigsz >= $min);
   }

   close IN;
   return $top_hits;
}


