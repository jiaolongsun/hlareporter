#!/bin/bash

chmod 777 HLAreporter.sh

chmod 777 bin/HLAminer_ssake.pl
chmod 777 bin/HLAminer.pl
chmod 777 bin/HLAminer_classII.pl
chmod 777 bin/HLAminer_classII_flanking.pl
chmod 777 bin/blastall
chmod 777 bin/formatdb
chmod 777 bin/parseXMLblast.pl
chmod 777 bin/TASR

chmod 777 ssake_v3-8-tar/ssake_v3-8/SSAKE

chmod 777 Hydra-Version-0.5.3/bin/bamToFastq
chmod 777 bam2fastq-1.1.0/bam2fastq

chmod 777 mytest/HPTASRwgs_classI-II_main.sh
chmod 777 mytest/HPTASRwgs_classI-II_flanking.sh
chmod 777 mytest/HPTASRwgs_classI-II_addition.sh
chmod 777 mytest/HPTASRwgs_classI-II_mainII.sh
chmod 777 mytest/HPTASRwgs_classI-II_flankingII.sh
chmod 777 mytest/HPTASRwgs_classI-II_additionII.sh
chmod 777 mytest/HPTASRwgs_classI-II_ssake.sh
chmod 777 mytest/HPTASRwgs_classI-II_flanking_recursive.sh
chmod 777 mytest/HPTASRwgs_classI-II_flanking_recursive1.sh
chmod 777 mytest/*.sh

chmod 777 mytest/four_digit_prediction/4digit_map_HLA.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classII_main.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classII_flanking.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classII_addition.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classI_ssake.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classI_flanking.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classI_addition.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classI_flanking_recursive.sh
chmod 777 mytest/four_digit_prediction/4digit_predict_classI_flanking_recursive1.sh
chmod 777 mytest/four_digit_prediction/format_classI_flanking_bed.sh
chmod 777 mytest/four_digit_prediction/format_classI_flanking_bed_recursive.sh
chmod 777 mytest/four_digit_prediction/format_classI_flanking_bed_recursive1.sh
chmod 777 mytest/four_digit_prediction/format_classI_addition_bed.sh
chmod 777 mytest/four_digit_prediction/format_classII_flanking_bed.sh
chmod 777 mytest/four_digit_prediction/format_classII_addition_bed.sh
chmod 777 mytest/four_digit_prediction/phase_analyser.sh
chmod 777 mytest/four_digit_prediction/phase_miner.sh
chmod 777 mytest/four_digit_prediction/multigap_miner.sh
chmod 777 mytest/four_digit_prediction/recursive_phase_miner.sh
chmod 777 mytest/four_digit_prediction/recursive_multigap_miner.sh
chmod 777 mytest/four_digit_prediction/purify_efgh.sh
chmod 777 mytest/four_digit_prediction/purify_abc.sh
chmod 777 mytest/four_digit_prediction/purify_efgh_flanking.sh
chmod 777 mytest/four_digit_prediction/purify_abc_flanking.sh
chmod 777 mytest/four_digit_prediction/purify_drb345.sh
chmod 777 mytest/four_digit_prediction/purify_drdqdp.sh
chmod 777 mytest/four_digit_prediction/report_hla.sh
chmod 777 mytest/four_digit_prediction/report_ggroup_hla.sh
chmod 777 mytest/four_digit_prediction/establish_ssake_db.sh
chmod 777 mytest/four_digit_prediction/*.sh

perlloc=`whereis perl | cut -d ' ' -f 2`
sed -i "s:\/usr\/bin\/perl:`echo $perlloc`:g" bin/HLAminer_ssake.pl
sed -i "s:\/usr\/bin\/perl:`echo $perlloc`:g" bin/HLAminer.pl
sed -i "s:\/usr\/bin\/perl:`echo $perlloc`:g" bin/HLAminer_classII.pl
sed -i "s:\/usr\/bin\/perl:`echo $perlloc`:g" bin/HLAminer_classII_flanking.pl

echo "Program initialized!"

