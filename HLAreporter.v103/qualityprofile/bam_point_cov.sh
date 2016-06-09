#!/bin/bash

if [ $# -ne 2 ]; then
	echo Usage: command ID gene
	echo Gene option: HLA_A HLA_B HLA_C HLA_DRB1 ...
	exit
fi

if [[ $2 != "HLA_A" && $2 != "HLA_B" && $2 != "HLA_C" && $2 != "HLA_DRB1" && $2 != "HLA_DRB3" && $2 != "HLA_DRB4" && $2 != "HLA_DRB5" && $2 != "HLA_DQB1" && $2 != "HLA_DQA1" && $2 != "HLA_DPB1" ]]; then
	echo "Gene entered out of scope"
	exit
fi

if [[ $2 = "HLA_A" || $2 = "HLA_B" || $2 = "HLA_C" ]]; then
	cp -p ../mytest/four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref/$1_*_exon23_high_resolution_multi_ref_mappedreads_target.fq ./
	bwa aln exon23_high_resolution_multi_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_1_$2_qc.sai
	bwa aln exon23_high_resolution_multi_ref.fa $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_2_$2_qc.sai
	bwa sampe exon23_high_resolution_multi_ref.fa $1_1_$2_qc.sai $1_2_$2_qc.sai $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_$2_qc.sam 
	samtools view -bS $1_$2_qc.sam > $1_$2_qc.bam
	samtools sort $1_$2_qc.bam -o $1_$2_qc_sorted.bam
	samtools index $1_$2_qc_sorted.bam
fi

if [[ $2 = "HLA_DRB1" || $2 = "HLA_DQB1" || $2 = "HLA_DRB3" || $2 = "HLA_DRB4" || $2 = "HLA_DRB5" || $2 = "HLA_DQA1" || $2 = "HLA_DPB1" ]]; then
	cp -p ../mytest/four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref/$1_*_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq ./
	cp -p ../mytest/four_digit_prediction/$2/$1_exon23_high_resolution_multi_ref/$1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.bam ./
	../Hydra-Version-0.5.3/bin/bamToFastq -bam $1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.bam -fq1 $1_1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.fq -fq2 $1_2_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.fq
	cat $1_1_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq > $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq
	cat $1_1_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.fq >> $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq
	cat $1_2_exon23_high_resolution_multi_ref_mappedreads_target_exon.fq > $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq
	cat $1_2_exon23_high_resolution_multi_ref_mappedreads_target_tempflanking.fq >> $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq
	bwa aln exon23_high_resolution_multi_ref.fa $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_1_$2_qc.sai
	bwa aln exon23_high_resolution_multi_ref.fa $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_2_$2_qc.sai
	bwa sampe exon23_high_resolution_multi_ref.fa $1_1_$2_qc.sai $1_2_$2_qc.sai $1_1_exon23_high_resolution_multi_ref_mappedreads_target.fq $1_2_exon23_high_resolution_multi_ref_mappedreads_target.fq > $1_$2_qc.sam 
	samtools view -bS $1_$2_qc.sam > $1_$2_qc.bam
	samtools sort $1_$2_qc.bam -o $1_$2_qc_sorted.bam
	samtools index $1_$2_qc_sorted.bam
fi

sed -i 's/\t.*\t.*/\t51\t52/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l > $1_$2_cov.out
sed -i 's/\t.*\t.*/\t52\t53/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t53\t54/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t54\t55/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t55\t56/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t56\t57/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t57\t58/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t58\t59/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t59\t60/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t60\t61/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t61\t62/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t62\t63/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t63\t64/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t64\t65/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t65\t66/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t66\t67/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t67\t68/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t68\t69/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t69\t70/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t70\t71/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t71\t72/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t72\t73/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t73\t74/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t74\t75/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t75\t76/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t76\t77/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t77\t78/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t78\t79/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t79\t80/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t80\t81/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t81\t82/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t82\t83/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t83\t84/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t84\t85/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t85\t86/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t86\t87/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t87\t88/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t88\t89/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t89\t90/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t90\t91/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t91\t92/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t92\t93/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t93\t94/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t94\t95/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t95\t96/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t96\t97/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t97\t98/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t98\t99/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t99\t100/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t100\t101/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t101\t102/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t102\t103/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t103\t104/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t104\t105/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t105\t106/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t106\t107/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t107\t108/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t108\t109/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t109\t110/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t110\t111/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t111\t112/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t112\t113/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t113\t114/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t114\t115/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t115\t116/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t116\t117/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t117\t118/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t118\t119/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t119\t120/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t120\t121/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t121\t122/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t122\t123/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t123\t124/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t124\t125/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t125\t126/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t126\t127/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t127\t128/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t128\t129/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t129\t130/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t130\t131/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t131\t132/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t132\t133/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t133\t134/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t134\t135/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t135\t136/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t136\t137/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t137\t138/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t138\t139/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t139\t140/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t140\t141/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t141\t142/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t142\t143/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t143\t144/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t144\t145/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t145\t146/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t146\t147/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t147\t148/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t148\t149/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t149\t150/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t150\t151/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t151\t152/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t152\t153/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t153\t154/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t154\t155/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t155\t156/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t156\t157/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t157\t158/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t158\t159/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t159\t160/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t160\t161/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t161\t162/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t162\t163/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t163\t164/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t164\t165/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t165\t166/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t166\t167/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t167\t168/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t168\t169/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t169\t170/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t170\t171/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t171\t172/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t172\t173/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t173\t174/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t174\t175/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t175\t176/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t176\t177/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t177\t178/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t178\t179/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t179\t180/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t180\t181/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t181\t182/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t182\t183/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t183\t184/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t184\t185/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t185\t186/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t186\t187/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t187\t188/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t188\t189/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t189\t190/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t190\t191/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t191\t192/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t192\t193/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t193\t194/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t194\t195/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t195\t196/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t196\t197/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t197\t198/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t198\t199/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t199\t200/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t200\t201/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t201\t202/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t202\t203/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t203\t204/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t204\t205/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t205\t206/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t206\t207/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t207\t208/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t208\t209/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t209\t210/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t210\t211/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t211\t212/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t212\t213/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t213\t214/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t214\t215/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t215\t216/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t216\t217/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t217\t218/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t218\t219/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t219\t220/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t220\t221/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t221\t222/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t222\t223/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t223\t224/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t224\t225/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t225\t226/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t226\t227/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t227\t228/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t228\t229/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t229\t230/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t230\t231/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t231\t232/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t232\t233/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t233\t234/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t234\t235/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t235\t236/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t236\t237/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t237\t238/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t238\t239/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t239\t240/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t240\t241/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t241\t242/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t242\t243/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t243\t244/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t244\t245/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t245\t246/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t246\t247/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t247\t248/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t248\t249/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t249\t250/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t250\t251/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t251\t252/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t252\t253/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t253\t254/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t254\t255/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t255\t256/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t256\t257/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t257\t258/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t258\t259/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t259\t260/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t260\t261/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t261\t262/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t262\t263/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t263\t264/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t264\t265/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t265\t266/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t266\t267/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t267\t268/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t268\t269/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t269\t270/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t270\t271/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t271\t272/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t272\t273/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t273\t274/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t274\t275/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t275\t276/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t276\t277/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t277\t278/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t278\t279/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t279\t280/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t280\t281/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t281\t282/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t282\t283/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t283\t284/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t284\t285/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t285\t286/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t286\t287/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t287\t288/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t288\t289/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t289\t290/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t290\t291/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t291\t292/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t292\t293/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t293\t294/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t294\t295/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t295\t296/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t296\t297/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t297\t298/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t298\t299/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t299\t300/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t300\t301/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t301\t302/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t302\t303/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t303\t304/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t304\t305/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t305\t306/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t306\t307/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t307\t308/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t308\t309/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t309\t310/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t310\t311/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t311\t312/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t312\t313/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t313\t314/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t314\t315/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t315\t316/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t316\t317/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t317\t318/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t318\t319/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t319\t320/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t320\t321/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out


if [[ $2 = "HLA_A" || $2 = "HLA_B" || $2 = "HLA_C" ]]; then

sed -i 's/\t.*\t.*/\t421\t422/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t422\t423/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t423\t424/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t424\t425/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t425\t426/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t426\t427/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t427\t428/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t428\t429/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t429\t430/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t430\t431/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t431\t432/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t432\t433/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t433\t434/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t434\t435/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t435\t436/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t436\t437/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t437\t438/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t438\t439/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t439\t440/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t440\t441/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t441\t442/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t442\t443/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t443\t444/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t444\t445/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t445\t446/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t446\t447/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t447\t448/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t448\t449/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t449\t450/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t450\t451/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t451\t452/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t452\t453/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t453\t454/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t454\t455/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t455\t456/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t456\t457/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t457\t458/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t458\t459/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t459\t460/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t460\t461/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t461\t462/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t462\t463/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t463\t464/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t464\t465/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t465\t466/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t466\t467/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t467\t468/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t468\t469/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t469\t470/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t470\t471/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t471\t472/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t472\t473/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t473\t474/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t474\t475/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t475\t476/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t476\t477/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t477\t478/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t478\t479/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t479\t480/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t480\t481/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t481\t482/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t482\t483/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t483\t484/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t484\t485/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t485\t486/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t486\t487/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t487\t488/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t488\t489/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t489\t490/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t490\t491/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t491\t492/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t492\t493/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t493\t494/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t494\t495/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t495\t496/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t496\t497/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t497\t498/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t498\t499/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t499\t500/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t500\t501/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t501\t502/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t502\t503/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t503\t504/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t504\t505/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t505\t506/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t506\t507/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t507\t508/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t508\t509/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t509\t510/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t510\t511/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t511\t512/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t512\t513/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t513\t514/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t514\t515/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t515\t516/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t516\t517/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t517\t518/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t518\t519/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t519\t520/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t520\t521/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t521\t522/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t522\t523/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t523\t524/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t524\t525/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t525\t526/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t526\t527/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t527\t528/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t528\t529/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t529\t530/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t530\t531/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t531\t532/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t532\t533/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t533\t534/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t534\t535/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t535\t536/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t536\t537/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t537\t538/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t538\t539/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t539\t540/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t540\t541/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t541\t542/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t542\t543/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t543\t544/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t544\t545/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t545\t546/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t546\t547/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t547\t548/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t548\t549/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t549\t550/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t550\t551/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t551\t552/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t552\t553/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t553\t554/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t554\t555/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t555\t556/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t556\t557/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t557\t558/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t558\t559/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t559\t560/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t560\t561/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t561\t562/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t562\t563/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t563\t564/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t564\t565/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t565\t566/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t566\t567/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t567\t568/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t568\t569/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t569\t570/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t570\t571/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t571\t572/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t572\t573/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t573\t574/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t574\t575/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t575\t576/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t576\t577/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t577\t578/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t578\t579/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t579\t580/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t580\t581/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t581\t582/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t582\t583/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t583\t584/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t584\t585/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t585\t586/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t586\t587/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t587\t588/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t588\t589/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t589\t590/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t590\t591/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t591\t592/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t592\t593/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t593\t594/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t594\t595/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t595\t596/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t596\t597/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t597\t598/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t598\t599/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t599\t600/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t600\t601/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t601\t602/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t602\t603/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t603\t604/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t604\t605/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t605\t606/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t606\t607/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t607\t608/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t608\t609/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t609\t610/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t610\t611/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t611\t612/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t612\t613/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t613\t614/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t614\t615/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t615\t616/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t616\t617/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t617\t618/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t618\t619/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t619\t620/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t620\t621/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t621\t622/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t622\t623/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t623\t624/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t624\t625/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t625\t626/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t626\t627/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t627\t628/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t628\t629/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t629\t630/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t630\t631/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t631\t632/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t632\t633/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t633\t634/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t634\t635/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t635\t636/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t636\t637/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t637\t638/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t638\t639/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t639\t640/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t640\t641/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t641\t642/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t642\t643/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t643\t644/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t644\t645/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t645\t646/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t646\t647/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t647\t648/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t648\t649/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t649\t650/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t650\t651/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t651\t652/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t652\t653/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t653\t654/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t654\t655/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t655\t656/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t656\t657/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t657\t658/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t658\t659/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t659\t660/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t660\t661/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t661\t662/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t662\t663/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t663\t664/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t664\t665/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t665\t666/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t666\t667/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t667\t668/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t668\t669/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t669\t670/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t670\t671/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t671\t672/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t672\t673/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t673\t674/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t674\t675/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t675\t676/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t676\t677/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t677\t678/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t678\t679/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t679\t680/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t680\t681/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t681\t682/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t682\t683/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t683\t684/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t684\t685/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t685\t686/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t686\t687/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t687\t688/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t688\t689/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t689\t690/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t690\t691/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t691\t692/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t692\t693/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t693\t694/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t694\t695/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t695\t696/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out
sed -i 's/\t.*\t.*/\t696\t697/g' $2.bed
samtools view -L $2.bed $1_$2_qc_sorted.bam | wc -l >> $1_$2_cov.out

fi


#summerize the percentage

DEPTH10=0
DEPTH20=0
DEPTH30=0
DEPTH50=0
TOTAL=0

exec 3<> $1_$2_cov.out
while read STRING <&3
do

	TOTAL=`expr $TOTAL + 1`

	if [ $STRING -ge 10 ]; then

		DEPTH10=`expr $DEPTH10 + 1`

	fi

	if [ $STRING -ge 20 ]; then

		DEPTH20=`expr $DEPTH20 + 1`

	fi

	if [ $STRING -ge 30 ]; then

		DEPTH30=`expr $DEPTH30 + 1`

	fi

	if [ $STRING -ge 50 ]; then


		DEPTH50=`expr $DEPTH50 + 1`

	fi

done
exec 3>&-
#above close the info discription space 3

#expr does not handle fraction part, so add two zero first
PERCENT10=`expr ${DEPTH10}00 / $TOTAL`
PERCENT20=`expr ${DEPTH20}00 / $TOTAL`
PERCENT30=`expr ${DEPTH30}00 / $TOTAL`
PERCENT50=`expr ${DEPTH50}00 / $TOTAL`

echo -e "10xcov%\t$PERCENT10\t20xcov%\t$PERCENT20\t30xcov%\t$PERCENT30\t50xcov%\t$PERCENT50" > $1_$2_covsummary.out

echo "HLA data quality profile"
cat $1_$2_covsummary.out

#relocate temp files
mkdir $2
mkdir $2/$1
mv $1* $2/$1
echo "Done generating HLA data quality profile"
echo "HLA report created in report.out"

