#!/bin/bash

module load bedtools/2.29.0
ln -s /project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref/genenames.txt .
perl /project/PHG/PHG_Clinical/devel/clinseq_workflows/process_scripts/genect_rnaseq/concat_cts.pl -o ./ */*/*.cts
perl /project/PHG/PHG_Clinical/devel/clinseq_workflows/process_scripts/genect_rnaseq/concat_fpkm.pl -o ./ */*/*.fpkm.txt
cut -f 2,4- countTable.fpkm.txt |perl -pi -e 's/SYMBOL/Hugo_Symbol/g' > fpkm.txt

ls ../*/CNV/*.txt | awk -F '/' '{print "cut -f 1-3,5",$0,"|bedtools intersect -wao -a stdin -b tempus.genes.hg19.bed | cut -f 1-3,4,8 >",$2".cnv_continuous.txt"}' |sh
ls ../*/CNV/*.txt | awk -F '/' '{print "cut -f 1-3,12",$0,"|bedtools intersect -wao -a stdin -b tempus.genes.hg19.bed | cut -f 1-3,4,8 >",$2".cnv_discreet.txt"}' |sh
 perl /project/PHG/PHG_Clinical/devel/clinseq_workflows/process_scripts/cbioportal/concat_cnvs.pl
