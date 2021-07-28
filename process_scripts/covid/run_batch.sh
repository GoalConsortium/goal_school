#!/bin/bash

index_path=/project/shared/bicf_workflow_ref/virus/covid

while read i; do
    line=($i)
    caseID=${line[0]}
    fq1=${line[1]}
    fq2=${line[2]}
    cd /project/bioinformatics/Cantarel_lab/shared/covid/jsorelle/batch1/$caseID
    sbatch -p 32GB,super /project/bioinformatics/Cantarel_lab/shared/covid/analyze_sample.sh -a $fq1 -b $fq2 -n $caseID -r $index_path
done <design.txt
