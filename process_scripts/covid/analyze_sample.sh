#!/bin/bash

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-n  --NuCLIA CaseID"
  echo "-b  --targetbed"
  echo "-a  --archive"
  echo "Example: bash unify_case.sh -p prefix -r /path/GRCh38"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :n:a:b:r:h opt
do
    case $opt in
        n) caseID=$OPTARG;;
	a) fqr1=$OPTARG;;
 	b) fqr2=$OPTARG;;
        r) index_path=$OPTARG;;
	h) usage;;
    esac
done

shift $(($OPTIND -1))
#module load trimgalore/0.6.4 cutadapt/2.5 picard/2.10.3 samtools/gcc/1.10 bcftools/gcc/1.10.2 bedtools/2.29.2  snpeff/4.3q

module load spades/gcc/3.13.0 singularity/3.5.3 VarScan/2.4.2

if [[ -z $index_path ]]
then
    index_path=/project/shared/bicf_workflow_ref/virus/covid
fi

export PATH=$PATH:/project/shared/bicf_workflow_ref/seqprg/minimap2-2.16_x64-linux/

NPROC=`nproc`

#trim step
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/trim_galore.sif trim_galore -q 25 --illumina --gzip --length 35 --paired ${fqr1} ${fqr2}
#move files to new names
mv *_val_1.fq.gz ${caseID}.trim.R1.fastq.gz
mv *_val_2.fq.gz ${caseID}.trim.R2.fastq.gz

metaspades.py -1 ${caseID}.trim.R1.fastq.gz -2 ${caseID}.trim.R2.fastq.gz -o ${caseID}_assembly

minimap2 -a -o ${caseID}.sam -R "@RG\tID:${caseID}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${caseID}" --MD -t $NPROC -x sr ${index_path}/genome.mmi ${caseID}.trim.R1.fastq.gz ${caseID}.trim.R2.fastq.gz

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools view -1 -o output.unsort.bam ${caseID}.sam
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools sort -n --threads $NPROC -o ${caseID}.namesort.bam output.unsort.bam
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools view -h ${caseID}.namesort.bam > ${caseID}.namesort.sam

/project/shared/bicf_workflow_ref/seqprg/bin/primerclip /project/bioinformatics/Cantarel_lab/shared/covid/sarscov2_v2_masterfile.txt ${caseID}.namesort.sam ${caseID}.pclip.sam

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools view -1 -o output.unsort.bam ${caseID}.pclip.sam

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools sort --threads $NPROC -o ${caseID}.bam output.unsort.bam

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img java -XX:ParallelGCThreads=$NPROC -Djava.io.tmpdir=./ -Xmx16g  -jar /usr/local/bin/picard.jar MarkDuplicates I=${caseID}.bam O=${caseID}.dedup.bam M=${caseID}.dedup.stat.txt

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bedtools genomecov -bga -split -ibam ${caseID}.dedup.bam > ${caseID}.covhist.txt
awk '$4 < 5'  ${caseID}.covhist.txt |cut -f 1,2,3 > ${caseID}.lowdepth.bed

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img samtools mpileup -B -f ${index_path}/genome.fasta ${caseID}.dedup.bam > ${caseID}.pileup

VarScan mpileup2snp ${caseID}.pileup --strand-filter 0 --output-vcf 1 > ${caseID}.snp.vcf
VarScan mpileup2indel ${caseID}.pileup --strand-filter 0 --output-vcf 1 > ${caseID}.indel.vcf
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bgzip ${caseID}.snp.vcf
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bgzip ${caseID}.indel.vcf
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img tabix ${caseID}.snp.vcf.gz
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img tabix ${caseID}.indel.vcf.gz
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools concat ${caseID}.snp.vcf.gz ${caseID}.indel.vcf.gz | singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools sort -Oz -o ${caseID}.varscan.vcf.gz 
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img tabix ${caseID}.varscan.vcf.gz

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools mpileup --threads $NPROC -a 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/DP,FORMAT/SP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' -Ou -o ${caseID}.mpileup -A -d 1000000 -C50 -f ${index_path}/genome.fasta ${caseID}.bam
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools call -A --threads 10 -vmO z -o ${caseID}.vcf.gz ${caseID}.mpileup
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools view --threads $NPROC --min-af 0.5 -i 'FORMAT/DP[0] > 5' -o  ${caseID}.pass.vcf.gz -Oz ${caseID}.vcf.gz
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img tabix ${caseID}.pass.vcf.gz

singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img bcftools consensus -f ${index_path}/genome.fasta -m ${caseID}.lowdepth.bed -o ${caseID}.fasta ${caseID}.varscan.vcf.gz
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img java -Xmx10g -jar /usr/local/bin/snpEff/snpEff.jar -no-downstream -no-upstream -no-intergenic -lof -c /project/shared/bicf_workflow_ref/virus/snpeff/snpEff.config covid.wuhanHu1 ${caseID}.varscan.vcf.gz > ${caseID}.annot.vcf
singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/goalconsortium-variantcalling-1.0.7.img java -Xmx10g -jar /usr/local/bin/snpEff/SnpSift.jar extractFields ${caseID}.annot.vcf CHROM POS REF ALT ANN[*].GENE ANN[*].HGVS_C ANN[*].HGVS_P GEN[*].DP GEN[*].FREQ > ${caseID}.annot.txt


#singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/pangolin_v2.3.8.img pangolin jsorelle/batch1/batch1.fasta -t 2
