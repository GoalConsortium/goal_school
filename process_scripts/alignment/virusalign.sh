#!/bin/bash
#abra.sh

usage(){
  echo "-h Help documentation for gatk4runner.sh"
  echo "-b --BAM File"
  echo "-r --Reference path e.g. GRCh38"
  echo "-p --Sample/Project ID"
}

OPTIND=1 #Reset OPTIN
while getopts :b:r:p:fh opt
do
  case $opt in
    b) bam=$OPTARG;;
    r) ref=$OPTARG;;
    p) pairid=$OPTARG;;
    f) filter=1;;
    h) usage;;
  esac
done

shift $(($OPTIND -1))

#Check for mandatory options 
if [[ -z $bam ]] ||  [[ -z $pairid ]] 
then  
  usage
fi

if [[ -z $SLURM_CPUS_ON_NODE ]]
then 
  SLURM_CPUS_ON_NODE=1
fi

reffa=${ref}/idt_virus_reference.fa
if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load bwa/intel/0.7.17 picard/2.10.3 samtools/1.6 
fi
baseDir="`dirname \"$0\"`"

samtools view -@ 8 -b -u -F 2 ${bam} |samtools sort -n - >unmapped.bam
java -Djava.io.tmpdir=./ -Xmx4g -jar $PICARD/picard.jar SamToFastq I=unmapped.bam FASTQ=unmapped.R1.fastq SECOND_END_FASTQ=unmapped.R2.fastq UNPAIRED_FASTQ=unmapped.unpaired.fastq

bwa mem -M -t 4 -R "@RG\tID:${pairid}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pairid}" ${reffa} unmapped.R1.fastq unmapped.R2.fastq > out.sam

samtools view -h -F 256 -b out.sam -o out.bam
samtools sort out.bam -o ${pairid}.viral.bam
samtools index ${pairid}.viral.bam
samtools idxstats ${pairid}.viral.bam >${pairid}.viral.idxstats.txt
if [[ $filter == 1 ]]
then
    perl $baseDir/filter_viral_idxstats.pl -p ${pairid} ${pairid}.viral.idxstats.txt
fi
