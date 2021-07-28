#!/bin/bash
#rnaseqalign.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-b  --BAM File"
  echo "-g  --GTF File"
  echo "-s  --stranded"
  echo "-p  --Prefix for output file name"
  echo "Example: bash geneabudance.sh genect_rnaseq/geneabundance.sh -s stranded -g gtf_file -p pair_id -b bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :b:g:p:r:f:h opt
do
    case $opt in
        b) sbam=$OPTARG;;
        g) gtf=$OPTARG;;
        p) pair_id=$OPTARG;;
        r) index_path=$OPTARG;;
        s) stranded=$OPTARG;;
	f) filter=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]
then
    usage
fi
if [[ -z $stranded ]]
then
    stranded=0
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ $NPROC > 64 ]]
then
    NPROC=64
fi
if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH
    module load samtools/gcc/1.10 bedtools/2.29.0 
fi
baseDir="`dirname \"$0\"`"

samtools index -@ $NPROC ${sbam}
regtools junctions extract -o ${pair_id}.junc -s ${stranded} ${sbam}
regtools junctions annotate -o ${pair_id}.juncannot.txt ${pair_id}.junc ${index_path}/genome.fa ${index_path}/gencode.gtf
awk ' $14 == 1' ${pair_id}.juncannot.txt | awk '{print $1"\t"$2"\t"$3"\t"$6":"$5":"$14":"$9}' | grep -v chrom > ${pair_id}.knownjuct.bed
awk ' $14 == 0 && $5 > 9' ${pair_id}.juncannot.txt | awk '{print $1"\t"$2"\t"$3"\t"$6":"$5":"$14":"$9}' | grep -v chrom > ${pair_id}.exonskip.bed
bedtools intersect -wao -a ${pair_id}.exonskip.bed -b ${pair_id}.knownjuct.bed > ${pair_id}.skip2known.txt
bedtools intersect -wao -a ${pair_id}.exonskip.bed -b ${index_path}/gencode.exon.bed > ${pair_id}.exonskip.genes.txt
perl ${baseDir}/filter_exonskipping.pl -p $pair_id -r ${index_path} -s ${pair_id}.exonskip.genes.txt -k  ${pair_id}.skip2known.txt
