#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-v  --VCF File" 
  echo "Example: bash union.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:v:sh opt
do
    case $opt in
        p) pair_id=$OPTARG;;
	v) vcf=$OPTARG;;
	r) index_path=$OPTARG;;
	s) skipnorm=1;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load bedtools/2.26.0 samtools/gcc/1.8 bcftools/gcc/1.8 snpeff/4.3q 
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi

perl $baseDir\/uniform_vcf_gt.pl $pair_id $vcf
bgzip -f ${pair_id}.uniform.vcf
j=${pair_id}.uniform.vcf.gz
tabix -f $j
if [[ skipnorm==1 ]]
then
    cp $j ${pair_id}.norm.vcf.gz
else
    bcftools norm --fasta-ref $reffa -m - -Oz $j -o ${pair_id}.norm.vcf.gz
fi
