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
while getopts :r:p:g:v:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	    v) vcf=$OPTARG;;
	    g) snpeffgeno=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi

perl $baseDir\/uniform_vcf_gt.pl $pair_id $vcf
mv ${vcf} ${pair_id}.ori.vcf.gz
bgzip -f ${pair_id}.uniform.vcf
j=${pair_id}.uniform.vcf.gz
tabix -f $j
bcftools norm --fasta-ref $reffa -m - -Oz $j -o ${pair_id}.norm.vcf.gz
bash $baseDir/annotvcf.sh -p ${pair_id} -r $index_path -v ${pair_id}.norm.vcf.gz -g $snpeffgeno
vt decompose_blocksub ${pair_id}.annot.vcf.gz -p -a -o ${pair_id}.vcf
bgzip -f ${pair_id}.vcf
