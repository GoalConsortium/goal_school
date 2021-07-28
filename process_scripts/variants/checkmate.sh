#!/bin/bash
#checkmate.sh

usage() {
  echo "-h Help documentation for checkmate"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-c  --Path to Reference GenomeNGS Checkmate SNP BED file"
  echo "Example: bash checkmate.sh -p prefix -r /path/GRCh38 -c /path/NGSCheckMate.bed -f"
  exit 1
}

OPTIND=1 # Reset OPTIND
while getopts :r:l:n:c:b:p:fh opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        b) sbam=$OPTARG;;
        n) normal=$OPTARG;;
	f) filter=1;;
	c) capture=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh	
    module load samtools/gcc/1.8 bcftools/gcc/1.8
    ncm=/project/shared/bicf_workflow_ref/seqprg/NGSCheckMate/ncm.py
    export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH
fi

if [[ -f /usr/local/bin/ncm.py ]]
then
    ncm=/usr/local/bin/ncm.py
fi

if [[ -z $capture ]]
then
    capture="${index_path}/NGSCheckMate.bed"
fi
if [[ -f "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
fi

for i in *.bam; do
    prefix="${i%.bam}"
    echo ${prefix}
    bcftools mpileup -A -d 1000000 -C50 -Ou --gvcf 0 -f ${reffa} -T ${capture} $i | bcftools call -m --gvcf 0 -Ov | bcftools convert --gvcf2vcf -f ${reffa} -Ov -o ${prefix}.vcf
done

python $ncm -V -d ./ -bed $capture -O ./ -N ${pair_id}
perl $baseDir/sequenceqc_somatic.pl -i ${pair_id}_all.txt -o ${pair_id}.sequence.stats.txt
