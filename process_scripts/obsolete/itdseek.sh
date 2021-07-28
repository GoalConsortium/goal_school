#!/bin/bash
#svcalling.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "Example: bash svcalling.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:l:p:fh opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;	
        p) pair_id=$OPTARG;;
	l) itdbed=$OPTARG;;
	g) snpeffgeno=$OPTARG;;
	f) filter=1;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"


# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi

source /etc/profile.d/modules.sh	
export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH

module load samtools/gcc/1.8 snpeff/4.3q vcftools/0.1.14 htslib/gcc/1.8 bcftools/gcc/1.8 bedtools/2.26.0 
stexe=`which samtools`

samtools view -@ $NPROC -L ${itdbed} ${sbam} | itdseek.pl --refseq ${reffa} --samtools ${stexe} --bam ${sbam} | vcf-sort | bedtools intersect -header -b ${itdbed} -a stdin | bgzip > ${pair_id}.itdseek.vcf.gz
tabix ${pair_id}.itdseek.vcf.gz
bcftools norm --fasta-ref $reffa -c w -m - -Ov ${pair_id}.itdseek.vcf.gz | java -Xmx30g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config $snpeffgeno - |bgzip > ${pair_id}.itdseek_tandemdup.vcf.gz

if [[ $filter == 1 ]]
then
    perl $baseDir/filter_itdseeker.pl -t ${pair_id} -d ${pair_id}.itdseek_tandemdup.vcf.gz
    mv ${pair_id}.itdseek_tandemdup.vcf.gz ${pair_id}.itdseek_tandemdup.unfilt.vcf.gz
    mv ${pair_id}.itdseek_tandemdup.pass.vcf ${pair_id}.itdseek_tandemdup.vcf
    bgzip ${pair_id}.itdseek_tandemdup.pass.vcf
fi
