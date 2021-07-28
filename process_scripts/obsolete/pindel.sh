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
while getopts :r:l:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	l) idtbed=$OPTARG;;
	g) snpeffgeno=$OPTARG;;
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

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi

source /etc/profile.d/modules.sh	

genomefiledate=`find ${reffa} -maxdepth 0 -printf "%TY%Tm%Td\n"`
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi

module load samtools/gcc/1.8 bcftools/gcc/1.8 htslib/gcc/1.8 pindel/0.2.5-intel snpeff/4.3q bedtools/2.26.0
touch ${pair_id}.pindel.config
for i in *.bam; do
    sname=`echo "$i" |cut -f 1 -d '.'`
    echo -e "${i}\t400\t${sname}" >> ${pair_id}.pindel.config
done

pindel -T $NPROC -f ${reffa} -i ${pair_id}.pindel.config -o ${pair_id}.pindel_out --RP
pindel2vcf -P ${pair_id}.pindel_out -r ${reffa} -R HG38 -d ${genomefiledate} -v pindel.vcf
cat pindel.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].AD[1] >= 10 )" | bgzip > pindel.vcf.gz
tabix pindel.vcf.gz
bash $baseDir/norm_annot.sh -r ${index_path} -p pindel -v pindel.vcf.gz
perl $baseDir/parse_pindel.pl ${pair_id} pindel.norm.vcf.gz
java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.indel.vcf |bgzip > ${pair_id}.pindel_indel.vcf.gz

if [[ -a $idtbed ]]
then
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.dup.vcf | bedtools intersect -header -b ${idtbed} -a stdin | bgzip > ${pair_id}.pindel_tandemdup.vcf.gz
else
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.dup.vcf | bgzip > ${pair_id}.pindel_tandemdup.vcf.gz
fi

java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.sv.vcf | bgzip > ${pair_id}.pindel_sv.vcf.gz
