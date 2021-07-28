#!/bin/bash
#gatkrunner.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-b  --BAM File"
  echo "-p  --Prefix for output file name"
  echo "-a  --Algorithm/Command"
  echo "Example: bash hisat.sh -p prefix -r GRCh38 -b File.bam"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]] || [[ -z $index_path ]]
then
    usage
fi

NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

dbsnp="${index_path}/dbSnp.gatk4.vcf.gz"
if [[ ! -f "${index_path}/dbSnp.gatk4.vcf.gz" ]]
then
    echo "Missing dbSNP File: ${index_path}/dbSnp.gatk4.vcf.gz"
    usage
fi
reffa="${index_path}/genome.fa"
if [[ ! -f "${index_path}/genome.fa" ]]
then
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load gatk/4.1.4.0 samtools/gcc/1.8
fi
which samtools
samtools index -@ $NPROC ${sbam}

gatk --java-options "-Xmx32g" BaseRecalibrator -I ${sbam} --known-sites $dbsnp -R ${reffa} -O ${pair_id}.recal_data.table --use-original-qualities
gatk --java-options "-Xmx32g" ApplyBQSR -I ${sbam} -R ${reffa} -O ${pair_id}.final.bam --use-original-qualities -bqsr ${pair_id}.recal_data.table
samtools index -@ $NPROC ${pair_id}.final.bam
