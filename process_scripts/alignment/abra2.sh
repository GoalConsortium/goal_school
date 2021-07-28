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
while getopts :r:a:c:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	c) tbed=$OPTARG;;
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

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load abra2/2.18 samtools/gcc/1.8
    abrajar=/cm/shared/apps/abra2/lib/abra2.jar
else
    abrajar=/usr/local/bin/abra2.jar
fi
ioopt="--in ${sbam} --out ${pair_id}.abra2.bam"
opt=''
if [ -n "$tbed" ]
then
    opt="--targets $tbed"
fi

which samtools
samtools index -@ $NPROC ${sbam}
mkdir tmpdir
java -Xmx16G -jar ${abrajar} ${ioopt} --ref ${index_path}/genome.fa --threads $NPROC $opt --tmpdir tmpdir --mbq 150 --mnf 5 --mer 0.05 > abra.log
samtools index -@ $NPROC ${pair_id}.abra2.bam
