#!/bin/bash
#indexbams.sh

usage() {
  echo "-h  --Help documentation for markdups.sh"
  echo "Example: bash indexbams.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"

  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:p:h opt
do
    case $opt in
        h) usage;;
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	b) bam=$OPTARG;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options

NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

baseDir="`dirname \"$0\"`"
if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load igvtools/2.3.71 samtools/1.6
exit
samtools index  -@ $NPROC $bam
igvtools count -z 5 $bam ${pair_id}.tdf ${index_path}/igv/human.genome
