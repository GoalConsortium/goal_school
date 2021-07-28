#!/bin/bash
#indexbams.sh

usage() {
  echo "-h  --Help documentation for markdups.sh"
  echo "Example: bash indexbams.sh"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :h opt
do
    case $opt in
        h) usage;;
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
    module load samtools/1.6
fi

for i in *.bam; do
    samtools index -@ $NPROC ${i}
done 
