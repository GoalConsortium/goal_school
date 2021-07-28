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
while getopts :r:l:n:c:b:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        b) sbam=$OPTARG;;
        n) normal=$OPTARG;;
	c) capture=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

# Check for mandatory options
if [[ -z $sbam ]] || [[ -z $index_path ]]; then
    usage
fi

if [[ -z $isdocker ]]
then
    export PATH=/project/shared/bicf_workflow_ref/seqprg/bin:$PATH
fi
bedopt=''
if [[ -n $capture ]]
then
    bedopt="-e $capture"
fi
   
if [[ -n $normal ]]
then
    msisensor-pro msi -d ${index_path}/microsatellites.list -n $normal -t $sbam -o ${pair_id}.msi $bedopt
else
    # -M ${index_path}/msi_tumoronly_model 
    msisensor-pro pro -d ${index_path}/microsatellites.list_baseline -t ${sbam} -o ${pair_id}.msi $bedopt
fi
