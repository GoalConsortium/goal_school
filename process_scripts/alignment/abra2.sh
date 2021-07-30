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
while getopts :r:b:p:t:c:m:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	    t) tbed=$OPTARG;;
	    c) cpus=$OPTARG;;
	    m) memory=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]] || [[ -z $index_path ]]
then
    usage
fi

abrajar=/usr/local/bin/abra2.jar

opt=''
if [ -n "$tbed" ]
then
    opt="--targets $tbed"
fi

samtools index -@ $cpus ${sbam}
mkdir tmpdir
java -Xmx${memory}G -jar ${abrajar} --in ${sbam} --out ${pair_id}.abra2.bam --ref ${index_path}/genome.fa --threads $cpus $opt --tmpdir tmpdir --mbq 150 --mnf 5 --mer 0.05 > abra.log
samtools index -@ $cpus ${pair_id}.abra2.bam
