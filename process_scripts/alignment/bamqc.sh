#!/bin/bash
#trimgalore.sh

usage() {
    echo "-h Help documentation for hisat.sh"
    echo "-r  --Reference Genome: GRCh38 or GRCm38"
    echo "-b  --BAM File"
    echo "-n  --NucType"
    echo "-p  --Prefix for output file name"
    echo "-c  --Capture Bedfile"
    echo "-d  --RemoveDuplicates 1=yes, 0=no default=no"
    echo "Example: bash bamqc.sh -p prefix -r /project/shared/bicf_workflow_ref/human/GRCh38 -b SRR1551047.bam  -n dna -c target.bed"
    exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:b:c:n:p:u:e:s:d:x:y:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        b) sbam=$OPTARG;;
        c) bed=$OPTARG;;
        n) nuctype=$OPTARG;;
        p) pair_id=$OPTARG;;
	    u) user=$OPTARG;;
	    e) version=$OPTARG;;
	    s) skiplc=1;;
	    d) dedup=$OPTARG;;
	    x) cpus=$OPTARG;;
	    y) memory=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
#if [[ -z $pair_id ]] || [[ -z $sbam ]]; then
#    usage
#fi

if [[ -z $version ]]
then
    version='NA'
fi

parseopt=""
if [[ -n $user ]]
then
    parseopt=" -u $USER"
fi
if [[ -n $version ]]
then
    parseopt="$parseopt -e $version"
fi
if [[ -f $index_path/reference_info.txt ]]
then
    parseopt="$parseopt -r $index_path"
fi


tmpdir=`pwd`

samtools flagstat ${sbam} > ${pair_id}.flagstat.txt
fastqc -f bam ${sbam}
baseDir="`dirname \"$0\"`"

if [[ $dedup == 1 ]]
then
    mv $sbam ori.bam
    samtools view -@ $cpus -F 1024 -b -o ${sbam} ori.bam
fi

if [[ $nuctype == 'dna' ]]
then
    bedtools sort -faidx $index_path/genome.fa.fai -i ${bed} > panel.sorted.bed
    bedtools coverage -sorted -g  ${index_path}/genomefile.txt -a panel.sorted.bed -b ${sbam} -hist > ${pair_id}.covhist.txt
    grep ^all ${pair_id}.covhist.txt >  ${pair_id}.genomecov.txt
    perl $baseDir/calculate_depthcov.pl ${pair_id}.covhist.txt
    if [[ -z $skiplc ]]
    then
	samtools view -@ $cpus -b -L ${bed} -o ${pair_id}.ontarget.bam ${sbam}
	samtools index -@ $cpus ${pair_id}.ontarget.bam
	samtools flagstat  ${pair_id}.ontarget.bam > ${pair_id}.ontarget.flagstat.txt
	java -Xmx${memory}g -Djava.io.tmpdir=${tmpdir} -XX:ParallelGCThreads=$cpus -jar $PICARD/picard.jar EstimateLibraryComplexity BARCODE_TAG=RG I=${sbam} OUTPUT=${pair_id}.libcomplex.txt TMP_DIR=${tmpdir}
    fi
    perl $baseDir/sequenceqc_dna.pl $parseopt ${pair_id}.genomecov.txt
else
    perl $baseDir/sequenceqc_rna.pl $parseopt ${pair_id}.flagstat.txt
fi
