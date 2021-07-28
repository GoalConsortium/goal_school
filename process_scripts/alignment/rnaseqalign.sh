#!/bin/bash
#rnaseqalign.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-x  --FastQ R1"
  echo "-y  --FastQ R2"
  echo "-a  --Method: hisat or star"
  echo "-p  --Prefix for output file name"
  echo "-u  --UMI sequences are in FQ Read Name"
  echo "Example: bash rnaseqalign.sh -a hisat -p prefix -u -r /project/shared/bicf_workflow_ref/human/GRCh38 -x SRR1551047_1.fastq.gz  -y SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:x:y:p:hu opt
do
    case $opt in
        r) index_path=$OPTARG;;
        x) fq1=$OPTARG;;
        y) fq2=$OPTARG;;
	a) algo=$OPTARG;;
	u) umi=1;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]]
then
    usage
fi

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load  samtools/1.6 picard/2.10.3
fi
baseDir="`dirname \"$0\"`"
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ -f $fq1 ]]
then
    fqs="$fq1"
    if [[ -f $fq2 ]]
    then
	diff $fq1 $fq2 > difffile
	if [[ -s difffile ]]
	then
	    fqs+=" $fq2"
	fi
    fi
else
    fqs=''
    i=0
    numfq=${#fqs[@]}
    while [[ $i -le $numfq ]]
    do
	fqs="$fqs $1"
	i=$((i + 1))
	shift 1
    done
fi
numfq=0
for k in $fqs
do
    numfq=$((numfq + 1))
done

hisat_opt=''
fqarray=($fqs)

if [[ $numfq == 2 ]]
then
    hisat_opt="-1 ${fqarray[0]} -2 ${fqarray[1]}"
else
    hisat_opt="-U ${fqarray[0]}"
fi

if [[ -z $isdocker ]]
then
    module load hisat2/2.1.0-intel
fi

idx="${index_path}/genome"
if [[ -d ${index_path}/hisat_index ]]
then
    idx ="${index_path}/hisat_index/genome"
fi

hisat2 -p $NPROC --rg-id ${pair_id} --rg LB:tx --rg PL:illumina --rg PU:barcode --rg SM:${pair_id} --dta -x ${idx} $hisat_opt -S out.sam --summary-file ${pair_id}.alignerout.txt

if [[ $umi == 1 ]]
then
    python ${baseDir}/add_umi_sam.py -s out.sam -o output.bam
else
    samtools view -1 --threads $NPROC -o output.bam out.sam
fi
samtools sort -n -@ $NPROC -O BAM -o output.dups.bam output.bam
java -Djava.io.tmpdir=./ -Xmx4g  -jar $PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.dups.bam O=${pair_id}.bam
samtools index -@ $NPROC ${pair_id}.bam
