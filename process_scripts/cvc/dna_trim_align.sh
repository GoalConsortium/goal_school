#!/bin/bash
#trimgalore.sh

usage() {
  echo "-h Help documentation for hisat.sh"
  echo "-a  --FastQ R1"
  echo "-b  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash trimgalore.sh -p prefix -a SRR1551047_1.fastq.gz  -b SRR1551047_2.fastq.gz"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:g:q:b:p:ufh opt
do
    case $opt in
        p) pair_id=$OPTARG;;
	f) filter=1;;
        r) index_path=$OPTARG;;
	u) umi='umi';;
	g) read_group=$OPTARG;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"

fqs=''
i=0
numfq=$#

while [[ $i -le $numfq ]]
do
    fqs="$fqs $1"
    i=$((i + 1))
    shift 1
done

# Check for mandatory options
if [[ -z $pair_id ]]
then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi
if [[ -z $read_group ]]
then
    read_group=$pair_id
fi
if [[ $index_path == *project* ]]
then
    testexe='/project/shared/bicf_workflow_ref/seqprg/bin'
else
    testexe='/usr/local/bin'
fi

source /etc/profile.d/modules.sh
module load trimgalore/0.6.4 cutadapt/1.9.1 bwakit/0.7.15 samtools/gcc/1.8 picard/2.10.3

threads=`expr $NPROC / 2`

trim_galore --cores 4 --paired -q 25 --illumina --gzip --length 35 ${fqs}

if [[ ${filter} == 1 ]]
then
      perl $baseDir/parse_trimreport.pl ${pair_id}.trimreport.txt *trimming_report.txt
fi

bwa mem -M -t $threads -R "@RG\tID:${read_group}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${read_group}" ${index_path}/genome.fa *.fq.gz > out.sam

if [[ $umi == 'umi' ]] && [[ -f "${index_path}/genome.fa.alt" ]]
then
    k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam | python ${baseDir}/add_umi_sam.py -s - -o output.unsort.bam
elif [[ -f "${index_path}/genome.fa.alt" ]]
then
    k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam| samtools view -1 - > output.unsort.bam
elif [[ $umi == 'umi' ]]
then
    python ${baseDir}/add_umi_sam.py -s out.sam -o output.unsort.bam
else
    samtools view -1 -o output.unsort.bam out.sam
fi

samtools sort -n --threads $threads -o output.dups.bam output.unsort.bam
java -Djava.io.tmpdir=./ -Xmx4g  -jar $PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.dups.bam O=${pair_id}.bam
samtools index ${pair_id}.bam
