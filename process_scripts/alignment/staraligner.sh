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


star_opt=$fqs
fqarray=($fqs)


if [[ -z $isdocker ]]
then
    module load star/2.7.3a
fi

STAR --genomeDir ${index_path}/star_index/ --readFilesIn $star_opt --readFilesCommand zcat --genomeLoad NoSharedMemory --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMheaderCommentFile COfile.txt --outSAMheaderHD @HD VN:1.4 SO:coordinate --outSAMunmapped Within --outFilterType BySJout --outSAMattributes NH HI AS NM MD --outSAMstrandField intronMotif --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM --sjdbScore 1 --limitBAMsortRAM 60000000000 --outFileNamePrefix out

mv outLog.final.out ${pair_id}.alignerout.txt
mv outAligned.sortedByCoord.out.bam ${pair_id}.bam
samtools index -@ $NPROC ${pair_id}.bam
