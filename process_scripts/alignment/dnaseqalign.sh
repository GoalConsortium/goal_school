#!/bin/bash
#dnaseqalign.sh

set -e

usage() {
  echo "-h Help documentation for dnaseqalign.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-x  --FastQ R1"
  echo "-y  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "-u  --UMI"
  echo "Example: bash dnaseqalign.sh -p prefix -u 1 -r /project/shared/bicf_workflow_ref/human/GRCh38 -x SRR1551047_1.fastq.gz  -y SRR1551047_2.fastq.gz"
  exit 1
}

OPTIND=1 # Reset OPTIND
#while getopts :r:x:y:g:p:a:uh opt
while getopts :r:x:y:g:p:a:c:m:uh opt
do
    case $opt in
        r) index_path=$OPTARG;;
        x) fq1=$OPTARG;;
        y) fq2=$OPTARG;;
	    g) read_group=$OPTARG;;
        p) pair_id=$OPTARG;;
	    a) aligner=$OPTARG;;
	    c) cpus=$OPTARG;;
	    m) memory=$OPTARG;;
	    u) umi='umi';;
        h) usage;;
    esac
done

shift $(($OPTIND -1))


# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $fq1 ]]; then
    usage
fi


if [[ -z $read_group ]]
then
    read_group=$pair_id
fi

testexe='/usr/local/bin'

scriptDir="`dirname \"$0\"`"

if cmp -s "$fq1" "$fq2"; then
    file_opt="${fq1}"
else
    file_opt="${fq1} ${fq2}"
fi

echo bwa mem -M -t $cpus -R "@RG\tID:${read_group}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${read_group}" ${index_path}/genome.fa $file_opt 
ls -rlth ${index_path}/genome.fa*


bwa mem -M -t $cpus -R "@RG\tID:${read_group}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${read_group}" ${index_path}/genome.fa $file_opt > out.sam

if [[ $umi == 'umi' ]] && [[ -f "${index_path}/genome.fa.alt" ]]
then
    k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam | python ${scriptDir}/add_umi_sam.py -s - -o output.unsort.bam
elif [[ -f "${index_path}/genome.fa.alt" ]]
then
    k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam| samtools view -1 - > output.unsort.bam
elif [[ $umi == 'umi' ]]
then
    python ${scriptDir}/add_umi_sam.py -s out.sam -o output.unsort.bam
else
    samtools view -1 -o output.unsort.bam out.sam
fi

which samtools

echo samtools sort -n --threads $cpus -o output.dups.bam output.unsort.bam

samtools sort -n --threads $cpus -o output.dups.bam output.unsort.bam
java -Djava.io.tmpdir=./ -Xmx${memory}g  -jar $PICARD/picard.jar FixMateInformation ASSUME_SORTED=TRUE SORT_ORDER=coordinate ADD_MATE_CIGAR=TRUE I=output.dups.bam O=${pair_id}.bam
samtools index ${pair_id}.bam
