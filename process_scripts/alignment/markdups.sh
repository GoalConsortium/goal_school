#!/bin/bash
#hisat.sh

usage() {
  echo "-h  --Help documentation for markdups.sh"
  echo "-m  --Mark Duplication Method: sambamba, samtools, picard, picard_umi, fgbio_umi, null; default is null"
  echo "-b  --BAM file"
  echo "-p  --Prefix for output file name"
  echo "Example: bash markdups.sh -p prefix -b file.bam -a picard"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :a:b:p:r:c:m:h opt
do
    case $opt in
        a) algo=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	    r) index_path=$OPTARG;;
	    c) cpus=$OPTARG;;
	    m) memory=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]; then
    usage
fi

if [[ -z $index_path ]]
then
    index_path="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
fi

testexe='/usr/local/bin'

if [ $algo == 'samtools' ]
then
    samtools markdup -s --output-fmt BAM -@ $cpus sort.bam ${pair_id}.dedup.bam
    touch ${pair_id}.dedup.stat.txt
elif [ $algo == 'picard' ]
then
    java -XX:ParallelGCThreads=$cpus -Djava.io.tmpdir=./ -Xmx${memory}g  -jar $PICARD/picard.jar MarkDuplicates I=${sbam} O=${pair_id}.dedup.bam M=${pair_id}.dedup.stat.txt
elif [ $algo == 'picard_umi' ]
then
    java -XX:ParallelGCThreads=$cpus -Djava.io.tmpdir=./ -Xmx${memory}g  -jar $PICARD/picard.jar MarkDuplicates BARCODE_TAG=RX I=${sbam} O=${pair_id}.dedup.bam M=${pair_id}.dedup.stat.txt
elif [ $algo == 'fgbio_umi' ]   
then
    samtools index -@ $cpus ${sbam}
    fgbio --tmp-dir ./ GroupReadsByUmi -s identity -i ${sbam} -o group.bam --family-size-histogram ${pair_id}.umihist.txt -e 0 -m 0
    fgbio --tmp-dir ./ CallMolecularConsensusReads -i group.bam -p consensus -M 1 -o ${pair_id}.consensus.bam -S ':none:'
    samtools index -@ $cpus ${pair_id}.consensus.bam
    samtools fastq -1 ${pair_id}.consensus.R1.fastq -2 ${pair_id}.consensus.R2.fastq ${pair_id}.consensus.bam
    gzip ${pair_id}.consensus.R1.fastq
    gzip ${pair_id}.consensus.R2.fastq
    bwa mem -M -C -t $cpus -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${pair_id}.consensus.R1.fastq.gz ${pair_id}.consensus.R2.fastq.gz > out.sam
    if [[ -f ${index_path}/genome.fa.alt ]]
    then
	k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam | samtools view -1 - > ${pair_id}.consensus.bam
    else
	samtools view -1 out.sam > ${pair_id}.consensus.bam
    fi
    samtools sort --threads $cpus -o ${pair_id}.dedup.bam ${pair_id}.consensus.bam
    samtools sort --threads $cpus -o ${pair_id}.group.bam group.bam
    samtools index -@ $cpus ${pair_id}.group.bam
else
    cp ${sbam} ${pair_id}.dedup.bam    
fi
samtools index -@ $cpus ${pair_id}.dedup.bam
