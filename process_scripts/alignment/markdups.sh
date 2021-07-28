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
while getopts :a:b:r:p:h opt
do
    case $opt in
        a) algo=$OPTARG;;
        b) sbam=$OPTARG;;
        p) pair_id=$OPTARG;;
	r) index_path=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $sbam ]]; then
    usage
fi

NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

if [[ -z $index_path ]]
then
    index_path="/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref"
fi

if [[ $index_path == *project* ]]
then
    testexe='/project/shared/bicf_workflow_ref/seqprg/bin'
else
    testexe='/usr/local/bin'
fi

baseDir="`dirname \"$0\"`"

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load picard/2.10.3
fi

if [ $algo == 'samtools' ]
then
    if [[ -z $isdocker ]]
    then
	module load samtools/gcc/1.8
    fi
    samtools markdup -s --output-fmt BAM -@ $NPROC sort.bam ${pair_id}.dedup.bam
    touch ${pair_id}.dedup.stat.txt
elif [ $algo == 'picard' ]
then
    java -XX:ParallelGCThreads=$NPROC -Djava.io.tmpdir=./ -Xmx16g  -jar $PICARD/picard.jar MarkDuplicates I=${sbam} O=${pair_id}.dedup.bam M=${pair_id}.dedup.stat.txt
elif [ $algo == 'picard_umi' ]
then
    java -XX:ParallelGCThreads=$NPROC -Djava.io.tmpdir=./ -Xmx16g  -jar $PICARD/picard.jar MarkDuplicates BARCODE_TAG=RX I=${sbam} O=${pair_id}.dedup.bam M=${pair_id}.dedup.stat.txt
elif [ $algo == 'fgbio_umi' ]   
then
    if [[ -z $isdocker ]]
    then
	module load fgbio bwakit/0.7.15 bwa/intel/0.7.17 samtools/gcc/1.8
    fi
    samtools index -@ $NPROC ${sbam}
    fgbio --tmp-dir ./ GroupReadsByUmi -s identity -i ${sbam} -o group.bam --family-size-histogram ${pair_id}.umihist.txt -e 0 -m 0
    fgbio --tmp-dir ./ CallMolecularConsensusReads -i group.bam -p consensus -M 1 -o ${pair_id}.consensus.bam -S ':none:'
    samtools index -@ $NPROC ${pair_id}.consensus.bam
    samtools fastq -1 ${pair_id}.consensus.R1.fastq -2 ${pair_id}.consensus.R2.fastq ${pair_id}.consensus.bam
    gzip ${pair_id}.consensus.R1.fastq
    gzip ${pair_id}.consensus.R2.fastq
    bwa mem -M -C -t $NPROC -R "@RG\tID:${pair_id}\tLB:tx\tPL:illumina\tPU:barcode\tSM:${pair_id}" ${index_path}/genome.fa ${pair_id}.consensus.R1.fastq.gz ${pair_id}.consensus.R2.fastq.gz > out.sam
    if [[ -f ${index_path}/genome.fa.alt ]]
    then
	k8 ${testexe}/bwa-postalt.js -p tmphla ${index_path}/genome.fa.alt out.sam | samtools view -1 - > ${pair_id}.consensus.bam
    else
	samtools view -1 out.sam > ${pair_id}.consensus.bam
    fi
    samtools sort --threads $NPROC -o ${pair_id}.dedup.bam ${pair_id}.consensus.bam
    samtools sort --threads $NPROC -o ${pair_id}.group.bam group.bam
    samtools index -@ $NPROC ${pair_id}.group.bam
else
    cp ${sbam} ${pair_id}.dedup.bam    
fi
if [[ -z $isdocker ]]
then
    module load samtools/gcc/1.8
fi
samtools index -@ $NPROC ${pair_id}.dedup.bam
