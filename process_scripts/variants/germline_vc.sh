#!/bin/bash
#germline_vc.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-a  --Algorithm/Command: gatk, mpileup, speedseq, platypus"
  echo "-t  --RNASeq Data"
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:a:b:q:p:c:m:th opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        a) algo=$OPTARG;;
	    t) rna=1;;
	    b) tbed=$OPTARG;;
	    q) pon=$OPTARG;; 
	    c) cpus=$OPTARG;;
	    m) memory=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi

if [[ -s "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage
fi
if [[ -f $pon ]]
then
    ponopt="--pon $pon"
else
    ponopt='';
fi

fbsplit="${index_path}/genomefile.5M.txt"
cat ${reffa}.fai |cut -f 1 |grep -v decoy |grep -v 'HLA' |grep -v alt |grep -v 'chrUn' |grep -v 'random' > intervals.txt
interval=`cat intervals.txt | perl -pe 's/\n/ -L /g' |perl -pe 's/-L $//'`

if [[ -n $tbed ]]
then
    awk '{print $1":"$2"-"$3}' $tbed > intervals.txt
    interval=$tbed
fi

for i in *.bam; do
    if [[ ! -f ${i}.bai ]]
    then
	samtools index -@ $cpus $i
    fi
done

if [[ $algo == 'mpileup' ]]
then
    threads=`expr $cpus - 10`
    bcftools mpileup --threads $threads -a 'INFO/AD,INFO/ADF,INFO/ADR,FORMAT/DP,FORMAT/SP,FORMAT/AD,FORMAT/ADF,FORMAT/ADR' -Ou -A -d 1000000 -C50 -f ${reffa} *.bam | bcftools call -A --threads 10 -vmO z -o ${pair_id}.vcf.gz
    vcf-annotate -n --fill-type ${pair_id}.vcf.gz | bcftools norm -c s -f ${reffa} -w 10 -O v -o sam.vcf
    java -jar $PICARD/picard.jar SortVcf I=sam.vcf O=${pair_id}.sam.vcf R=${reffa} CREATE_INDEX=TRUE
    bgzip ${pair_id}.sam.vcf
elif [[ $algo == 'fb' ]]
then
    paropt="--delay 1 --jobs 0 --memfree 2G"
    bamlist=''
    for i in *.bam; do
    bamlist="$bamlist --bam ${PWD}/${i}"
    done
    cut -f 1 $fbsplit | parallel ${paropt} "freebayes -f ${index_path}/genome.fa  --min-mapping-quality 0 --min-base-quality 20 --min-coverage 10 --min-alternate-fraction 0.01 -C 3 --use-best-n-alleles 3 -r {} ${bamlist} > fb.{}.vcf"
    vcf-concat fb.*.vcf | vcf-sort | vcf-annotate -n --fill-type | bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.fb.vcf.gz -
elif [[ $algo == 'platypus' ]]
then
    bamlist=`join_by , *.bam`
    Platypus.py callVariants --minMapQual=0 --minReads=3 --mergeClusteredVariants=1 --nCPU=$cpus --bamFiles=${bamlist} --refFile=${reffa} --output=platypus.vcf
    for i in *.bam
    do
	prefix="${i%.bam}"
	sid=`samtools view -H ${i} |grep '^@RG' |perl -pe 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
	perl -pi -e "s/$prefix/$sid/g" platypus.vcf
    done
    vcf-sort platypus.vcf |vcf-annotate -n --fill-type -n |bgzip > platypus.vcf.gz
    tabix platypus.vcf.gz
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.platypus.vcf.gz platypus.vcf.gz
elif [[ $algo == 'gatk' ]]
then
    gatk4_dbsnp="${index_path}/dbSnp.gatk4.vcf.gz"
    if [[ ! -f "${index_path}/dbSnp.gatk4.vcf.gz" ]]
    then
	echo "Missing dbSNP File: ${index_path}/dbSnp.vcf.gz"
	usage
    fi
    user=$USER
    gvcflist=''
    for i in *.bam; do
	prefix="${i%.bam}"
	echo ${prefix}
	gatk --java-options "-Xmx${memory}g" HaplotypeCaller -R ${reffa} -I ${i} -A FisherStrand -A QualByDepth -A DepthPerAlleleBySample -A TandemRepeat --emit-ref-confidence GVCF -G StandardAnnotation -G AS_StandardAnnotation -G StandardHCAnnotation -O haplotypecaller.vcf.gz -L $interval
	java -jar $PICARD/picard.jar SortVcf I=haplotypecaller.vcf.gz O=${prefix}.gatk.g.vcf R=${reffa} CREATE_INDEX=TRUE
	
	gvcflist="$gvcflist -V ${prefix}.gatk.g.vcf"
    done
    gatk --java-options "-Xmx${memory}g" GenomicsDBImport $gvcflist --genomicsdb-workspace-path gendb -L $interval --reader-threads $cpus 
    gatk --java-options "-Xmx${memory}g" GenotypeGVCFs -V gendb://gendb -R ${reffa} -D ${gatk4_dbsnp} -O gatk.vcf -L $interval
    bcftools norm -c s -f ${reffa} -w 10 -O v gatk.vcf | vcf-annotate -n --fill-type | bgzip > ${pair_id}.gatk.vcf.gz
    tabix ${pair_id}.gatk.vcf.gz
elif [ $algo == 'mutect' ]
then
    threads=`expr $cpus / 2`
    bamlist=''
    for i in *.bam; do
	bamlist+="-I ${i} "
    done
    gatk --java-options "-Xmx${memory}g" Mutect2 $ponopt -R ${reffa} ${bamlist} --output ${pair_id}.mutect.vcf -RF AllowAllReadsReadFilter --independent-mates  --tmp-dir `pwd` -L $interval
    vcf-sort ${pair_id}.mutect.vcf | vcf-annotate -n --fill-type | java -jar $SNPEFF_HOME/SnpSift.jar filter -p '(GEN[*].DP >= 10)' | bgzip > ${pair_id}.mutect.vcf.gz
elif [[ $algo == 'strelka2' ]]
then
    opt=''
    if [[ -n $tbed ]]
    then
	if [[ -f "${tbed}.gz" ]]
	then
	    opt="--callRegions ${tbed}.gz"
	else
	    cp $tbed panel.bed
	    bgzip panel.bed
	    tabix panel.bed.gz
	    opt="--callRegions panel.bed.gz"
	fi
    fi
    if [[ $rna == 1 ]]
    then
	mode="--rna"
    else
	mode="--exome"
    fi
    mkdir manta strelka
    gvcflist=''
    for i in *.bam; do
	gvcflist="$gvcflist --bam ${i}"
    done
    configManta.py $gvcflist --referenceFasta ${reffa} $mode --runDir manta
    manta/runWorkflow.py -m local -j $cpus
    if [[ -f manta/results/variants/candidateSmallIndels.vcf.gz ]]
    then
	configureStrelkaGermlineWorkflow.py $gvcflist --referenceFasta ${reffa} $mode --indelCandidates manta/results/variants/candidateSmallIndels.vcf.gz --runDir strelka
    else
	configureStrelkaGermlineWorkflow.py $gvcflist --referenceFasta ${reffa} $mode --runDir strelka
    fi
    strelka/runWorkflow.py -m local -j $cpus
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.strelka2.vcf.gz strelka/results/variants/variants.vcf.gz
fi
