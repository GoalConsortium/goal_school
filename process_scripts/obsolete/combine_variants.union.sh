#!/bin/bash
#union.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }
shift $(($OPTIND -1))
baseDir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
module load gatk/3.7 python/2.7.x-anaconda bedtools/2.26.0 snpeff/4.3q samtools/1.6 vcftools/0.1.14

HS=*.hotspot.vcf.gz
list1=`ls *vcf.gz|grep -v hotspot`
list2=`ls *vcf.gz|grep -v hotspot`
varlist=''
calllist=''
for i in *.vcf.gz; do
    EXT="${i#*.}"
    CALL="${EXT%%.*}"
    calllist="$calllist $CALL"
    tabix $i
    if [[ $i == $HS ]]
    then
	bedtools multiinter -i $list1 |cut -f 1,2,3 |bedtools intersect -header -v -a $i -b stdin |bgzip > hotspot.nooverlap.vcf.gz
	tabix hotspot.nooverlap.vcf.gz
	list2="$list2 hotspot.nooverlap.vcf.gz"
	varlist="$varlist --variant:$CALL hotspot.nooverlap.vcf.gz"
    else 
	varlist="$varlist --variant:$CALL $i"
    fi
done 

bedtools multiinter -i $list2 -names $calllist | cut -f 1,2,3,5 | bedtools sort -i stdin | bedtools merge -c 4 -o distinct >  ${pair_id}_integrate.bed

priority='ssvar'
if [[ *.platypus.vcf.gz ]]
then
    priority="$priority,platypus"
fi
priority="$priority,sam,gatk"
if [[ -f *.hotspot.vcf.gz ]]
then
    priority="$priority,hotspot"
fi

java -Xmx32g -jar $GATK_JAR -R ${index_path}/genome.fa -T CombineVariants --filteredrecordsmergetype KEEP_UNCONDITIONAL $varlist -genotypeMergeOptions PRIORITIZE -priority $priority -o ${pair_id}.int.vcf

perl $baseDir/uniform_integrated_vcf.pl ${pair_id}.int.vcf
bgzip ${pair_id}_integrate.bed
tabix ${pair_id}_integrate.bed.gz
bgzip ${pair_id}.uniform.vcf
tabix ${pair_id}.uniform.vcf.gz
bcftools annotate -a ${pair_id}_integrate.bed.gz --columns CHROM,FROM,TO,CallSet -h ${index_path}/CallSet.header ${pair_id}.uniform.vcf.gz | bgzip > ${pair_id}.union.vcf.gz
