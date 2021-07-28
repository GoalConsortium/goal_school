#!/bin/bash
#annotvcf.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-p  --Prefix for output file name"
  echo "-v  --VCF File"
  echo "Example: bash hisat.sh -p prefix -r /path/GRCh38 -v vcffile"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:v:g:p:h opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
	v) unionvcf=$OPTARG;;
	g) snpeffgeno=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load bedtools/2.26.0  snpeff/4.3q
fi

if [[ $index_path == '/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref' ]]
then
    index_path='/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref'
fi
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi

if  [[ -f "${index_path}/gnomad.txt.gz" ]] 
then
    tabix -f ${unionvcf}
    bcftools annotate -Oz -a ${index_path}/gnomad.txt.gz -h ${index_path}/gnomad.header -c CHROM,POS,REF,ALT,GNOMAD_HOM,GNOMAD_AF,AF_POPMAX,GNOMAD_HG19_VARIANT -o ${pair_id}.gnomad.vcf.gz ${unionvcf}
    tabix -f ${pair_id}.gnomad.vcf.gz
    unionvcf=${pair_id}.gnomad.vcf.gz
fi
if  [[ -f "${index_path}/oncokb_hotspot.txt.gz" ]] 
then
    bcftools annotate -Oz -a ${index_path}/oncokb_hotspot.txt.gz -o ${pair_id}.oncohotspot.vcf.gz -h ${index_path}/oncokb_hotspot.header -c CHROM,FROM,TO,OncoKB_REF,OncoKB_ALT,Gene,OncoKB_ProteinChange,OncoKB_AF,OncoTree_Tissue,OncoTree_MainType,OncoTree_Code,OncoKBHotspot $unionvcf
    tabix -f ${pair_id}.oncohotspot.vcf.gz
    unionvcf=${pair_id}.oncohotspot.vcf.gz
fi
if [[ -f "${index_path}/repeat_regions.bed.gz" ]]
then
    bcftools annotate -Oz -a ${index_path}/repeat_regions.bed.gz -o ${pair_id}.repeat.vcf.gz --columns CHROM,FROM,TO,RepeatType -h ${index_path}/RepeatType.header ${unionvcf}
    unionvcf=${pair_id}.repeat.vcf.gz
fi

acom="java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-downstream -no-upstream -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config $snpeffgeno ${unionvcf}"

if [[ -f ${index_path}/dbSnp.vcf.gz ]]
then
    acom+=" | java -jar $SNPEFF_HOME/SnpSift.jar annotate -id ${index_path}/dbSnp.vcf.gz -"
fi
if [[ -f ${index_path}/clinvar.vcf.gz ]]
then
    acom+=" | java -jar $SNPEFF_HOME/SnpSift.jar annotate -info CLNSIG,CLNDSDB,CLNDSDBID,CLNDBN,CLNREVSTAT,CLNACC ${index_path}/clinvar.vcf.gz -"
fi
if [[ -f ${index_path}/cosmic.vcf.gz ]]
then
    acom+=" | java -jar $SNPEFF_HOME/SnpSift.jar annotate -info LEGACY_ID,CNT ${index_path}/cosmic.vcf.gz -"
fi
if [[ -f ${index_path}/dbNSFP.txt.gz ]]
then
    acom+=" | java -Xmx10g -jar $SNPEFF_HOME/SnpSift.jar dbnsfp -v -db ${index_path}/dbNSFP.txt.gz -"
fi

acom+=" | bgzip > ${pair_id}.annot.vcf.gz "

eval $acom 

tabix ${pair_id}.annot.vcf.gz
