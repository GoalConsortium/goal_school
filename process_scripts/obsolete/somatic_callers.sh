#!/bin/bash
#run_somatic_caller.sh


usage(){
  echo "-h --Help documentation for run_somatic_caller.sh"
  echo "-a --Somatic Workflow Method: strelka2, virmid, speedseq, mutect, varscan, shimmer"
  echo "-n --Normal"
  echo "-t --Tumor"
  echo "Example: bash somatic_callers.sh -a strelka2 -n ORD1_N_panel1385 -t ORD1_T_panel1385"
  exit 1
}

OPTIND=1 # Reset OPTIND

while getopts :n:t:a:h opt
do
    case $opt in 
      n) normal=$OPTARG;;
      t) tumor=$OPTARG;;
      a) algo=$OPTARG;;
      h) usage;;
    esac
done

shift $(($OPTIND -1))


#Check for mandatory options
if [[ -z $normal ]] || [[ -z $tumor ]] || [[ -z $algo ]]; then
  usage
fi 

if [[ -z $SLURM_CPUS_ON_NODE ]] 
  then 
    SLURM_CPUS_ON_NODE=1
fi

index_path=/project/shared/bicf_workflow_ref/human/GRCh38

genome_reference=${index_path}/genome.fa
cosmic_reference=${index_path}/cosmic.vcf.gz
dbSnp_reference=${index_path}/dbSnp.vcf.gz
target_bed=${index_path}/UTSWV2.bed

if [ $algo == 'strelka2' ]
  then
    module load strelka/2.8.3 manta/1.2.0
    manta_analysisPath=MantaAnalysisPath 
    strelka_analysisPath=StrelkaAnalysisPath
    mkdir ${manta_analysisPath}
    mkdir ${strelka_analysisPath}
    /cm/shared/apps/manta/1.2.0/bin/configManta.py --normalBam ${normal}.final.bam --tumorBam ${tumor}.final.bam --referenceFasta ${genome_reference} --runDir ${manta_analysisPath}
    ${manta_analysisPath}/runWorkflow.py -m local -j 8
    /cm/shared/apps/strelka/2.8.3/bin/configureStrelkaSomaticWorkflow.py --normalBam ${normal}.final.bam --tumorBam ${tumor}.final.bam --referenceFasta ${genome_reference} --targeted --indelCandidates ${manta_analysisPath}/results/variants/candidateSmallIndels.vcf.gz --runDir ${strelka_analysisPath}
    ${strelka_analysisPath}/runWorkflow.py -m local -j 8
    bcftools norm -c s -f ${reffa} -w 10 -O z -o ${pair_id}.strelka2.vcf.gz strelka/results/variants/variants.vcf.gz
fi

if [ $algo == 'virmid' ]
  then 
    module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 virmid/1.2 vcftools/0.1.14
    virmid -R ${genome_reference} -D ${tumor}.final.bam -N ${normal}.final.bam -s ${cosmic_reference} -t $SLURM_CPUS_ON_NODE -M 2000 -c1 10 -c2 10
    perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/addgt_virmid.pl ${tumor}.final.bam.virmid.som.passed.vcf
    perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/addgt_virmid.pl ${tumor}.final.bam.virmid.loh.passed.vcf
    vcf-concat *gt.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((NDP >= 10) & (DDP >= 10))' | perl -pe 's/TUMOR/'${tumor}'/' | perl -pe 's/NORMAL/'${normal}'/g' |bedtools intersect -header -a stdin -b ${target_bed} |bgzip > ${tumor}_${normal}.virmid.vcf.gz
fi

if [ $algo == 'speedseq' ]
  then 
    module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 speedseq/20160506 bcftools/intel/1.3 vcftools/0.1.14
    speedseq somatic -q 10 -w ${target_bed} -t $SLURM_CPUS_ON_NODE -o ${tumor}.sssom ${genome_reference} ${normal}.final.bam ${tumor}.final.bam
    vcf-annotate -H -n --fill-type ${tumor}.sssom.vcf.gz | java -jar $SNPEFF_HOME/SnpSift.jar filter --pass '((QUAL >= 10) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/'${tumor}'/' | perl -pe 's/NORMAL/'${normal}'/g' |bgzip > ${tumor}_${normal}.sspanel.vcf.gz
fi

if [ $algo == 'mutect' ]
  then
    module load parallel python/2.7.x-anaconda gatk/3.8  bcftools/intel/1.3 bedtools/2.25.0 snpeff/4.2 vcftools/0.1.14
    cut -f 1 /project/shared/bicf_workflow_ref/human/GRCh38/genomefile.5M.txt | parallel --delay 2 -j 10 "java -Xmx20g -jar $GATK_JAR -R ${genome_reference} -D ${dbSnp_reference} -T MuTect2 -stand_call_conf 30 -stand_emit_conf 10.0 -A FisherStrand -A QualByDepth -A VariantType -A DepthPerAlleleBySample -A HaplotypeScore -A AlleleBalance -I:tumor ${tumor}.final.bam -I:normal ${normal}.final.bam --cosmic ${cosmic} -o ${tumor}.{}.mutect.vcf -L {}"
    vcf-concat ${tumor}*.vcf | vcf-sort | vcf-annotate -n --fill-type | java -jar $SNPEFF_HOME/SnpSift.jar filter -p '((FS <= 60) & GEN[*].DP >= 10)' | perl -pe 's/TUMOR/'${tumor}'/' | perl -pe 's/NORMAL/'${normal}'/g' |bgzip > ${tumor}_${normal}.mutect.vcf.gz
fi

if [ $algo == 'varscan' ]
  then
    module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3 samtools/intel/1.3 VarScan/2.4.2 speedseq/20160506 vcftools/0.1.14
    sambamba mpileup -L ${target_bed} -t $SLURM_CPUS_ON_NODE ${tumor}.final.bam --samtools "-C 50 -f ${genome_reference}"  > t.mpileup
    sambamba mpileup -L ${target_bed} -t $SLURM_CPUS_ON_NODE ${normal}.final.bam --samtools "-C 50 -f ${genome_reference}"  > n.mpileup
    VarScan somatic n.mpileup t.mpileup ${tumor}.vscan --output-vcf 1
    VarScan copynumber n.mpileup t.mpileup ${tumor}.vscancnv 
    vcf-concat ${tumor}.vscan*.vcf | vcf-sort | vcf-annotate -n --fill-type -n | java -jar $SNPEFF_HOME/SnpSift.jar filter '((exists SOMATIC) & (GEN[*].DP >= 10))' | perl -pe 's/TUMOR/'${tumor}'/' | perl -pe 's/NORMAL/'${normal}'/g' |bedtools intersect -header -a stdin -b ${target_bed} |bgzip > ${tumor}_${normal}.varscan.vcf.gz
fi

if [ $algo == 'shimmer' ]
  then
    module load python/2.7.x-anaconda bedtools/2.25.0 snpeff/4.2 bcftools/intel/1.3  shimmer/0.1.1 vcftools/0.1.14
    shimmer.pl --minqual 25 --ref ${genome_reference} ${normal}.final.bam ${tumor}.final.bam --outdir shimmer 2> shimmer.err
    perl /project/PHG/PHG_Clinical/clinseq_workflows/scripts/add_readct_shimmer.pl
    vcf-annotate -n --fill-type shimmer/somatic_diffs.readct.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter '(GEN[*].DP >= 10)' | perl -pe 's/TUMOR/'${tumor}'/' | perl -pe 's/NORMAL/'${normal}'/g' | bedtools intersect -header -a stdin -b ${target_bed} | bgzip > ${tumor}_${normal}.shimmer.vcf.gz
fi
