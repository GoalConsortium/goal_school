#!/bin/bash
#svcalling.sh

usage() {
  echo "-h Help documentation for gatkrunner.sh"
  echo "-r  --Path to Reference Genome with the file genome.fa"
  echo "-p  --Prefix for output file name"
  echo "-b  --Bam File"
  echo "-n  --Reference Bam File"
  echo "Example: bash svcalling.sh -p prefix -r /path/GRCh38 -a gatk"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :r:p:b:t:x:c:g:y:n:l:a:hf opt
do
    case $opt in
        r) index_path=$OPTARG;;
        p) pair_id=$OPTARG;;
        t) tumor=$OPTARG;;
        n) normal=$OPTARG;;
	a) method=$OPTARG;;
        x) tid=$OPTARG;;
        y) nid=$OPTARG;;
	f) filter=1;;
	g) snpeffgeno=$OPTARG;;
        b) sbam=$OPTARG;;
	c) tbed=$OPTARG;;
	l) itdbed=$OPTARG;;
        h) usage;;
    esac
done
function join_by { local IFS="$1"; shift; echo "$*"; }

shift $(($OPTIND -1))
baseDir="`dirname \"$0\"`"


# Check for mandatory options
if [[ -z $pair_id ]] || [[ -z $index_path ]]; then
    usage
fi
NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

if [[ -a "${index_path}/genome.fa" ]]
then
    reffa="${index_path}/genome.fa"
    dict="${index_path}/genome.dict"
else 
    echo "Missing Fasta File: ${index_path}/genome.fa"
    usage

fi
if [[ -z $snpeffgeno ]]
then
    snpeffgeno='GRCh38.86'
fi
mkdir -p temp

if [[ -z $tid ]] && [[ -f ${sbam} ]]
then
    tid=`samtools view -H ${sbam} |grep '^@RG' |perl -pe 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
fi

bams=''
touch ${pair_id}.pindel.config
for i in *.bam; do
    bams="$bams $i"
    sid=`samtools view -H ${i} |grep '^@RG' |perl -pe 's/\t/\n/g' |grep ID |cut -f 2 -d ':'`
    echo -e "${i}\t400\t${sid}" >> ${pair_id}.pindel.config
    if [[ $sid =~ "_T_" ]]
    then
	tid=$sid
    fi
    samtools index -@ $NPROC $i
done
bamlist=''
for i in *.bam; do
    bamlist="$bamlist -t ${i}"
done

#RUN DELLY

if  [[ $method == 'delly' ]] || [[ $method == 'svaba' ]] || [[ $method == 'gridss' ]]
then
    if [[ $method == 'delly' ]]
    then
	delly2 call -t BND -o ${pair_id}.delly_translocations.bcf -q 30 -g ${reffa} ${bams} 
	delly2 call -t DUP -o ${pair_id}.delly_duplications.bcf -q 30 -g ${reffa} ${bams}
	delly2 call -t INV -o ${pair_id}.delly_inversions.bcf -q 30 -g ${reffa} ${bams}
	delly2 call -t DEL -o ${pair_id}.delly_deletion.bcf -q 30 -g ${reffa} ${bams}
	delly2 call -t INS -o ${pair_id}.delly_insertion.bcf -q 30 -g ${reffa} ${bams}
	#MERGE DELLY AND MAKE BED
	bcftools concat -a -O v ${pair_id}.delly_duplications.bcf ${pair_id}.delly_inversions.bcf ${pair_id}.delly_translocations.bcf ${pair_id}.delly_deletion.bcf ${pair_id}.delly_insertion.bcf | vcf-sort -t temp | bgzip > ${pair_id}.${method}.svar.vcf.gz
    elif [[ $method == 'svaba' ]]
    then
	svaba run -p $NPROC -G ${reffa} -a ${pair_id} $bamlist
	vcf-concat ${pair_id}.svaba.unfiltered*sv.vcf | perl -pe 's/\.consensus|\.bam//g' | vcf-sort| bgzip > ${pair_id}.${method}.var.vcf.gz
	vcf-concat ${pair_id}.svaba.unfiltered*indel.vcf | perl -pe 's/\.consensus|\.bam//g' | vcf-sort | java -jar $SNPEFF_HOME/SnpSift.jar filter "( SPAN >= 20)" - |bgzip > ${pair_id}.svaba.indel.vcf.gz
	vcf-concat ${pair_id}.${method}.var.vcf.gz ${pair_id}.svaba.indel.vcf.gz | vcf-sort| bgzip > ${pair_id}.${method}.svar.vcf.gz
	rm ${pair_id}.contigs.bam
    elif [[ $method == 'gridss' ]]
    then
	singularity exec --no-home --cleanenv /project/shared/bicf_workflow_ref/seqprg/singularity/gridss_v2.11.1.img gridss.sh --reference ${reffa} -o  ${pair_id}.${method}.svar.vcf.gz -t 8 -a gridss.assembly.bam --workingdir temp --steps All ${bams} --jvmheap 31g
	rm gridss.assembly.bam
    fi
    bash $baseDir/norm_annot.sh -r ${index_path} -p ${pair_id}.${method}.sv -v ${pair_id}.${method}.svar.vcf.gz -s    
    java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].DP >= 20 ) & ( FILTER = 'PASS'  | GEN[*].AO > 20)" ${pair_id}.${method}.sv.norm.vcf.gz | java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} - | bgzip > ${pair_id}.${method}.vcf.gz
    
    if [[ $filter == 1 ]]
    then
	mv ${pair_id}.${method}.vcf.gz ${pair_id}.${method}.ori.vcf.gz
	perl $baseDir/filter_sv.pl -t $tid -p ${pair_id}.${method} -i ${pair_id}.${method}.ori.vcf.gz
	bgzip -f ${pair_id}.${method}.vcf
    fi
elif [[ $method == 'pindel' ]]
then
    genomefiledate=`find ${reffa} -maxdepth 0 -printf "%TY%Tm%Td\n"`
    bedopt=''
    if [[ -f $tbed ]]
    then
	bedopt="-j $tbed"
    fi
    pindel -T $NPROC -f ${reffa} -i ${pair_id}.pindel.config -o ${pair_id}.pindel_out --report_inversions false $bedopt
    pindel2vcf -P ${pair_id}.pindel_out -r ${reffa} -R HG38 -d ${genomefiledate} -v pindel.vcf
    cat pindel.vcf | java -jar $SNPEFF_HOME/SnpSift.jar filter "( GEN[*].AD[1] >= 10 )" | bgzip > pindel.vcf.gz
    tabix pindel.vcf.gz
    bash $baseDir/norm_annot.sh -r ${index_path} -p pindel -v pindel.vcf.gz
    perl $baseDir/parse_pindel.pl ${pair_id} pindel.norm.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.indel.vcf |bgzip > ${pair_id}.pindel_indel.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.dup.vcf | bedtools intersect -header -b ${itdbed} -a stdin | bgzip > ${pair_id}.pindel_tandemdup.vcf.gz
    java -Xmx10g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} ${pair_id}.sv.vcf | bgzip > ${pair_id}.pindel.sv.vcf.gz
    if [[ $filter == 1 ]]
    then
	perl $baseDir/filter_pindel.pl -d ${pair_id}.pindel_tandemdup.vcf.gz -s ${pair_id}.pindel.sv.vcf.gz -i ${pair_id}.pindel_indel.vcf.gz
	mv ${pair_id}.pindel_tandemdup.vcf.gz ${pair_id}.pindel_tandemdup.unfilt.vcf.gz
	mv ${pair_id}.pindel_tandemdup.pass.vcf ${pair_id}.pindel_tandemdup.vcf
	bgzip ${pair_id}.pindel_tandemdup.vcf
	mv ${pair_id}.pindel_indel.pass.vcf ${pair_id}.pindel.vcf
	bgzip ${pair_id}.pindel.vcf
	mv ${pair_id}.pindel.sv.vcf.gz ${pair_id}.pindel.sv.unfilt.vcf.gz
	mv ${pair_id}.pindel.sv.pass.vcf ${pair_id}.pindel.sv.vcf
	bgzip ${pair_id}.pindel.sv.vcf
	zgrep '#CHROM' ${pair_id}.pindel.sv.vcf.gz > ${pair_id}.pindel.genefusion.txt
	zcat ${pair_id}.pindel.sv.vcf.gz | $SNPEFF_HOME/scripts/vcfEffOnePerLine.pl |java -jar $SNPEFF_HOME/SnpSift.jar extractFields - CHROM POS CHROM END ANN[*].EFFECT ANN[*].GENE ANN[*].BIOTYPE  FILTER FORMAT GEN[*] |grep -E 'gene_fusion|feature_fusion' | sort -u >> ${pair_id}.pindel.genefusion.txt
    fi
elif [[ $method == 'itdseek' ]]
then
    stexe=`which samtools`
    echo $stexe
    samtools view -@ $NPROC -L ${itdbed} ${sbam} | itdseek.pl --refseq ${reffa} --samtools ${stexe} --bam ${sbam} | vcf-sort | bedtools intersect -header -b ${itdbed} -a stdin | java -Xmx30g -jar $SNPEFF_HOME/SnpSift.jar filter "( LEN < 10000 )" | bgzip > ${pair_id}.itdseek.vcf.gz
    
    tabix ${pair_id}.itdseek.vcf.gz
    
    bcftools norm --fasta-ref $reffa -m - -Ov ${pair_id}.itdseek.vcf.gz | java -Xmx30g -jar $SNPEFF_HOME/snpEff.jar -no-intergenic -lof -c $SNPEFF_HOME/snpEff.config ${snpeffgeno} - |bgzip > ${pair_id}.itdseek_tandemdup.vcf.gz
    if [[ $filter == 1 ]]
    then
	perl $baseDir/filter_itdseeker.pl -t ${pair_id} -d ${pair_id}.itdseek_tandemdup.vcf.gz
	mv ${pair_id}.itdseek_tandemdup.vcf.gz ${pair_id}.itdseek_tandemdup.unfilt.vcf.gz
	mv ${pair_id}.itdseek_tandemdup.pass.vcf ${pair_id}.itdseek_tandemdup.vcf
	bgzip ${pair_id}.itdseek_tandemdup.vcf
    fi
fi
