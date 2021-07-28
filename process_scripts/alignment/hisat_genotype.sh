#!/bin/bash
#hlatyping.sh

usage() {
  echo "-h Help documentation for dnaseqalign.sh"
  echo "-r  --Reference Genome: GRCh38 or GRCm38"
  echo "-x  --FastQ R1"
  echo "-y  --FastQ R2"
  echo "-p  --Prefix for output file name"
  echo "Example: bash hisat_genotype.sh -p prefix"
  exit 1
}
OPTIND=1 # Reset OPTIND
while getopts :x:y:p:uh opt
do
    case $opt in
        x) fq1=$OPTARG;;
        y) fq2=$OPTARG;;
	u) umi='umi';;
        p) pair_id=$OPTARG;;
        h) usage;;
    esac
done

shift $(($OPTIND -1))

# Check for mandatory options
if [[ -z $pair_id ]]; then
    usage
fi

NPROC=$SLURM_CPUS_ON_NODE
if [[ -z $NPROC ]]
then
    NPROC=`nproc`
fi

if [[ -z $isdocker ]]
then
    source /etc/profile.d/modules.sh
    module load hisat-genotype/1.0.1
fi
diff $fq1 $fq2 > difffile

if [ -s difffile ]
then
hisatgenotype_extract_reads.py  -p 16  --database-list hla --base /project/shared/bicf_workflow_ref/human/hisat_genotype_hla/genotype_genome -1 $fq1 -2 $fq2 --out-dir hisatgeno_out
else
hisatgenotype_extract_reads.py  -p 16  --database-list hla --base /project/shared/bicf_workflow_ref/human/hisat_genotype_hla/genotype_genome -U $fq1 --out-dir hisatgeno_out
fi
hisatgenotype_locus_samples.py -p 16 --region-list hla.A,hla.B,hla.C,hla.DQA1,hla.DQB1,hla.DRB1,hla.DPB1 --assembly --read-dir hisatgeno_out --out-dir ${pair_id}.hisat_hla > ${pair_id}.hisat_hla.txt

tar cf ${pair_id}.hisat_hla.tar ${pair_id}.hisat_hla
gzip ${pair_id}.hisat_hla.tar
