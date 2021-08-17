# Software for Clinical Health Omics Oncology Laboratories
School is a collection of genomics analysis workflows that are used for detecting single nucleotide variants (SNVs), insertions/deletions (indels), copy number variants (CNVs) and translocations from RNA and DNA sequencing. These workflows were first developed and validated in a CLIA laboratory at UTSW, and will continue to be developed and maintained by the Genomics Organization for Academic Laboratories (GOAL) Consortium.

## Prerequisites

These bioinformatics pipelines use **`nextflow`**, a framework for defining and executing a directed acyclic graph (DAG) of interdependent steps. Also required is either **`singularity`** or **`docker`**, executors for tools that are [containerized](https://www.docker.com/resources/what-container) for portability across computing environments.

[Follow these instructions](https://gist.github.com/ckandoth/982ce140b4dd9d6bf72a780c05a549a3) to install Nextflow and Singularity in a Linux environment. For better portability across computing environments (Linux, macOS, Windows), [follow these instructions](https://docs.docker.com/get-docker) to install Docker. Docker requires administrative rights, which you normally have on a personal laptop/workstation. But in shared computers like HPC clusters, there are [valid concerns](https://duo.com/decipher/docker-bug-allows-root-access-to-host-file-system) against installing Docker, and then Singularity makes more sense.

Some HPC clusters will have these tools pre-installed as [environment modules](https://modules.readthedocs.io/en/latest/). Use command `module avail` to see what's available, and `module load` to load them into your `$PATH`. But make sure you have Nextflow 20.07.1 or newer and Singularity 3.5.2 or newer.

## Quick Start

Clone a branch of this repo that you want to test and `cd` into it:
```bash
git clone -b UFHPL_branch_1 --single-branch https://github.com/goalconsortium/goal_school.git
cd goal_school
```

Download and unzip resource files needed by the pipeline (39GB download that unzips to 43GB):
```bash
curl -LO https://reference-files-bucket.s3.amazonaws.com/Reference_Files.zip
unzip Reference_Files.zip
```

In a folder named `fastq`, download small FASTQs created from DNA-seq of a tumor (or use your own FASTQs):
```bash
mkdir fastq
wget -P fastq https://github.com/mskcc/roslin-variant/raw/2.4.x/setup/examples/data/fastq/DU874145-T/DU874145-T_IGO_00000_TEST_L001_R{1,2}_001.fastq.gz
```

Prepare a design file for the DNA-seq variant calling pipeline:
```bash
echo -e "SampleID\tCaseID\tTumorID\tNormalID\tFqR1\tFqR2" > fastq/design.txt
echo -e "DU874145-T\tDU874145\tDU874145-T\t\tDU874145-T_IGO_00000_TEST_L001_R1_001.fastq.gz\tDU874145-T_IGO_00000_TEST_L001_R2_001.fastq.gz" >> fastq/design.txt
```

Run the `goalConsensus.nf` pipeline using the `standard` profile that uses Nextflow with singularity on the local machine:
```bash
nextflow run -work-dir .nextflow_work -profile standard goalConsensus.nf --input fastq --output analysis --repoDir ${PWD} --seqrunid H7YRLADXX --genome Reference_Files
```

To run this workflow on a different computing environment, lookup the institute-specific profiles in `nextflow.config` and/or create your own.

# Run Nextflow Workflows

This workflow can run either DNA or RNA sequencing. Please determine the desired configuration to achieve the proper analysis run.

## SlideRule: DNA Workflow

### DNA Design File

The design file must named design.txt and be in tab seperated format for the workflows. This workflow can be run with tumor-only or with tumor and normal pairs. If running tumor-only then do not include the NormalID collumn in the design file.

| SampleID | CaseID | TumorID | NormalID | FqR1 | FqR2 |
|---|---|---|---|---|---|
| Sample1 | Fam1 | Sample1 | Sample2 | Sample1.R1.fastq.gz | Sample1.R2.fastq.gz |
| Sample2 | Fam1 | Sample1 | Sample2 | Sample2.R1.fastq.gz | Sample2.R2.fastq.gz |
| Sample3 | Fam2 | Sample3 | Sample4 | Sample3.R1.fastq.gz | Sample3.R2.fastq.gz |
| Sample4 | Fam2 | Sample3 | Sample4 | Sample4.R1.fastq.gz | Sample4.R2.fastq.gz |

### DNA Parameters
* **--input**
  * directory containing the design file and fastq files
  * default is set to *'${basedir}/fastq'*
  * eg: **--input '/project/shared/bicf_workflow_ref/workflow_testdata/germline_variants/fastq'**
* **--output**
  * directory for the analysis output
  * default is set to *'${basedir}/analysis'*
  * eg: **--output '${basedir}/output'**
* **--seqrunid**
  * Illumina Run ID, used to distinguish repetitive runs
  * default set  to *'runtest'*
  * eg: **--seqrunid 'run'**
* **--genome**
  * directory containing all reference files for the various tools. This includes the genome.fa, dbsnp.vcf.gz, GoldIndels.vcf.gz, ncm.conf, ect.
  * default is set for use on UTSW BioHPC.
  * eg: **--genome '/project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref'**
* **--virus_genome**
  * directory containing viral genome reference files
  * default is set for use on UTSW BioHPC
  * eg: **--viral_genome '/project/shared/bicf_workflow_ref/human_virus_genome/clinlab_idt_genomes'**
* **--pon**
  * pon file for mutect capture kit bed file
  * default is set to not included. If a pon file is included, then mutect will run with reference to the inputted file
  * eg: **--pon '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/mutect2.pon.vcf.gz'** 
* **--capture**
  * bed file containing information on gene capture kit
  * eg: **--capture '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/targetpanel.bed'**
* **--capturedir**
  * directory containing capture bed and supporting files
  * eg: **--capturedir '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme'**
* **--platypus**
  * select whether to *'skip'* platypus algorithm or *'detect'* using platypus
  * default is set to *'skip'*
  * eg: **--platypus 'detect'**
* **--markdups**
  * select the method for marking duplicates from *'picard'*, *'samtools'*, *'fgbio_umi'*, *'picard_umi'*, or *'none'*
  * default is set to *'fgbio_umi'*
  * eg: **--markdups 'picard_umi'**
* **--version**
  * version of workflow analysis pipeline
  * default is set to *'v4'*
  * For current git version run *'gittag=$(git describe --abbrev=0 --tags)'*
  * eg: **--version $gittag**
* **--snpeff_vers**
  * version of reference genome for snpeff tool
  * default is set to *'GRCh38.86'*
  * eg: **--snpeff_vers 'GRCh38.86'**

### DNA Run Workflow

```
nextflow run -w $workdir ${baseDir}/dna.nf --input /project/shared/bicf_workflow_ref/workflow_testdata/germline_variants/fastq --output ${basedir}/output --seqrunid 'SHI1333-27' --pon /project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/mutect2.pon.vcf.gz --capture /project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/targetpanel.bed --capturedir /project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme --version $gittag --genome /project/shared/bicf_workflow_ref/human/grch38_cloud/dnaref -resume
```

## Abbacus: RNASeq Workflow

The RNA workflow can be run in with the whole genome, or with a specific list of genes of interest.

### RNA Design File

The design file must named design.txt and be in tab seperated format for the workflows. All RNA workflows can be run usin the same design file format.

| SampleID | CaseID | FqR1 | FqR2 |
|---|---|---|---|
| Sample1 | Fam1 | Sample1.R1.fastq.gz | Sample1.R2.fastq.gz |
| Sample2 | Fam1 | Sample2.R1.fastq.gz | Sample2.R2.fastq.gz |
| Sample3 | Fam2 | Sample3.R1.fastq.gz | Sample3.R2.fastq.gz |
| Sample4 | Fam2 | Sample4.R1.fastq.gz | Sample4.R2.fastq.gz |

### RNA Parameters
* **--input**
  * directory containing the design file and fastq files
  * default is set to *'${basedir}/fastq'*
  * eg: **--input '/project/shared/bicf_workflow_ref/workflow_testdata/rnaseq/fastq'**
* **--output**
  * directory for the analysis output
  * default is set to *'${basedir}/analysis'*
  * eg: **--output '${basedir}/output'**
* **--seqrunid**
  * Illumina Run ID, used to distrinquished repetative runs
  * default set  to *'runtest'*
  * eg: **--seqrunid 'run'**
* **--genome**
  * directory containing all reference files for the various tools. This includes the genome.fa, gencode.gtf, genenames.txt, ect.
  * default is set for use on UTSW BioHPC.
  * eg: **--genome '/project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref'**
* **--geneinfo**
  * path to textfile with the geneinfo for the human genome with reference to the genome being used
  * default is set for use on UTSW BioHPC
  * eg: **--geneinfo '/project/shared/bicf_workflow_ref/human/gene_info.human.txt'**
* **--version**
  * version of workflow analysis pipeline
  * default is set to *'v5'*
  * For current git version run *'gittag=$(git describe --abbrev=0 --tags)'*
  * eg: **--version $gittag**
* **--snpeff_vers**
  * version of reference genome for snpeff tool
  * default is set to *'GRCh38.86'*
  * eg: **--snpeff_vers 'GRCh38.86'**
* **--stranded**
  * option for -s flag in featurecount used in geneabundance calculations
  * default is set to *'0'*
  * eg: **--stranded '0'**
* **--pairs**
  * select either 'pe' (paired-end) or 'se' (single-end) based on read inputs. Select 'pe' when both R1 and R2 are present. If only R1, then select 'se'.
  * default is set to *'pe'*
  * eg: **--pairs 'pe'**
* **--align**
  * select the algorithm/tool for alignment from 'hisat' or 'star'
  * default is set to *'hisat'*
  * eg: **--align 'hisat'**
* **--bamct**
  * choose to either 'detect' or 'skip' bam read count process
  * default is set to *'detect'*
  * eg: **--bamct 'detect'**
* **--glist**
  * this is an optional parameter. If a genelist is included, then the whole rnaseq is not run but instead is restricted to the list containing genes of interest.
  * default is set to run whole rnaseq
  * eg: **--glist '/project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/genelist.txt'**
* **--umi**
  * This is an optional paramter. If included, then the -u option is added then UMI sequences are in the fastq read names.
  * default is set to not include umi
  * eg: **--umi 'true'**

### RNA Run Workflow

Whole rnaseq example

```
nextflow run -w $workdir ${baseDir}/rna.nf --input /project/shared/bicf_workflow_ref/workflow_testdata/rnaseq/fastq --output ${basedir}/output --seqrunid 'test' --version $gittag --genome /project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref --geneinfo /project/shared/bicf_workflow_ref/human/gene_info.human.txt -resume
```

GOI rnaseq example

```
nextflow run -w $workdir ${baseDir}/rna.nf --input /project/shared/bicf_workflow_ref/workflow_testdata/rnaseq/fastq --output ${basedir}/output --seqrunid 'test' --version $gittag --genome /project/shared/bicf_workflow_ref/human/grch38_cloud/rnaref --geneinfo /project/shared/bicf_workflow_ref/human/gene_info.human.txt --glist /project/shared/bicf_workflow_ref/human/grch38_cloud/panels/UTSW_V4_heme/genelist.txt -resume
```
