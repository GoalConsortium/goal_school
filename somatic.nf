#!/usr/bin/env nextflow

// Default parameter values to run tests
params.input = './'
params.output = './'

params.bams="$params.input/*.bam"
params.design="$params.input/design.txt"
params.callsvs="skip"
params.genome="/project/shared/bicf_workflow_ref/GRCh38"
params.targetpanel="$params.genome/UTSWV2.bed"

dbsnp="$params.genome/dbSnp.vcf.gz"
cosmic="$params.genome/cosmic.vcf.gz"

design_file = file(params.design)
bams=file(params.bams)

index_path = file(params.genome)
capture_bed = file(params.targetpanel)
ncmconf = file("$params.genome/ncm.conf")
dbsnp=file(dbsnp)

snpeff_vers = 'GRCh38.86';

def fileMap = [:]

bams.each {
    final fileName = it.getFileName().toString()
    prefix = fileName.lastIndexOf('/')
    fileMap[fileName] = it
}
def oribam = []
def tarbam = []
new File(params.design).withReader { reader ->
    def hline = reader.readLine()
    def header = hline.split("\t")
    tidx = header.findIndexOf{it == 'TumorID'};
    nidx = header.findIndexOf{it == 'NormalID'};
    oneidx = header.findIndexOf{it == 'TumorBAM'};
    twoidx = header.findIndexOf{it == 'NormalBAM'};
    totidx = header.findIndexOf{it == 'TumorFinalBAM'};
    notidx = header.findIndexOf{it == 'NormalFinalBAM'};

    if (twoidx == -1) {
       twoidx = oneidx
       }      
    while (line = reader.readLine()) {
    	   def row = line.split("\t")
	   if (fileMap.get(row[oneidx]) != null) {
	      oribam << tuple(row[tidx],row[nidx],fileMap.get(row[oneidx]),fileMap.get(row[twoidx]))
	      tarbam << tuple(row[tidx],row[nidx],fileMap.get(row[totidx]),fileMap.get(row[notidx]))
	   }
} 
}

if( ! oribam) { error "Didn't match any input files with entries in the design file" }
if( ! tarbam) { error "Didn't match any input files with entries in the design file" }

process indexoribams {
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal) from oribam
  output:
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into dellybam
  set val("${tid}_${nid}"),tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mantrabam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into checkbams

  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}
process indextarbams {
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal) from tarbam
  output:
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into mutectbam
  set val("${tid}_${nid}"),tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into strelkabam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into ssbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into shimmerbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into vscanbam
  set tid,nid,file(tumor),file(normal),file("${tumor}.bai"),file("${normal}.bai") into virmidbam
  set val("${tid}_${nid}"),tid,nid into pairnames
  script:
  """
  bash $baseDir/process_scripts/alignment/indexbams.sh 
  """
}
process checkmates {
  publishDir "$params.output", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from checkbams
  file(conf) from ncmconf
  output:
  file("${tid}_${nid}*") into checkmateout
  script:
  """
  source /etc/profile.d/modules.sh
  module load python/2.7.x-anaconda
  python /project/shared/bicf_workflow_ref/seqprg/NGSCheckMate/ncm.py -B -d ./ -bed ${index_path}/NGSCheckMate.bed -O ./ -N ${tid}_${nid}
  perl $baseDir/scripts/sequenceqc_somatic.pl -r ${index_path} -i ${tid}_${nid}_all.txt -o ${tid}_${nid}.sequence.stats.txt
  """
}

process svcall {
  publishDir "$params.output", mode: 'copy'
  errorStrategy 'ignore'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from dellybam

  output:
  file("${tid}_${nid}.delly.vcf.gz") into dellyvcf
  file("${tid}_${nid}.sssv.sv.vcf.gz") into lumpyvcf
  file("${tid}_${nid}.sv.annot.vcf.gz") into svintvcf
  file("${tid}_${nid}.sv.annot.txt") into svannot
  file("${tid}_${nid}.sv.annot.genefusion.txt") into gfannot
  when:
  params.callsvs == "detect"
  script:
  """
  source /etc/profile.d/modules.sh
  perl $baseDir/scripts/make_delly_sample.pl ${tid} ${nid}
  bash $baseDir/process_scripts/variants/svcalling.sh -r ${index_path} -p ${tid}_${nid} -b ${tumor} -n ${normal} -k ${tid}
  """
}
process sstumor {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from ssbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.sssom.vcf.gz") into ssvcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a speedseq
  """
}
process mutect {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'

  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from mutectbam

  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.mutect.vcf.gz") into mutectvcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a mutect2
  """
}
Channel
  .empty()
  .mix(mantrabam,strelkabam)
  .groupTuple(by:0)
  .into { illuminabams }

process strelka {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'

  input:
  set pair_id,mtid,mnid,file(mtumor),file(mnormal),file(mtidx),file(mnidx),tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from illuminabams

  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.strelka2.vcf.gz") into strelkavcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a strelka2
  """
}
process varscan {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from vscanbam
  output:
  set val("${tid}_${nid}"),file("${tid}_${nid}.varscan.vcf.gz") into varscanvcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a varscan
  """
}

process shimmer {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from shimmerbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.shimmer.vcf.gz") into shimmervcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a shimmer
  """
}

process virmid {
  errorStrategy 'ignore'
  //publishDir "$params.output", mode: 'copy'
  input:
  set tid,nid,file(tumor),file(normal),file(tidx),file(nidx) from virmidbam
  output:
  set val("${tid}_${nid}"), file("${tid}_${nid}.virmid.vcf.gz") into virmidvcf
  script:
  """
  bash $baseDir/process_scripts/variants/somatic_vc.sh -r $index_path -x $tid -y $nid -n $normal -t $tumor -a virmid
  """
}

Channel
  .empty()
  .mix(ssvcf,mutectvcf,varscanvcf,virmidvcf,shimmervcf)
  .groupTuple(by:0)
  .into { vcflist}


process integrate {
  errorStrategy 'ignore'
  publishDir "$params.output", mode: 'copy'
  input:
  set subjid,file(vcf) from vcflist
  output:
  set subjid,file("${subjid}.union.vcf.gz") into union
  file("${subjid}.annot.vcf.gz") into annotvcf
  script:
  """
  bash $baseDir/process_scripts/variants/union.sh -r $index_path -p $subjid
  bash $baseDir/process_scripts/variants/annotvcf.sh -p $subjid -r $index_path -v ${subjid}.union.vcf.gz
  """
}
