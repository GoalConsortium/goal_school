#!/cm/shared/apps/R/gcc/4.0.2/bin/Rscript

# module load python/3.7.x-anaconda R/4.0.2-gccmkl rstudio-desktop/1.1.456
# rstudio

#Set MiB limit needed for SCT Integration Step
Mb = 128000
options(future.globals.maxSize= Mb*1024^2)
options(java.parameters = "-Xmx128g")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
#loompy <- reticulate::import('loompy')
args <- commandArgs(trailingOnly=TRUE)
sample=args[1]

root='/project/BICF/BICF_Core/shared/bicf_helpdesk/issue605_Tabrizi/New_Data/analysis/cloud_10x'
setwd(root)

##############################

sObj <- readRDS(paste0(sample,"-annot.rds",sep=''))

convertFormat(sObj,from="seurat",to="anndata",outFile=paste0("/project/BICF/shared/cellxgene/DataShare/ReddyLab/HYFJMDSXY_20210401/",sample,".h5ad",sep=''),assay="SCT",main_layer="data")
