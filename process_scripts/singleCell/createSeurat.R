# module load python/3.7.x-anaconda R/4.0.2-gccmkl rstudio-desktop/1.1.456
# rstudio

#Set MiB limit needed for SCT Integration Step
Mb = 64000
options(future.globals.maxSize= Mb*1024^2)
options(java.parameters = "-Xmx64g")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
args <- commandArgs(trailingOnly=TRUE)

sample=args[1]

root='/project/BICF/BICF_Core/shared/bicf_helpdesk/issue605_Tabrizi/New_Data/analysis/cloud_10x'
setwd(root)

data <- Read10X(data.dir=paste0(sample,"_protgex","/filtered_feature_bc_matrix",sep=''))
sObj <- CreateSeuratObject(counts=data$`Gene Expression`, project=sample)
sObj[['Protein']] = CreateAssayObject(counts = data$`Antibody Capture`)

saveRDS(sObj,paste(sample,"rds",sep='.'))
