# module load python/3.6.4-anaconda R/4.0.2-gccmkl rstudio-desktop/1.1.456
# rstudio

#Set MiB limit needed for SCT Integration Step
Mb = 64000
options(future.globals.maxSize= Mb*1024^2)
options(java.parameters = "-Xmx64g")

library(Seurat)
library(sceasy)
library(reticulate)
use_condaenv('sceasy')
#py_install('anndata')
#loompy <- reticulate::import('loompy')
library(Matrix)
library(ggplot2)
library(Ecfun)
library(gridExtra)
library(viridis)
library(SingleR)
library(genefu)
library(tidyverse)
library(scater)
library(clustree)
library(celldex)
library(BiocParallel)

args <- commandArgs(trailingOnly=TRUE)

sample=args[1]
root='/project/BICF/BICF_Core/shared/bicf_helpdesk/issue605_Tabrizi/New_Data/analysis/cloud_10x'
setwd(root)

scObj <- readRDS(paste0(sample,"-filtered.rds",sep=''))

##############################
##### Annotate Clusters using automated naming with Human Primary Cell Atlas (HPCA)
#SingleR Annotations
SCE.scObj <- as.SingleCellExperiment(scObj)

### HumanPrimaryCellAtlasData
hpca.se <- HumanPrimaryCellAtlasData()
common_hpca.scObj <- intersect(rownames(SCE.scObj), rownames(hpca.se))

hpca.sescObj <- hpca.se[common_hpca.scObj,]
hpca.scObj <- SCE.scObj[common_hpca.scObj,]
hpca.scObj <- logNormCounts(hpca.scObj)

# Main Label Groups based on 0.1 Resolution Clusterng
pred.hpca.scObj <- SingleR(test = hpca.scObj,
                           ref = hpca.sescObj,
                           labels = hpca.sescObj$label.main,
                           method = "cluster",
                           clusters = hpca.scObj$SCT_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["HPCA_Main_Labels"]] <- pred.hpca.scObj$labels[scObj$SCT_snn_res.1]

# Fine Label Groups based on 0.1 Resolution Clustering
pred.hpca.scObj <- SingleR(test = hpca.scObj,
                           ref = hpca.sescObj,
                           labels = hpca.sescObj$label.fine,
                           method = "cluster",
                           clusters = hpca.scObj$SCT_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["HPCA_Fine_Labels"]] <- pred.hpca.scObj$labels[scObj$SCT_snn_res.1]

##############################
##### Annotate Clusters using automated naming with Monaco Immune Data Set (MIDS)
#MonacoImmuneDataSet
micd.se <- MonacoImmuneData()
common_micd.scObj <- intersect(rownames(SCE.scObj), rownames(micd.se))

micd.sescObj <- micd.se[common_micd.scObj,]
micd.scObj <- SCE.scObj[common_micd.scObj,]
micd.scObj <- logNormCounts(micd.scObj)

# Main Label Groups based on 1.0 Resolution Clusterng
pred.micd.scObj <- SingleR(test = micd.scObj,
                           ref = micd.sescObj,
                           labels = micd.sescObj$label.main,
                           method = "cluster",
                           clusters = micd.scObj$SCT_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["MICD_Main_Labels"]] <- pred.micd.scObj$labels[scObj$SCT_snn_res.1]


# Fine Label Groups based on 1.0 Resolution Clusterng
pred.micd.scObj <- SingleR(test = micd.scObj,
                           ref = micd.sescObj,
                           labels = micd.sescObj$label.fine,
                           method = "cluster",
                           clusters = micd.scObj$SCT_snn_res.1,
                           BPPARAM = MulticoreParam(workers=20))
scObj[["MICD_Fine_Labels"]] <- pred.micd.scObj$labels[scObj$SCT_snn_res.1]
saveRDS(scObj, file=paste0(sample,"-annot.rds",sep=''))
