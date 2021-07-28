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
library(DoubletFinder)


args <- commandArgs(trailingOnly=TRUE)

sample=args[1]
root='/project/BICF/BICF_Core/shared/bicf_helpdesk/issue605_Tabrizi/New_Data/analysis/cloud_10x'
setwd(root)
tcellloc <- paste0(sample,"_T/",sep='')
bcellloc <- paste0(sample,"_B/",sep='')

##############################

add_clonetype <- function(tcr_folder, seurat_obj){
    tcr <- read.csv(paste(tcr_folder,"filtered_contig_annotations.csv", sep=""))
    # Remove the -1 at the end of each barcode.
    # Subsets so only the first line of each barcode is kept,
    # as each entry for given barcode will have same clonotype.
    tcr$barcode <- gsub("-1", "", tcr$barcode)
    tcr <- tcr[!duplicated(tcr$barcode), ]
    # Only keep the barcode and clonotype columns. 
    # We'll get additional clonotype info from the clonotype table.
    tcr <- tcr[,c("barcode", "raw_clonotype_id")]
    names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
    # Clonotype-centric info.
    clono <- read.csv(paste(tcr_folder,"clonotypes.csv", sep=""))
    # Slap the AA sequences onto our original table by clonotype_id.
    tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")])
    # Reorder so barcodes are first column and set them as rownames.
    tcr <- tcr[, c(2,1,3)]
    rownames(tcr) <- tcr[,1]
    tcr[,1] <- NULL
    # Add to the Seurat object's metadata.
    clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
    return(clono_seurat)
    }

sObj <- readRDS(paste(sample,"rds",sep='.'))
sObj <- PercentageFeatureSet(sObj, pattern = "^MT-", col.name = "percent.mito")

sObj <- SCTransform(sObj, method = "glmGamPoi", vars.to.regress = "percent.mito", verbose = FALSE, return.only.var.genes = FALSE, assay = "RNA")
sObj <- RunPCA(sObj, features = VariableFeatures(object = sObj), assay = "SCT", verbose = FALSE)

########## Run clustering dim of pc.use and resolution 0.1
res <- c(0.1,0.5,1.0,1.5,2.0,2.5,3.5)
sObj <- FindNeighbors(object = sObj, dims = 1:50, assay = "SCT")
sObj <- FindClusters(object = sObj, resolution = res, assay = "SCT")

### TSNE and UMAP, add PCA
sObj <- RunTSNE(sObj,dims=1:50,check_duplicates=FALSE, assay = "SCT")
sObj <- RunUMAP(sObj,dims=1:50, assay="SCT")

### Add cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
sObj <- CellCycleScoring(sObj, s.features=s.genes, g2m.features=g2m.genes, set.ident=T, assay = "SCT")

### Add B-Cell and T-Cell Annotation

sObj <- add_clonetype(tcellloc, sObj)
sObj <- add_clonetype(bcellloc, sObj)

### Filter QC ###

numcells <- nrow(sObj@meta.data)
expdouble <- 0.008*(numcells/1000)
nExp_poi <- round(expdouble*numcells)

if (nExp_poi > 20) {
   sObj <- doubletFinder_v3(sObj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = TRUE)
   pann.column <- str_subset(names(sObj@meta.data),'pANN')
   homotypic.prop <- modelHomotypic(sObj@meta.data$seurat_clusters)
   nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
   sObj <- doubletFinder_v3(sObj, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = pann.column, sct = FALSE)
   dfclass.name <- str_subset(names(sObj@meta.data),'DF.classifications')[2]
   sObj$doubletClass <- sObj[[dfclass.name]]
   sObj <- subset(sObj, subset = nFeature_RNA > 200 & doubletClass == 'Singlet' & percent.mito < 15)
} else {
   sObj <- subset(sObj, subset = nFeature_RNA > 200 & percent.mito < 15)
}

saveRDS(sObj, file=paste0(sample,"-filtered.rds",sep=''))
