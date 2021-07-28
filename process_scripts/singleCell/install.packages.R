install.packages('Seurat')
install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("LoomExperiment","SingleR","genefu"))

remotes::install_github("cellgeni/sceasy",force=TRUE)
install.packages(c("reticulate","Matrix","ggplot2","Ecfun",
                   "gridExtra","viridis","tidyverse"))
install.packages(c("NMI","scater","clustree",'xlsx'))
BiocManager::install(c("celldex","BiocParallel","scater"))
