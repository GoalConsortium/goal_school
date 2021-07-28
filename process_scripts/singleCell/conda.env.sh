#!/bin/bash

module load python/3.7.x-anaconda R/4.0.2-gccmkl rstudio-desktop/1.1.456

conda activate sceasy

conda install -c bioconda r-sceasy
conda install anndata==0.6.19 scipy==1.2.1 -c bioconda
conda install loompy -c bioconda
