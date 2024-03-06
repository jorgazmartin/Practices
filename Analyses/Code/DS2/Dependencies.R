library(renv)

library(tidyverse)
library(readxl)
library(writexl)
library(reshape2)
library(ggrepel)
library(cowplot)
library(rgl)


## Genomic annotations:
library(biomaRt)

## Unsupervised methods for dimensional reduction & clustering
library(Rtsne)
library(umap)
library(dendextend)

## RNA-seq quantification & Differential expression libraries.
library(tximport)
library(edgeR)
library(limma)
library(DESeq2)

library(mashr)

## Useful libraries for running enrichment analyses and interpreting stats results
library(fgsea)
library(clusterProfiler)

library(org.Hs.eg.db)
library(AnnotationDbi)
library(msigdbr)
library(enrichplot)
library(ggnewscale)

###############################
#### ENVIRONMENT VARIABLES ####
###############################

rootDir <- "C:/Users/jorge/OneDrive/Master BQB/S2_Big Data/Practices/Analysis/DS2/"
inputDir <- paste0(rootDir,"1.Input/")
outputDir <- paste0(rootDir,"3.Output/")
