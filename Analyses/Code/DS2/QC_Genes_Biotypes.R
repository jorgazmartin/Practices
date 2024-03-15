# Load dependencies.
source("Analyses/Code/Dependencies.R")

Validation <- function (lengths, metadata, reads, tpms) {
  
  fail <- TRUE
  
  # Check that all the objects have the same number of genes (rows)
  if (dim(tpms)[1] != dim(reads)[1] || dim(lengths)[1] != dim(reads)[1]) {
    fail <- FALSE
    cat("VALIDATION | NOK - Tables do not have the same number of genes\n",
        "Featured genes: ",dim(genes_fly)[1],"\n",
        "Reads: ",dim(reads)[1],"\n",
        "Lengths: ",dim(lengths)[1],"\n",
        "TPMs: ",dim(tpms)[1],"\n",sep="")
  }
  
  ## VALIDATION: Check that the orders match (we know there are the same number, but are they the same?)
  if (length(which(rownames(reads) != rownames(lengths))) != 0) {
    fail <- FALSE
    cat("VALIDATION | NOK - Reads and lengths are not aligned\n")
  }
  if (length(which(rownames(reads) != rownames(tpms))) != 0) {
    fail <- FALSE
    cat("VALIDATION | NOK - Reads and tpms are not aligned\n")
  }
  
  ## VALIDATION: Meta-data Check order between samples (column-wise) and meta_data annotations (row-wise)
  if (length(which(rownames(metadata) != colnames(reads))) != 0){
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched metadata in reads table\n")
  }
  if (length(which(rownames(metadata) != colnames(lengths))) != 0){
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched metadata in lengths table\n")
  }
  if (length(which(rownames(metadata) != colnames(tpms))) != 0){
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched metadata in lengths table\n")
  }
  
  return (fail)
  
}

## Load Lengths, reads, tpms, feature_data and meta_data from the annotated folder:
reads <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/reads.txt")
lengths <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/lengths.txt")
tpms <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/tpms.txt")
genes_fly <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/feature_data.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/metadata.txt")

# Validate that the input data is OK
if(Validation(lengths, metadata, reads, tpms)) {
  
  ## Exploration: number of genes per chromosome.
  unique(genes_fly$chromosome_name)
  summary(factor(genes_fly$chromosome_name))
  
  ## Exploration: number of genes by biotype.
  unique(genes_fly$gene_biotype)
  summary(factor(genes_fly$gene_biotype))
  
  ## Filter: We will focus our analysis in some selected chromosomes and only protein coding genes.
  indexes_to_keep=which((genes_fly$chromosome_name %in% c("2L","2R","3L","3R","4")) & (genes_fly$gene_biotype == "protein_coding"))
  reads=reads[indexes_to_keep,]
  lengths=lengths[indexes_to_keep,]
  tpms=tpms[indexes_to_keep,]
  genes_fly=genes_fly[indexes_to_keep,]
  
}

save(reads,lengths,tpms,genes_fly,metadata,file="Analyses/Code/DS2/filteredBiotypes.RData")
