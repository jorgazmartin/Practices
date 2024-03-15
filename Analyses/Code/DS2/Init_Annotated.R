# Load the dependencies.
source("Analyses/Code/Dependencies.R")

Validation <- function (genes_fly, lengths, metadata, reads, tpms) {
  
  fail <- TRUE
  
  ## VALIDATION: Is the mapping leading to the feature data annotations unequivocal?
  if (dim(genes_fly)[1] != length(unique(genes_fly$flybase_gene_id))) {
    fail <- FALSE
    cat("VALIDATION | NOK - Mapping leadin to feature data annotations not unequivocal\n")
  }
  
  if (dim(reads)[1] != dim(genes_fly)[1] || dim(tpms)[1] != dim(genes_fly)[1] || dim(lengths)[1] != dim(genes_fly)[1]) {
    fail <- FALSE
    cat("VALIDATION | NOK - Tables do not have the same number of genes\n",
        "Featured genes: ",dim(genes_fly)[1],"\n",
        "Reads: ",dim(reads)[1],"\n",
        "Lengths: ",dim(lengths)[1],"\n",
        "TPMs: ",dim(tpms)[1],"\n",sep="")
  }
  
  ## VALIDATION: Check that the orders match (we know there are the same number, but are they the same?)
  if (length(which(rownames(reads) != rownames(genes_fly))) != 0) {
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched genes in reads table\n")
  }
  if (length(which(rownames(lengths) != rownames(genes_fly))) != 0) {
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched genes in lengths table\n")
  }
  if (length(which(rownames(tpms) != rownames(genes_fly))) != 0) {
    fail <- FALSE
    cat("VALIDATION | NOK - Unmatched genes in tpms table\n")
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

# Load the tables (reads, lengths, tpms, metadata).
reads <- read.table("Analyses/Inputs/Raw_datasets/D2/reads.txt")
lengths <- read.table("Analyses/Inputs/Raw_datasets/D2/lengths.txt")
tpms <- read.table("Analyses/Inputs/Raw_datasets/D2/tpms.txt")
metadata <- read.table("Analyses/Inputs/Raw_datasets/D2/meta_data.txt")

# Query ensemble to load the mart or load the local data.
first_time <- TRUE
if(first_time)
{
  fly <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="dmelanogaster_gene_ensembl")
  dir.create("Analyses/Inputs/Raw_datasets/D2/marts",recursive=TRUE)
  save(fly,file="Analyses/Inputs/Raw_datasets/D2/marts/fly_mart.Rdata")
} else {
  base::load(file="Analyses/Inputs/Raw_datasets/D2/marts/fly_mart.Rdata")
}

# possible_attributes <- listAttributes(fly)
# head(possible_attributes)

# Get the desired gene attributes only for the genes in reads.
genes_fly <- getBM(attributes = c("flybase_gene_id","start_position","end_position","chromosome_name","gene_biotype"), filters = "flybase_gene_id", values = rownames(reads) , mart = fly)
rownames(genes_fly) <- genes_fly$flybase_gene_id

## How many of our genes we caught?
num_unfeatured_genes <- length(which(!(rownames(reads) %in% rownames(genes_fly))))
cat("Number of unfeatured genes: ",num_unfeatured_genes,"\n",sep="")

# Filter reads, tpms and lengths: we just take the featured genes.
reads <- reads[which(rownames(reads) %in% rownames(genes_fly)),]
lengths <- lengths[which(rownames(lengths) %in%rownames(reads)),]
tpms <- tpms[which(rownames(tpms) %in% rownames(reads)),]

# Order by gene (alphabetically)
genes_fly <- genes_fly[order(rownames(genes_fly)),]
reads <- reads[order(rownames(reads)),]
tpms <- tpms[order(rownames(tpms)),]
lengths <- lengths[order(rownames(lengths)),]

if (Validation(genes_fly, lengths, metadata, reads, tpms)) {
  ## SAVE: Now save the reads table as well as the feature table after removing those 43 un-annotated genes.
  dir.create("Analyses/Inputs/Processed_datasets/D2/annotated/",recursive=TRUE)
  write.table(reads,"Analyses/Inputs/Processed_datasets/D2/annotated/reads.txt")
  write.table(genes_fly,"Analyses/Inputs/Processed_datasets/D2/annotated/feature_data.txt")
  write.table(lengths,"Analyses/Inputs/Processed_datasets/D2/annotated/lengths.txt")
  write.table(tpms,"Analyses/Inputs/Processed_datasets/D2/annotated/tpms.txt")
  write.table(metadata,"Analyses/Inputs/Processed_datasets/D2/annotated/metadata.txt")
}


