source("Analyses/Code/Dependencies.R")

## Lengths, reads, tpms, feature_data and meta_data from the annotated folder:
reads=read.table("Analyses/Inputs/Processed_datasets/D2/annotated/reads.txt")
lengths=read.table("Analyses/Inputs/Processed_datasets/D2/annotated/lengths.txt")
tpms=read.table("Analyses/Inputs/Processed_datasets/D2/annotated/tpms.txt")
genes_fly=read.table("Analyses/Inputs/Processed_datasets/D2/annotated/feature_data.txt")
metadata=read.table("Analyses/Inputs/Processed_datasets/D2/annotated/metadata.txt")


## let us check the dimensions of these objects and their conformation:
dim(reads)
dim(lengths)
dim(tpms)
dim(metadata)
## compare to: 
dim(genes_fly)
## And just to be extra sure:
length(which(rownames(reads)==rownames(lengths)))
length(which(rownames(reads)==rownames(tpms)))
length(which(rownames(reads)==rownames(genes_fly)))
## As well as:
length(which(colnames(reads)==colnames(tpms)))
length(which(colnames(reads)==colnames(lengths)))
length(which(colnames(reads)==rownames(metadata)))

## first let us check things are aligned:
unique(genes_fly$chromosome_name)
summary(factor(genes_fly$chromosome_name))

## first let us check things are aligned:
unique(genes_fly$gene_biotype)
summary(factor(genes_fly$gene_biotype))

indexes_to_keep=which((genes_fly$chromosome_name %in% c("2L","2R","3L","3R","4")) & (genes_fly$gene_biotype == "protein_coding"))
reads=reads[indexes_to_keep,]
lengths=lengths[indexes_to_keep,]
tpms=tpms[indexes_to_keep,]
genes_fly=genes_fly[indexes_to_keep,]