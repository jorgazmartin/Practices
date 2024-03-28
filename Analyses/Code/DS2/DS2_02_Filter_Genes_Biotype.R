# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_02_Filter_Genes_Biotype ****\n")
rm(list = ls())

###################
#### LOAD DATA ####
###################

reads <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/reads.txt")
lengths <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/lengths.txt")
tpms <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/tximport_tpms.txt")
genes_fly <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/feature_data.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D2/annotated/metadata.txt")

########################
#### BIOTYPE FILTER ####
########################

## Exploration: number of genes per chromosome.
cat("Genes per chromosome:\n")
summary(factor(genes_fly$chromosome_name))

## Exploration: number of genes by biotype.
cat("Genes by biotype:\n")
summary(factor(genes_fly$gene_biotype))

## Filter: We will focus our analysis in some selected chromosomes and only protein coding genes.
indexes_to_keep <- which((genes_fly$chromosome_name %in% c("2L","2R","3L","3R","4")) & (genes_fly$gene_biotype == "protein_coding"))
reads <- reads[indexes_to_keep,]
lengths <- lengths[indexes_to_keep,]
tpms <- tpms[indexes_to_keep,]
genes_fly <- genes_fly[indexes_to_keep,]

## Save results.
dir.create("Analyses/Inputs/Processed_datasets/D2/filtered_biotype")
write.table(genes_fly,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/feature_data.txt")
write.table(reads,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/reads.txt")
write.table(lengths,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/lengths.txt")
write.table(tpms,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/tximport_tpms.txt")
write.table(metadata,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/metadata.txt")

cat("**** END OF JOB DS2_02_Filter_Genes_Biotype ****\n")
