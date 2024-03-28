# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_04_Filter_Genes_Expression ****\n")
rm(list = ls())

###################
#### FUNCTIONS ####
###################

## Filter with tentative threshold: at least median tpm>=1 in at least 1 group:
filter <- function(tpms,threshold,metadata,grouping_factor){
  number_of_groups <- length(unique(metadata[,grouping_factor]))
  # Create medians matrix as an empty container.
  medians <- matrix(0,nrow=dim(tpms)[1],ncol=number_of_groups)
  colnames(medians) <- unique(metadata[,grouping_factor])
  rownames(medians) <- rownames(tpms)
  # Fill medians data.
  for(i in 1:number_of_groups)
  {
    set=which(metadata[,grouping_factor]==colnames(medians)[i])
    chunk <- tpms[,set]
    medians[,i] <- apply(chunk,1,median)
  }
  # We only keep those whose maximum median is bigger than the threshold.
  medians <- cbind(medians,max_median=apply(medians,1,max))
  genes_to_keep=rownames(medians)[which(medians[,"max_median"]>threshold)]
  
  return(genes_to_keep)
}

###################
#### LOAD DATA ####
###################

reads <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/reads.txt")
lengths <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/lengths.txt")
tpms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/tximport_tpms.txt")
genes_fly <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/feature_data.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/metadata.txt")
tpms_1 <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/built_tpms.txt")
fpkms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/fpkms.txt")

################################
#### GENE EXPRESSION FILTER ####
################################

# Filter: with the criteria of the function defined above.
genes_to_keep <- filter(tpms=tpms,threshold=1,metadata=metadata,grouping_factor="Condition")
genes_fly <- genes_fly[genes_to_keep,]
filtered_reads <- reads[genes_to_keep,]
filtered_lengths <- lengths[genes_to_keep,]
filtered_tpms <- tpms[genes_to_keep,]
filtered_tpms_1 <- tpms_1[genes_to_keep,]

## Save results.
dir.create("Analyses/Inputs/Processed_datasets/D2/filtered_expression")
write.table(genes_fly,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/feature_data.txt")
write.table(filtered_reads,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/reads.txt")
write.table(filtered_lengths,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/lengths.txt")
write.table(filtered_tpms,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/tximport_tpms.txt")
write.table(filtered_tpms_1,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/built_tpms.txt")
write.table(metadata,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/metadata.txt")


cat("**** END OF JOB DS2_04_Filter_Genes_Expression ****\n")
