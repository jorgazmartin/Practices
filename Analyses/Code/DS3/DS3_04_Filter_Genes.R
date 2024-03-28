cat("**** START OF JOB DS3_04_Filter_Genes ****\n")
rm(list = ls())
source("Analyses/Code/00_Dependencies.R")

###################
#### FUNCTIONS ####
###################

## Adapted to specified requirements.
filter <- function(tpms,threshold,metadata,grouping_factor){
  
  number_of_groups <- length(unique(metadata[,grouping_factor]))
  # Create medians matrix as an empty container.
  medians <- matrix(0,nrow=dim(tpms)[1],ncol=number_of_groups)
  colnames(medians) <- unique(metadata[,grouping_factor])
  rownames(medians) <- rownames(tpms)
  # Fill medians data.
  for(i in 1:number_of_groups)
  {
    set <- which(metadata[,grouping_factor]==colnames(medians)[i])
    chunk <- tpms[,set]
    medians[,i] <- apply(chunk,1,median)
  }
  # We only keep those whose maximum median is bigger than the threshold.
  medians <- cbind(medians,num_success=apply(medians,1,function (x) length(which(x >= threshold))))
  genes_to_keep=rownames(medians)[which(medians[,"num_success"] == number_of_groups)]
  
  return(genes_to_keep)
}


###################
#### LOAD DATA ####
###################

# Load data from the harmonized files:
reads <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/reads.txt")
eff_lengths <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/lengths.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/metadata.txt")
tpms <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/tpms.txt")

####################
#### BY BIOTYPE ####
####################

## Duda.

############################
#### BY GENE EXPRESSION ####
############################

genes_to_keep <- filter(tpms=tpms,threshold=1.0,metadata=metadata,grouping_factor="species")
filtered_reads <- reads[genes_to_keep,]
filtered_eff_lengths <- eff_lengths[genes_to_keep,]
filtered_tpms <- tpms[genes_to_keep,]

# Save data.
dir.create("Analyses/Inputs/Processed_datasets/D3/filtered/",recursive=TRUE)
write.table(filtered_reads,"Analyses/Inputs/Processed_datasets/D3/filtered/reads.txt")
write.table(filtered_eff_lengths,"Analyses/Inputs/Processed_datasets/D3/filtered/lengths.txt")
write.table(filtered_tpms,"Analyses/Inputs/Processed_datasets/D3/filtered/tpms.txt")
write.table(metadata,"Analyses/Inputs/Processed_datasets/D3/filtered/metadata.txt")

cat("**** END OF JOB DS3_04_Filter_Genes ****\n")
