# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_03_Quantification ****\n")
rm(list = ls())

###################
#### FUNCTIONS ####
###################

get_effective_counts_from_counts <- function(counts, lengths, effective_lengths, offset=1E-6)
{
  
  # Zeros in the denominator make trouble.  
  issues_lengths <- length(which(effective_lengths==0))
  if(issues_lengths>0){
    print(paste(issues_lengths," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths <- effective_lengths+offset
  }
  
  # Definition of the metric.
  eff_counts <- counts*lengths/effective_lengths
  
  return(eff_counts)
}

get_cpm_from_counts <- function(counts,offset=1E-6)
{
  # Sum by columns (genes in a sample)
  depths <- apply(counts,2,sum)
  
  # Zeros in the denominator make trouble.
  issues_depths <- which(depths==0)
  if(length(issues_depths)>0){
    print("Empty libraries found; CPMs set to 0")
    depths[issues_depths] <- depths[issues_depths]+offset
  }
  
  # Store in a matrix where each column (sample) hast the corresponding depth repeated.
  depths_mat <- matrix(rep(depths,each=nrow(counts)),ncol=ncol(counts))
  
  # Definition of the metric.
  cpms <- 1E6*counts/(depths_mat)
  
  # Undo the offset transformation.
  if(length(issues_depths)>0){
    cpms[,issues_depths] <- 0
  }
  
  return(cpms)
}

get_fpkm_from_counts <- function(counts, effective_lengths,offset=1E-6)
{
  
  # Zeros in the denominator make trouble.
  issues_lengths <- which(effective_lengths==0)
  if(length(issues_lengths>0)){
    print(paste(lengths(issues_lengths)," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths <- effective_lengths+offset
  }
  
  depths <- apply(counts,2,sum)
  depths_mat <- matrix(rep(depths,each=nrow(counts)),ncol=ncol(counts))
  
  # Definition of the metric.
  fpkms <- 1E9*counts/(effective_lengths*depths_mat)
  
  # Undo the offset transformation.
  if(length(issues_lengths)>0){
    fpkms[,issues_lengths] <- 0
  }
  
  return(fpkms)
  
}

get_tpm_from_counts <- function(counts, effective_lengths,offset=1E-6)
{
  
  # Zeros in the denominator make trouble.
  issues_lengths <- which(effective_lengths==0)
  if(length(issues_lengths)>0){
    print(paste(length(issues_lengths)," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths <- effective_lengths+offset
  }
  
  # Metric definition.
  rates <- counts/effective_lengths
  denoms <- apply(rates,2,sum)
  denoms_mat <- matrix(rep(denoms,each=nrow(rates)),ncol=ncol(rates))
  tpms <- 1E6*rates/denoms_mat
  
  # Undo the offset transformation?
  
  return(tpms)
  
}

###################
#### LOAD DATA ####
###################

genes_fly <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/feature_data.txt")
reads <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/reads.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/metadata.txt")
lengths <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/lengths.txt")
tpms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/tximport_tpms.txt")

##############
#### MAIN ####
##############

## Calculate the length of each gene. Store in genes_fly featured data.
genes_fly$lengths <- genes_fly$end_position-genes_fly$start_position

## Quantification. Actually, from now on we will use just TPM.
ecs <- get_effective_counts_from_counts(counts=reads,lengths=genes_fly$lengths,effective_lengths=lengths)
cpms <- get_cpm_from_counts(counts=reads)
fpkms <- get_fpkm_from_counts(counts=reads,effective_lengths=lengths)
tpms_1 <- get_tpm_from_counts(counts=reads,effective_lengths=lengths)

## Save data
write.table(ecs,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/ecs.txt")
write.table(cpms,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/cpms.txt")
write.table(fpkms,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/fpkms.txt")
write.table(tpms_1,"Analyses/Inputs/Processed_datasets/D2/filtered_biotype/built_tpms.txt")

cat("**** START OF JOB DS2_03_Quantification ****\n")
