# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_06_Pr2_Ex1 ****\n")
rm(list = ls())

###################
#### FUNCTIONS ####
###################

get_tpm_from_fpkms <- function(fpkms)
{
  # In this case, the denominator is the sum of fpkms over each sample; it will be zero only if all the gene expressions are 0 (unlikely)
  
  # Build a matrix for the denominator (sum of fpkms in each sample).
  sampleSum <- apply(fpkms,2,sum)
  denoms_mat <- matrix(rep(sampleSum,each=nrow(fpkms)),ncol=ncol(fpkms))/1E6 # The million was not in the slides, but should be?
  
  tpm <- fpkms/denoms_mat
  
  return(tpm)
  
}

aux_validation_difference <- function (x,threshold=1E-6) {
  val <- 0
  if(x > threshold) {
    val <- 1
  }
  return(val)
}

######################
#### EXERCISE 1.1 ####
######################

## Load data
tpms_1 <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/built_tpms.txt")
fpkms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_biotype/fpkms.txt")

## Calculation
tpms_2 <- get_tpm_from_fpkms(fpkms)

## Validation
if(sum(apply(abs(tpms_2-tpms_1),c(1,2),aux_validation_difference))>0) {
  cat("VALIDATION | NOK - Ex1.1: tpms_1 != tpms_2\n")
} else{
  cat("VALIDATION | OK - Ex1.1: tpms_1 = tpms_2\n")
}

######################
#### EXERCISE 1.2 ####
######################

TMM_cpms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_expression/TMM_cpms.txt")
filtered_reads <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_expression/reads.txt")
base::load("Analyses/Code/DS2/_snapshots/dge.RData")

gene_name <- rownames(filtered_reads)[1]
sample_name <- colnames(filtered_reads)[1]

R_1_1 <- filtered_reads[gene_name,sample_name]
TMM_1 <- dge$samples[sample_name,"norm.factors"]
N_1 <- dge$samples[sample_name,"lib.size"]
log_cpm_1_1 <- log((R_1_1+0.5)*1E6/(TMM_1*N_1+1),base=2)

if(abs(TMM_cpms[gene_name,sample_name]-log_cpm_1_1)>1E-6) {
  cat("VALIDATION | NOK - Ex1.2:\n",
      "Black box value: ",TMM_cpms[gene_name,sample_name],"\n",
      "Estimated value: ",log_cpm_1_1,"\n",sep="")
} else {
  cat("VALIDATION | OK - Ex1.2:\n",
      "Black box value: ",TMM_cpms[gene_name,sample_name],"\n",
      "Estimated value: ",log_cpm_1_1,"\n",sep="")
}

cat("**** END OF JOB DS2_06_Pr2_Ex1 ****\n")
