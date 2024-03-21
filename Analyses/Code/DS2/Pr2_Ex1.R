# Load dependencies.
source("Analyses/Code/Dependencies.R")

######################
#### EXERCISE 1.1 ####
######################

# Load data
base::load("Analyses/Code/DS2/filteredBiotypes.RData")

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

get_tpm_from_fpkms <- function(fpkms)
{
  # In this case, the denominator is the sum of fpkms over each sample; it will be zero only if all the gene expressions are 0 (unlikely)
  
  # Build a matrix for the denominator (sum of fpkms in each sample).
  sampleSum <- apply(fpkms,2,sum)
  denoms_mat <- matrix(rep(sampleSum,each=nrow(fpkms)),ncol=ncol(fpkms))/1E6 # The million was not in the slides, but should be?
  
  tpm <- fpkms/denoms_mat
  
  return(tpm)
  
}

tpms_1 <- get_tpm_from_counts(counts=reads,effective_lengths=lengths)
fpkms <- get_fpkm_from_counts(counts=reads,effective_lengths=lengths)
tpms_2 <- get_tpm_from_fpkms(fpkms)

#Validate that tpms_1 = tpms_2
aux_validation_difference <- function (x,threshold=1E-6) {
  val <- 0
  if(x > threshold) {
    val <- 1
  }
  return(val)
}
if(sum(apply(abs(tpms_2-tpms_1),c(1,2),aux_validation_difference))>0) {
  cat("VALIDATION | NOK - Ex1.1: tpms_1 != tpms_2\n")
} else{
  cat("VALIDATION | OK - Ex1.1: tpms_1 = tpms_2\n")
}

rm(list = setdiff(ls(),c("tpms_1","tpms_2"))) # TERMINAR ESTO.

######################
#### EXERCISE 1.2 ####
######################



base::load("Analyses/Code/DS2/filtNormExpression.RData")

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
