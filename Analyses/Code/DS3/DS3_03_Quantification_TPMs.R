cat("**** START OF JOB DS3_03_Quantification_TPMs ****\n")
rm(list = ls())
source("Analyses/Code/00_Dependencies.R")

###################
#### FUNCTIONS ####
###################

# Function to compute the tpms:
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

# Load data from the harmonized files:
reads <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/reads.txt")
eff_lengths <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/lengths.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D3/annotated/metadata.txt")

##############
#### MAIN ####
##############

tpms <- get_tpm_from_counts(counts=reads,effective_lengths=eff_lengths)

# Save data
write.table(tpms,"Analyses/Inputs/Processed_datasets/D3/annotated/tpms.txt")

cat("**** END OF JOB DS3_03_Quantification_TPMs ****\n")