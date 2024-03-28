cat("**** START OF JOB DS3_05_CPMs_Normalized ****\n")
rm(list = ls())
source("Analyses/Code/00_Dependencies.R")

###################
#### LOAD DATA ####
###################

filtered_reads <- read.table("Analyses/Inputs/Processed_datasets/D3/filtered/reads.txt")

##############
#### MAIN ####
##############

## CPMs normalized
dge <- DGEList(counts = filtered_reads)
dge <- calcNormFactors(dge)

v <- voom(dge)
TMM_cpms <- v$E

## validation. Why no sample has a norm.factor of 1? Where is the reference? See calcNormFactors description for answer: "For symmetry, normalization factors are adjusted so their product becomes 1"
if(abs(prod(dge$samples$norm.factors)-1)>0.001){
  cat("VALIDATION | NOK - Normalized factors product not 1\n")
} else{
  # Save.
  write.table(TMM_cpms,"Analyses/Inputs/Processed_datasets/D3/filtered/TMM_cpms.txt")
}

cat("**** END OF JOB DS3_05_CPMs_Normalized ****\n")