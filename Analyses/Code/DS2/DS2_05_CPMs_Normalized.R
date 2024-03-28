# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_05_CPMs_Normalized ****\n")
rm(list = ls())

###################
#### LOAD DATA ####
###################

filtered_reads <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_expression/reads.txt")

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
  # Save result
  write.table(TMM_cpms,"Analyses/Inputs/Processed_datasets/D2/filtered_expression/TMM_cpms.txt")
  # Save snapshot to Pr2.Ex1.2.
  save(dge,file="Analyses/Code/DS2/_snapshots/dge.RData")
}

cat("**** END OF JOB DS2_05_CPMs_Normalized ****\n")

