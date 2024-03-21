# Load dependencies.
source("Analyses/Code/Dependencies.R")

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


# Load data
base::load("Analyses/Code/DS2/filteredBiotypes.RData")

# Calculate the length of each gene. Store in genes_fly featured data.
genes_fly$lengths <- genes_fly$end_position-genes_fly$start_position

# Quantification. Actually, from now on we will use just TPM.
ecs <- get_effective_counts_from_counts(counts=reads,lengths=genes_fly$lengths,effective_lengths=lengths)
cpms <- get_cpm_from_counts(counts=reads)
fpkms <- get_fpkm_from_counts(counts=reads,effective_lengths=lengths)
tpms_1 <- get_tpm_from_counts(counts=reads,effective_lengths=lengths)

# Filter: with the criteria of the function defined above.
genes_to_keep <- filter(tpms=tpms,threshold=1,metadata=metadata,grouping_factor="Condition")
genes_fly <- genes_fly[genes_to_keep,]
filtered_reads <- reads[genes_to_keep,]
filtered_lengths <- lengths[genes_to_keep,]
filtered_tpms <- tpms[genes_to_keep,]

# Extra:
filtered_ecs=ecs[genes_to_keep,]
filtered_cpms=cpms[genes_to_keep,]
filtered_fpkms=fpkms[genes_to_keep,]
filtered_tpms_1=tpms_1[genes_to_keep,]


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
  save(filtered_cpms,filtered_ecs,filtered_fpkms,filtered_lengths,filtered_reads,filtered_tpms,filtered_tpms_1,genes_fly,metadata,
       dge,TMM_cpms,file="Analyses/Code/DS2/filtNormExpression.RData")
}

