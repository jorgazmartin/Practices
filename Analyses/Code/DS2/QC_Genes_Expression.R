source("Analyses/Code/Dependencies.R")

get_effective_counts_from_counts <- function(counts, lengths,effective_lengths,offset=1E-6)
{
  local_lengths=effective_lengths
  for(i in 1:ncol(local_lengths))
    local_lengths[,i]=lengths
  
  issues_lengths=length(which(effective_lengths==0))
  if(issues_lengths>0){
    print(paste(issues_lengths," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths=effective_lengths+offset
  }
  eff_counts=counts*lengths/effective_lengths
  return(eff_counts)
}

get_cpm_from_counts <- function(counts,offset=1E-6)
{
  depths <- apply(counts,2,sum)
  depths_mat=matrix(rep(depths,each=nrow(counts)),ncol=ncol(counts))
  
  issues_depths=which(depths==0)
  if(length(issues_depths)>0){
    print(paste(issues_depths," empty libraries found; CPMs set to 0"))
  }
  
  cpms=1E6*counts/(depths_mat)
  
  if(length(issues_depths)>0)
    cpms[,issues_depths]=0
  
  return(cpms)
}

get_fpkm_from_counts <- function(counts, effective_lengths,offset=1E-6)
{
  
  issues_lengths=length(which(effective_lengths==0))
  if(issues_lengths>0){
    print(paste(issues_lengths," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths=effective_lengths+offset
  }
  
  depths <- apply(counts,2,sum)
  depths_mat=matrix(rep(depths,each=nrow(counts)),ncol=ncol(counts))
  
  issues_depths=which(depths==0)
  if(length(issues_depths)>0){
    print(paste(issues_depths," empty libraries found; FPKMs set to 0"))
  }
  
  fpkms=1E9*counts/(effective_lengths*depths_mat)
  
  if(length(issues_depths)>0)
    fpkms[,issues_depths]=0
  
  return(fpkms)
}

get_tpm_from_counts <- function(counts, effective_lengths,offset=1E-6)
{
  issues_lengths=length(which(effective_lengths==0))
  if(issues_lengths>0){
    print(paste(issues_lengths," features have zero effective length: adding an arbitrary offset of ",offset," bp"))
    effective_lengths=effective_lengths+offset
  }
  rates = counts/effective_lengths
  denoms = apply(rates,2,sum)
  denoms_mat=matrix(rep(denoms,each=nrow(rates)),ncol=ncol(rates))
  tpms=1E6*rates/denoms_mat
  return(tpms)
}

genes_fly$lengths=genes_fly$end_position-genes_fly$start_position


ecs=get_effective_counts_from_counts(counts=reads,lengths=genes_fly$lengths,effective_lengths=lengths)
cpms=get_cpm_from_counts(counts=reads)
fpkms=get_fpkm_from_counts(counts=reads,effective_lengths=lengths)
tpms_1=get_tpm_from_counts(counts=reads,effective_lengths=lengths)


## We will use TPM

## Tentative threshold: at least median tpm>=1 in at least 1 group:

filter=function(tpms,threshold,metadata,grouping_factor){
  number_of_groups=length(unique(metadata[,grouping_factor]))
  medians=tpms[,1:number_of_groups]
  colnames(medians)=unique(metadata[,grouping_factor])
  ## This is unnecessary, but just for clarifying the fact that at this point edians is just an empty container:
  for(i in 1:number_of_groups)
    medians[,i]=0
  for(i in 1:number_of_groups)
  {
    set=which(metadata[,grouping_factor]==colnames(medians)[i])
    chunk=tpms[,set]
    medians[,i]=apply(chunk,1,median)
  }
  medians$max_median=apply(medians,1,max)
  genes_to_keep=rownames(medians)[which(medians$max_median>threshold)]
  return(genes_to_keep)
}

genes_to_keep=filter(tpms=tpms,threshold=1,metadata=metadata,grouping_factor="Condition")
genes_fly=genes_fly[genes_to_keep,]

filtered_reads=reads[genes_to_keep,]
filtered_lengths=lengths[genes_to_keep,]
filtered_tpms=tpms[genes_to_keep,]

filtered_ecs=ecs[genes_to_keep,]
filtered_cpms=cpms[genes_to_keep,]
filtered_fpkms=fpkms[genes_to_keep,]
filtered_tpms_1=tpms_1[genes_to_keep,]


## CPMs normalized
dge <- DGEList(counts = filtered_reads)
dge <- calcNormFactors(dge)

v=voom(dge)
TMM_cpms=v$E

## Let us check the library sizes and normalization factors:
## dge$samples
## Why no sample has a norm.factor of 1? Where is the reference? See calcNormFactors description for answer: "For symmetry, normalization factors are adjusted so their product becomes 1"
prod(dge$samples$norm.factors)