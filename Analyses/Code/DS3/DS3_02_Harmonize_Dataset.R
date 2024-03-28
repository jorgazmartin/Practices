cat("**** START OF JOB DS3_02_Harmonize_Dataset ****\n")
rm(list = ls())
source("Analyses/Code/00_Dependencies.R")

validationOK <- TRUE

## Load the datasets
mouse_name <- base::load("Analyses/Inputs/Raw_datasets/D3/txi_mouse.Rdata",envir = parent.frame())
human_name <- base::load("Analyses/Inputs/Raw_datasets/D3/txi_human.Rdata")
gene_map <- read.table("Analyses/Inputs/Raw_datasets/D3/gene_mapping.txt")
metadata <- read.table("Analyses/Inputs/Raw_datasets/D3/meta_data.txt")

## Get only the variables of interest for each species
reads_human <- txi_human$counts
reads_mouse <- txi_mouse$counts
lengths_human <- txi_human$length
lengths_mouse <- txi_mouse$length

## Clean row names (removing ENSEMBL version)
remove_versions=function(x){
  extract=function(i){strsplit(rownames(x)[i], "[.]")[[1]][1]}
  genes=sapply(1:nrow(x),extract)
  rownames(x)=genes
  return(x)
}
reads_human <- remove_versions(txi_human$counts)
reads_mouse <- remove_versions(txi_mouse$counts)
lengths_human <- remove_versions(txi_human$length)
lengths_mouse <- remove_versions(txi_mouse$length)

## Filter only ortholog genes (listed in gene_map).
reads_human <- reads_human[which(rownames(reads_human) %in% gene_map$Ensembl_human),]
reads_mouse <- reads_mouse[which(rownames(reads_mouse) %in% gene_map$Ensembl_mouse),]
lengths_human <- reads_human[which(rownames(reads_human) %in% gene_map$Ensembl_human),]
lengths_mouse <- reads_mouse[which(rownames(reads_mouse) %in% gene_map$Ensembl_mouse),]

## Filter validations:
if(dim(reads_human)[1] != dim(gene_map)[1] || dim(reads_mouse)[1] != dim(gene_map)[1] ||
   dim(lengths_human)[1] != dim(gene_map)[1] || dim(lengths_mouse)[1] != dim(gene_map)[1]) {
  cat("VALIDATION | NOK - Different number of rows\n",
      "  Ortholog list: ",dim(gene_map)[1],"\n", "  Reads_human: ",dim(reads_human)[1],"\n",
      "  Lenghts_human: ",dim(lengths_human)[1],"\n", "  Reads_mouse: ",dim(reads_mouse)[1],"\n",
      "  Lengths_mouse: ",dim(lengths_mouse)[1],"\n",sep="")
  validationOK <- FALSE
} else {
  cat("VALIDATION | OK - Same number of rows\n")
}

## As there is a missing human gene, we eliminate it from the other objects.
gene_map <- gene_map[which(gene_map$Ensembl_human %in% rownames(reads_human)),]
reads_mouse <- reads_mouse[which(rownames(reads_mouse) %in% gene_map$Ensembl_mouse),]
lengths_mouse <- reads_mouse[which(rownames(reads_mouse) %in% gene_map$Ensembl_mouse),]

## Validate again:
if(dim(reads_human)[1] != dim(gene_map)[1] || dim(reads_mouse)[1] != dim(gene_map)[1] ||
   dim(lengths_human)[1] != dim(gene_map)[1] || dim(lengths_mouse)[1] != dim(gene_map)[1]) {
  cat("VALIDATION | NOK - Different number of rows\n",
      "  Ortholog list: ",dim(gene_map)[1],"\n", "  Reads_human: ",dim(reads_human)[1],"\n",
      "  Lenghts_human: ",dim(lengths_human)[1],"\n", "  Reads_mouse: ",dim(reads_mouse)[1],"\n",
      "  Lengths_mouse: ",dim(lengths_mouse)[1],"\n",sep="")
  validationOK <- FALSE
} else {
  cat("VALIDATION | OK - Same number of rows\n")
  validationOK <- TRUE
}

## Transform mouse genes (row names) in their human orthologs
rownames(reads_mouse) <- sapply(rownames(reads_mouse),
                                function (id_mouse,gene_map) gene_map[gene_map$Ensembl_mouse == id_mouse,"Ensembl_human"],
                                gene_map = gene_map)
rownames(lengths_mouse) <- sapply(rownames(lengths_mouse),
                                function (id_mouse,gene_map) gene_map[gene_map$Ensembl_mouse == id_mouse,"Ensembl_human"],
                                gene_map = gene_map)

## Order all objects by row name (huma gene id)
reads_mouse <- reads_mouse[order(rownames(reads_mouse)),]
lengths_mouse <- lengths_mouse[order(rownames(lengths_mouse)),]
reads_human <- reads_human[order(rownames(reads_human)),]
lengths_human <- lengths_human[order(rownames(lengths_human)),]

## Validate that the rownames of both human and mouse are the same.

if (length(which(rownames(reads_human) != rownames(reads_mouse))) != 0) {
  validationOK <- FALSE
  cat("VALIDATION | NOK - Unmatched genes in reads table\n")
} else{
  cat("VALIDATION | OK - All genes match in reads table\n")
}

if (length(which(rownames(lengths_human) != rownames(lengths_mouse))) != 0) {
  validationOK <- FALSE
  cat("VALIDATION | NOK - Unmatched genes in lengths table\n")
} else{
  cat("VALIDATION | OK - All genes match in lengths table\n")
}

# If everything OK, it is safe to perfrome the join.
if (validationOK) {
  ## Perform the join of two species.
  reads_joined <- cbind(reads_human,reads_mouse)
  lengths_joined <- cbind(lengths_human,lengths_mouse)
  
  ## Change the name of the rows by human HUGO names.
  rownames(reads_joined) <- sapply(rownames(reads_joined),
                                   function (ensembl_id,gene_map) gene_map[gene_map$Ensembl_human == ensembl_id,"Gene_name_human"],
                                   gene_map = gene_map)
  reads_joined <- reads_joined[order(rownames(reads_joined)),]
  rownames(lengths_joined) <- sapply(rownames(lengths_joined),
                                   function (ensembl_id,gene_map) gene_map[gene_map$Ensembl_human == ensembl_id,"Gene_name_human"],
                                   gene_map = gene_map)
  lengths_joined <- lengths_joined[order(rownames(lengths_joined)),]
  
  ## Order the name of the columns to validate if all the samples are in metadata file.
  reads_joined <- reads_joined[,order(colnames(reads_joined))]
  lengths_joined <- lengths_joined[,order(colnames(lengths_joined))]
  metadata <- metadata[order(rownames(metadata)),]
  
  ## Sample validations
  if (length(which(!(colnames(reads_joined) %in% metadata$Sample_ID))) != 0 ||
      length(which(colnames(reads_joined) != metadata$Sample_ID)) != 0) {
    cat("VALIDATION | NOK - Samples are missing in reads\n")
    validationOK <- FALSE
  } else {
    cat("VALIDATION | OK - Samples in reads (column names) aligned with metadata\n")
  }
  if (length(which(!(colnames(lengths_joined) %in% metadata$Sample_ID))) != 0 ||
      length(which(colnames(lengths_joined) != metadata$Sample_ID)) != 0) {
    cat("VALIDATION | NOK - Samples are missing in lengths\n")
    validationOK <- FALSE
  } else {
    cat("VALIDATION | OK - Samples in lengths (column names) aligned with metadata\n")
  }
  
  if(validationOK){
    dir.create("Analyses/Inputs/Processed_datasets/D3/annotated/",recursive=TRUE)
    write.table(reads_joined,"Analyses/Inputs/Processed_datasets/D3/annotated/reads.txt")
    write.table(lengths_joined,"Analyses/Inputs/Processed_datasets/D3/annotated/lengths.txt")
    write.table(metadata,"Analyses/Inputs/Processed_datasets/D3/annotated/metadata.txt")
    cat("**** END OF JOB DS3_02_Harmonize_Dataset ****\n")
  }
  
}
