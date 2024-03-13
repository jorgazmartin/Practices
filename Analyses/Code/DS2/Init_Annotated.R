source("Analyses/Code/Dependencies.R")

reads <- read.table("Analyses/Inputs/Raw_datasets/D2/reads.txt")

## This browses ensembl, create, and save the local database, if it is ran the first time. If not, it loads the previously saved object.  

first_time <- TRUE
if(first_time)
{
  fly <- useMart(host="https://feb2021.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="dmelanogaster_gene_ensembl")
  dir.create("Analyses/Inputs/Raw_datasets/D2/marts",recursive=TRUE)
  save(fly,file="Analyses/Inputs/Raw_datasets/D2/marts/fly_mart.Rdata")
} else {
  base::load(file="Analyses/Inputs/Raw_datasets/D2/marts/fly_mart.Rdata")
}

## Now run this to get a sense of the type of attributes available from the database
possible_attributes <- listAttributes(fly)
head(possible_attributes)

## After inspecting that, we identify the attributes we want: 
## Now, run getBM: This means: go to the fly database, select all the items for which "external_gene_id" is included in rownames(reads), and, for all of them, give me the following attributes: c("ensembl_gene_id","flybase_gene_id","start_position","end_position","chromosome_name","gene_biotype")

genes_fly <- getBM(attributes = c("flybase_gene_id","start_position","end_position","chromosome_name","gene_biotype"), filters = "flybase_gene_id", values = rownames(reads) , mart = fly)

## Inspect the header of the result:
head(genes_fly)

## How many of our genes we caught?
length(which(!(rownames(reads) %in% genes_fly$flybase_gene_id)))

# There is a small discrepancy between the database and the genes that were annotated from the kallisto+tximport steps, from the reference transcriptome, and transcript to gene mappings: 
## 43 genes in the reads matrix are absent from the metadata: for the sake of progressing quickly to the points we want to illustrate, we will just drop them from the reads matrix:
rownames(genes_fly) <- genes_fly$flybase_gene_id
reads <- reads[which(rownames(reads) %in% rownames(genes_fly)),]

## Is the mapping leading to the feature data annotations unequivocal?
dim(genes_fly)
length(unique(genes_fly$flybase_gene_id))

# It is ok: the same 14197 genes, let us declare the rownames as flybase_gene_id, and order the items alphabetically

genes_fly <- genes_fly[order(rownames(genes_fly)),]

## Order reads rows also alphabetically:
reads <- reads[order(rownames(reads)),]
## Check that the orders match.
length(which(rownames(reads) != rownames(genes_fly)))

## Now save the reads table as well as the feature table after removing those 43 un-annotated genes.
dir.create("Analyses/Inputs/Processed_datasets/D2/annotated/",recursive=TRUE)
write.table(reads,"Analyses/Inputs/Processed_datasets/D2/annotated/reads.txt")
write.table(genes_fly,"Analyses/Inputs/Processed_datasets/D2/annotated/feature_data.txt")

## Load lengths, tpms and metadata and save the corresponding annotated versions in the same place

lengths=read.table("Analyses/Inputs/Raw_datasets/D2/lengths.txt")
tpms=read.table("Analyses/Inputs/Raw_datasets/D2/tpms.txt")
metadata=read.table("Analyses/Inputs/Raw_datasets/D2/meta_data.txt")

lengths=lengths[which(rownames(lengths) %in%rownames(reads)),]
tpms=tpms[which(rownames(tpms) %in% rownames(reads)),]

dim(lengths)
dim(tpms)

tpms=tpms[order(rownames(tpms)),]
lengths=lengths[order(rownames(lengths)),]

length(which(rownames(reads) != rownames(lengths)))

length(which(rownames(reads) != rownames(tpms)))

## All is ok.

## Check order between samples (column-wise) and meta_data annotations (row-wise)
length(which(rownames(metadata) != colnames(reads)))
length(which(rownames(metadata) != colnames(lengths)))
length(which(rownames(metadata) != colnames(tpms)))

write.table(lengths,"Analyses/Inputs/Processed_datasets/D2/annotated/lengths.txt")
write.table(tpms,"Analyses/Inputs/Processed_datasets/D2/annotated/tpms.txt")
write.table(metadata,"Analyses/Inputs/Processed_datasets/D2/annotated/metadata.txt")
