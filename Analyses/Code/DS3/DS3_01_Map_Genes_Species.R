cat("**** START OF JOB DS3_01_Map_Genes_Species ****\n")
rm(list = ls())
source("Analyses/Code/Dependencies.R")

## Load the datasets
mouse_name <- base::load("Analyses/Inputs/Raw_datasets/D3/txi_mouse.Rdata",envir = parent.frame())
human_name <- base::load("Analyses/Inputs/Raw_datasets/D3/txi_human.Rdata")

## Get only the counts matrices for now:
counts_human <- txi_human$counts
counts_mouse <- txi_mouse$counts

## inspects the results: 
head(counts_human)
head(counts_mouse)

## We see that the gene ensembl IDs come with decimal numbers at the end. These indicate the ENSEMBL version (see https://www.ensembl.org/Help/Faq?id=488) and (as discussed here: https://www.biostars.org/p/430644/) it is better to remove them in many contexts. We will do that with the following function, home made:
remove_versions=function(x){
  extract=function(i){strsplit(rownames(x)[i], "[.]")[[1]][1]}
  genes=sapply(1:nrow(x),extract)
  rownames(x)=genes
  return(x)
}

counts_human <- remove_versions(counts_human)
counts_mouse <- remove_versions(counts_mouse)

## Load the marts:

first_time_D3 <- FALSE
if(first_time_D3)
{
  mouse <- useMart(host="https://jan2019.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
  human <- useMart(host="https://jan2019.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  
  ## And we will save them, just in case 
  
  dir.create("Analyses/Inputs/Raw_datasets/D3/marts",recursive=TRUE)
  save(mouse,file="Analyses/Inputs/Raw_datasets/D3/marts/mouse.Rdata")
  save(human,file="Analyses/Inputs/Raw_datasets/D3/marts/human.Rdata")
  
}else{
  base::load("Analyses/Inputs/Raw_datasets/D3/marts/mouse.Rdata")
  base::load("Analyses/Inputs/Raw_datasets/D3/marts/human.Rdata")
}

## Now, call getLDS to retrieve tehe (raw) ID mappings:

genes_map <- getLDS(attributes = c("ensembl_gene_id","mgi_symbol"), filters = "ensembl_gene_id", values=rownames(counts_mouse) , mart = mouse, attributesL = c("ensembl_gene_id","hgnc_symbol"), martL = human, uniqueRows=T)

## Check the looks of this file: are the mappings unique?
head(genes_map)

## Since data in both matrix is mapped against ENSEMBL-annotated features, the important mapping is the one between these stable IDs:
genes_map$map <- paste0(genes_map$Gene.stable.ID.1,"_",genes_map$Gene.stable.ID)
length(unique(genes_map$map))
dim(genes_map)
dupes_mouse <- unique(genes_map$Gene.stable.ID[which(duplicated(genes_map$Gene.stable.ID))])
length(dupes_mouse)
dupes_human <- unique(genes_map$Gene.stable.ID.1[which(duplicated(genes_map$Gene.stable.ID.1))])
length(dupes_human)

## That means that for 1408 mouse ensemble gene IDs, there's more than one human correspondence. And the inverse is true (that there is more than one mouse ensemble ID mapping) for 1087 human ensembl IDs

## We will not dig further looking for more precise approaches, ortholog-homolog mapping is tricky, and beyond the scope of the course for us: let us play it safe: keep only the mappings containing unequivocal human-to-mouse mappings.

genes_map <- genes_map[which(   !(genes_map$Gene.stable.ID %in% dupes_mouse)  &  !(genes_map$Gene.stable.ID.1 %in% dupes_human)),]
dim(genes_map)

## There are 15842 genes with bijective human-to-mouse mappings. It is quite less than the original number of features in each of the counts tables, but let us stick with that.

## Let us check that now the ensembl gene IDs are all different among them, in each species:

length(unique(genes_map$Gene.stable.ID))

length(unique(genes_map$Gene.stable.ID.1))

## Let us take a look to common gene names (MGI and HGNC) Here, it is also normal that a single ENSEMBl id maps to several common names: for those cases where the mapping is not unequivocal, we will just remove the common name and substitute it with the ensembl ID. We do that too when there is no common name reported.

length(unique(genes_map$MGI.symbol))

non_unequivocal_mouse_names <- genes_map$MGI.symbol[which(duplicated(genes_map$MGI.symbol))]
genes_map$MGI.symbol[which(genes_map$MGI.symbol %in% non_unequivocal_mouse_names)] <- genes_map$Gene.stable.ID[which(genes_map$MGI.symbol %in% non_unequivocal_mouse_names)]

length(unique(genes_map$MGI.symbol))

## Now all are different
length(which(genes_map$MGI.symbol==""))

## Also, there are no gaps.

length(unique(genes_map$HGNC.symbol))
non_unequivocal_human_names <- genes_map$HGNC.symbol[which(duplicated(genes_map$HGNC.symbol))]
genes_map$HGNC.symbol[which(genes_map$HGNC.symbol %in% non_unequivocal_human_names)] <- genes_map$Gene.stable.ID.1[which(genes_map$HGNC.symbol %in% non_unequivocal_human_names)]

length(unique(genes_map$HGNC.symbol))
length(which(genes_map$HGNC.symbol==""))
length(unique(genes_map$Gene.stable.ID))
length(unique(genes_map$Gene.stable.ID.1))
length(unique(genes_map$MGI.symbol))
length(unique(genes_map$HGNC.symbol))

## Let us save the resulting table:

genes_map <- genes_map[,1:4]
colnames(genes_map) <- c("Ensembl_mouse","Gene_name_mouse","Ensembl_human","Gene_name_human")
write.table(genes_map,"Analyses/Inputs/Raw_datasets/D3/gene_mapping.txt")

cat("**** END OF JOB DS3_01_Map_Genes_Species ****\n")
