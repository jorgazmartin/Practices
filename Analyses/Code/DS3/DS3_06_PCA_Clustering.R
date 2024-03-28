cat("**** START OF JOB DS3_06_PCA_Clustering ****\n")
rm(list = ls())
source("Analyses/Code/00_Dependencies.R")

validationOK <- TRUE

###################
#### LOAD DATA ####
###################

TMM_cpms <- read.table("Analyses/Inputs/Processed_datasets/D3/filtered/TMM_cpms.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D3/filtered/metadata.txt")

#############
#### PCA ####
#############

## To get PCs in the vector space where genes define the axis (and thus samples are datapoints), prcomp should see the input matrix with genes in columns and samples in the rows:
pca <- prcomp(t(TMM_cpms))
sum_pca <- data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos) <- paste0("PC",c(1:ncol(datos)))

if(length(which(rownames(datos)!=rownames(metadata)))!=0){
  cat("VALIDATION | NOK - Samples (rows) from PCA data are not aligned with metadata\n")
  validationOK <- FALSE
} else {
  cat("VALIDATION | OK - Samples (rows) from PCA data are aligned with metadata. Join is safe.\n")
}

## Join information of PCA results to the metadata.
datos <- cbind(metadata,datos[,1:6])

## Define the factors
datos$species <- factor(datos$species,levels=unique(datos$species))
species_shape <- c(16,17)

datos$tissue <- factor(datos$tissue,levels=unique(datos$tissue))
colors_tissue=c("brown4", "blueviolet","azure4","aquamarine3","bisque3","darkorange2",
              "darkolivegreen","cyan3","chartreuse2","hotpink2","lightpink2","red","yellow3")

## Plots
PC1PC2_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=tissue,shape=species),size=3)+
  scale_colour_manual(values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))

PC1PC2_plot

PC1PC3_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=tissue,shape=species),size=3)+
  scale_colour_manual(values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))

PC1PC3_plot

###########################
#### GRAPHICAL CLUSTER ####
###########################

# Cluster=2, let's see if they group by species?
nclusters <-  2
set.seed(123)

## Perform clustering (kmeans)
clustering <- kmeans(t(TMM_cpms),centers=nclusters,nstart=10)
#print(names(clustering))

## Save the result in the PCA+metadata object:
datos$cluster=paste0("Cluster_",clustering$cluster)

species_shape <- c(21,24)

## Plot.
PC1PC2_plot_w_kmeans_clusters=ggplot(datos)+
  geom_point(aes(x=PC1,y=PC2,shape=species,color=cluster,fill=tissue),size=3,stroke=1.5)+
  scale_fill_manual("tissues",values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))+
  guides(fill = guide_legend("tissues", override.aes = list(shape = 21)))


PC1PC2_plot_w_kmeans_clusters

PC1PC3_plot_w_kmeans_clusters=ggplot(datos)+
  geom_point(aes(x=PC1,y=PC3,shape=species,color=cluster,fill=tissue),size=3,stroke=1.5)+
  scale_fill_manual("tissues",values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))+
  guides(fill = guide_legend("tissues", override.aes = list(shape = 21)))


PC1PC3_plot_w_kmeans_clusters

##############################
#### HIERARCHICAL CLUSTER ####
##############################

# Cluster=2, let's see if they group by species?
clusters <- hclust(dist(t(TMM_cpms)))
clusterCut <- cutree(clusters, 2)

clusterCut

plot(clusters)

#####################
#### TSNE & UMAP ####
#####################

set.seed(123456)
tsne_object <- Rtsne(t(TMM_cpms), check_duplicates=FALSE, theta=0,perplexity=3,max_iter=100000,verbose=FALSE)

tsne_tab <- tsne_object$Y
colnames(tsne_tab) <- c("Tsne_1","Tsne_2")
datosb <- cbind(datos,tsne_tab)

pl_tsne_w_clusters <- ggplot(datosb)+
  geom_point(aes(x=Tsne_1,y=Tsne_2,shape=species,color=cluster,fill=tissue),size=3,stroke=1.5)+
  scale_fill_manual("tissues",values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab("tsne_1")+ylab("tsne_2")+
  guides(fill = guide_legend("tissues", override.aes = list(shape = 21)))

pl_tsne_w_clusters

set.seed(1234567)

umap_object <- umap(t(TMM_cpms),n_neighbors=4)
umap_tab <- umap_object$layout
colnames(umap_tab) <- c("umap_1","umap_2")
datos_umap <- cbind(datos,umap_tab)

pl_umap_w_clusters <- ggplot(datos_umap)+
  geom_point(aes(x=umap_1,y=umap_2,shape=species,color=cluster,fill=tissue),size=3,stroke=1.5)+
  scale_fill_manual("tissues",values=colors_tissue)+scale_shape_manual(values=species_shape)+
  xlab("umap_1")+ylab("umap_2")+
  guides(fill = guide_legend("tissues", override.aes = list(shape = 21)))

pl_umap_w_clusters

##########################
#### SAVE OUTPUT DATA ####
##########################

dir.create("Analyses/Outputs/D3/1_dim_reduction/",recursive=TRUE)

pdf("Analyses/Outputs/D3/1_dim_reduction/PC1_vs_PC2_w_kmeans_cluster.pdf",height=7,width=8)
plot(PC1PC2_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D3/1_dim_reduction/PC1_vs_PC3_w_kmeans_cluster.pdf",height=7,width=8)
plot(PC1PC3_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D3/1_dim_reduction/clusters.pdf")
plot(clusters)
dev.off()

pdf("Analyses/Outputs/D3/1_dim_reduction/clusters_cut.pdf")
plot(clusterCut)
dev.off()

pdf("Analyses/Outputs/D3/1_dim_reduction/tsne_w_kmeans_cluster.pdf",height=7,width=8)
plot(pl_tsne_w_clusters)
dev.off()

pdf("Analyses/Outputs/D3/1_dim_reduction/umap_w_kmeans_cluster.pdf",height=7,width=8)
plot(pl_umap_w_clusters)
dev.off()

cat("**** END OF JOB DS3_06_PCA_Clustering ****\n")
