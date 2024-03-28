# Load the dependencies.
source("Analyses/Code/00_Dependencies.R")

cat("**** START OF JOB DS2_07_PCA_Clustering ****\n")
rm(list = ls())

validationOK <- TRUE

###################
#### LOAD DATA ####
###################

TMM_cpms <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_expression/TMM_cpms.txt")
metadata <- read.table("Analyses/Inputs/Processed_datasets/D2/filtered_expression/metadata.txt")

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

# Join information of PCA results to the metadata.
datos <- cbind(metadata,datos[,1:6])

datos$Condition=factor(datos$Condition,levels=unique(datos$Condition))
colors_base=c("slategray2", "orange2","darkorchid2","springgreen2")


PC1PC2_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Condition),size=3)+
  scale_colour_manual(values=colors_base)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))

PC1PC2_plot

PC1PC3_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Condition),size=3)+
  scale_colour_manual(values=colors_base)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))

PC1PC3_plot

datos$Colorby=rep(colors_base,each=3)
plot3d(datos$PC1,datos$PC2,datos$PC3,col=datos$Colorby,size=15,
       xlab=paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"),
       ylab=paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"),
       zlab=paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))
legend3d("topright", legend = unique(datos$Condition), col = unique(datos$Colorby), pch=19)

## Check the projections that reproduce the previous plot, tinker the good one and save it with:
view <- par3d("userMatrix")

## Later you can restore it by running the plot and then: 
par3d(userMatrix = view)


## To save the view to a figure file, we do:
dir.create("Analyses/Outputs/D2/1_dim_reduction/",recursive=TRUE)

rgl.postscript("Analyses/Outputs/D2/1_dim_reduction/PCA_3d.pdf",fmt="pdf")
rgl.snapshot("Analyses/Outputs/D2/1_dim_reduction/PCA_3d.png",fmt="png")
# Let us save the view as well:
save(view,file="Analyses/Outputs/D2/1_dim_reduction/PCA_3d_view.Rdata")

###########################
#### GRAPHICAL CLUSTER ####
###########################

nclusters <- 4
set.seed(123)

## Perform clustering (kmeans)
clustering <- kmeans(t(TMM_cpms),centers=nclusters,nstart=10)
#print(names(clustering))

## Save the result in the PCA+metadata object:
datos$cluster=paste0("Cluster_",clustering$cluster)

species_shape <- c(21,24)

## Plot.
PC1PC2_plot_w_kmeans_clusters=ggplot(datos)+
  geom_point(aes(x=PC1,y=PC2,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))


PC1PC2_plot_w_kmeans_clusters

## Same for PC1 vs PC3
PC1PC3_plot_w_kmeans_clusters=ggplot(datos)+
  geom_point(aes(x=PC1,y=PC3,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+
  xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+
  ylab(paste0("PC2: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))


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

pl_tsne_w_clusters=ggplot(datosb)+
  geom_point(aes(x=Tsne_1,y=Tsne_2,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+xlab("tsne_1")+ylab("tsne_2")

pl_tsne_w_clusters

set.seed(1234567)

umap_object <- umap(t(TMM_cpms),n_neighbors=4)
umap_tab <- umap_object$layout
colnames(umap_tab) <- c("umap_1","umap_2")
datos_umap <- cbind(datos,umap_tab)

pl_umap_w_clusters=ggplot(datos_umap)+
  geom_point(aes(x=umap_1,y=umap_2,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+xlab("umap_1")+ylab("umap_2")

pl_umap_w_clusters

##########################
#### SAVE OUTPUT DATA ####
##########################

dir.create("Analyses/Outputs/D2/1_dim_reduction/",recursive=TRUE)

pdf("Analyses/Outputs/D2/1_dim_reduction/PC1_vs_PC2_w_kmeans_cluster.pdf",height=7,width=8)
plot(PC1PC2_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/PC1_vs_PC3_w_kmeans_cluster.pdf",height=7,width=8)
plot(PC1PC3_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/clusters.pdf")
plot(clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/clusters_cut.pdf")
plot(clusterCut)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/tsne_w_kmeans_cluster.pdf",height=7,width=8)
plot(pl_tsne_w_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/umap_w_kmeans_cluster.pdf",height=7,width=8)
plot(pl_umap_w_clusters)
dev.off()

cat("**** END OF JOB DS2_07_PCA_Clustering ****\n")