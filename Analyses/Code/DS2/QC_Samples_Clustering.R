source("Analyses/Code/Dependencies.R")

nclusters = 4
set.seed(123)
clustering = kmeans(t(TMM_cpms),centers=nclusters,nstart=10)
print(names(clustering))

datos$cluster=paste0("Cluster_",clustering$cluster)

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

clusters <- hclust(dist(t(TMM_cpms)))
clusterCut <- cutree(clusters, 4)

clusterCut

plot(clusters)

set.seed(123456)
tsne_object <- Rtsne(t(TMM_cpms), check_duplicates=FALSE, theta=0,perplexity=3,max_iter=100000,verbose=FALSE)

tsne_tab=tsne_object$Y
colnames(tsne_tab)=c("Tsne_1","Tsne_2")
datosb=cbind(datos,tsne_tab)
pl_tsne_w_clusters=ggplot(datosb)+geom_point(aes(x=Tsne_1,y=Tsne_2,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+xlab("tsne_1")+ylab("tsne_2")

pl_tsne_w_clusters


### check what happens changing parameters (or even random seeds!)
### https://cran.r-project.org/web/packages/umap/vignettes/umap.html

set.seed(1234567)

umap_object = umap(t(TMM_cpms),n_neighbors=4)
umap_tab=umap_object$layout
colnames(umap_tab)=c("umap_1","umap_2")
datos_umap=cbind(datos,umap_tab)

pl_umap_w_clusters=ggplot(datos_umap)+geom_point(aes(x=umap_1,y=umap_2,color=cluster,fill=Condition),shape=21,size=3,stroke=1.5)+
  scale_fill_manual(values=colors_base)+xlab("umap_1")+ylab("umap_2")

pl_umap_w_clusters
