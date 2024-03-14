source("Analyses/Code/Dependencies.R")

## To get PCs in the vector space where genes define the axis (and thus samples are datapoints), prcomp should see the input matrix with genes in columns and samples in the rows:
pca <-prcomp(t(TMM_cpms))
sum_pca=data.frame(summary(pca)$importance[,c(1:5)])

datos <- data.frame(pca$x)
colnames(datos)=paste0("PC",c(1:ncol(datos)))
length(which(rownames(datos)!=rownames(metadata)))

datos=cbind(metadata,datos[,1:6])

datos$Condition=factor(datos$Condition,levels=unique(datos$Condition))
colors_base=c("slategray2", "orange2","darkorchid2","springgreen2")


PC1PC2_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC2,color=Condition),size=3)+
  scale_colour_manual(values=colors_base)+xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+ylab(paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"))

PC1PC2_plot

PC1PC3_plot=ggplot(datos)+geom_point(aes(x=PC1,y=PC3,color=Condition),size=3)+
  scale_colour_manual(values=colors_base)+xlab(paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"))+ylab(paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))

PC1PC3_plot

datos$Colorby=rep(colors_base,each=3)
plot3d(datos$PC1,datos$PC2,datos$PC3,col=datos$Colorby,size=15,xlab=paste0("PC1: ",round(100*sum_pca$PC1[2],digits=2),"% variance"),ylab=paste0("PC2: ",round(100*sum_pca$PC2[2],digits=2),"% variance"),zlab=paste0("PC3: ",round(100*sum_pca$PC3[2],digits=2),"% variance"))
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