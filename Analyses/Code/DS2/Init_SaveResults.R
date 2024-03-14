## Save filtered results.

dir.create("Analyses/Inputs/Processed_datasets/D2/filtered/",recursive=TRUE)
## Between 1_mapping and 3_dim_reduction, you should have produced 2_exploration
dir.create("Analyses/Outputs/D2/1_dim_reduction/",recursive=TRUE)

write.table(genes_fly,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_feature_data.txt")
write.table(metadata,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_metadata.txt")
write.table(filtered_reads,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_reads.txt")
write.table(filtered_lengths,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_lengths.txt")
write.table(filtered_tpms,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_tximport_tpms.txt")
write.table(filtered_tpms_1,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_built_tpms.txt")
write.table(filtered_cpms,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_cpms.txt")
write.table(filtered_fpkms,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_fpkms.txt")
write.table(TMM_cpms,"Analyses/Inputs/Processed_datasets/D2/filtered/filtered_TMM_cpms.txt")

## Now, save the dimensionality reduction figures:

pdf("Analyses/Outputs/D2/1_dim_reduction/PC1_vs_PC2_w_kmeans_cluster.pdf",height=4,width=6)
plot(PC1PC2_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/PC1_vs_PC3_w_kmeans_cluster.pdf",height=4,width=6)
plot(PC1PC3_plot_w_kmeans_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/clusters.pdf")
plot(clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/clusters_cut.pdf")
plot(clusterCut)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/tsne_w_kmeans_cluster.pdf",height=4,width=6)
plot(pl_tsne_w_clusters)
dev.off()

pdf("Analyses/Outputs/D2/1_dim_reduction/umap_w_kmeans_cluster.pdf",height=4,width=6)
plot(pl_umap_w_clusters)
dev.off()