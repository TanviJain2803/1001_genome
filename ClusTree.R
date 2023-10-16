#ClusTree Script
library(Seurat)
library(Matrix)
library(stringr)
library(reshape2)
library(ggplot2)
library(mrtree)
library(gplots)
#load seurat obj with resolutions
load("/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster.RData")

#load csv file with run2.2 root colnames (from 10x root obj intersection)
root_cell_types <- read.csv("/gale/ddn/ddn_neomorph/tjain/1001_genome/tenx_root_cell_types.csv")
root_cell_types <- root_cell_types[,1]

#transfering run2.2 root nuclei into Genomes_cluster metadata
Genomes_Cluster[['run2.2_root']] <- as.integer(colnames(Genomes_Cluster) %in% root_cell_types)
Idents(Genomes_Cluster)<- Genomes_Cluster@meta.data$run2.2_root

#accounting for all root nuclei from root nuclei in run 2.2 and cluster res 1
temp_df <- Genomes_Cluster[[]][c('run2.2_root', "integrated_snn_res.1")]
num_cluster <- max(as.numeric(unlist(Genomes_Cluster[['integrated_snn_res.1']])))
root_nuclei <- matrix(NA, nrow=num_cluster+1,ncol=4)
for(i in 0:num_cluster){
  root_nuclei[i+1,1] <- i #cluster num
  root_nuclei[i+1,2] <-nrow(temp_df[temp_df$integrated_snn_res.1==i,]) #num total nuclei in current cluster
  root_nuclei[i+1,3]<- sum(temp_df[temp_df$integrated_snn_res.1==i,]$run2.2_root) #num run2.2 root nuclei in current cluster
}
root_nuclei[,4]<-root_nuclei[,3]/root_nuclei[,2] #ratio of total nuclei: run 2.2 nuclei in cluster
root_clusters_gc <- root_nuclei[root_nuclei[,4]>0,1]  #classifying any cluster with ratio >0 as a root cluster
Genomes_Cluster[['integrated_snn_res.0.02']] <- as.integer(Genomes_Cluster$integrated_snn_res.1 %in% root_clusters_gc) #creating metadata column (this is manual split for root and shoot)
Idents(Genomes_Cluster)<- Genomes_Cluster@meta.data$integrated_snn_res.0.02
Genomes_Cluster[['run2.2_root']] <- NULL #deleting previous run2.2_root metadata column
Idents(Genomes_Cluster)<- Genomes_Cluster@meta.data$integrated_snn_res.0.02

#save new data obj
save(Genomes_Cluster, file = "/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster_split_res.RData")

#creating UMAP with initial manual res split
jpeg(file="/gale/ddn/ddn_neomorph/tjain/root_shoot_split_UMAP.jpeg")
DimPlot(object = Genomes_Cluster, reduction = "umap", label=TRUE, raster=TRUE)
dev.off()

#creating clustree diagram with final resolution 1
out = mrtree(Genomes_Cluster, prefix = 'integrated_snn_res.', n.cores = 2, consensus=FALSE, augment.path=FALSE)
jpeg(file="/gale/ddn/ddn_neomorph/tjain/clustree.jpeg")
plot_tree(labelmat=out$labelmat.mrtree, plot.piechart = TRUE,
          node.size = 0.2, tip.label.dist = 10, bottom.margin=30)
dev.off()
