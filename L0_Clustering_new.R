library(Seurat)
library(Matrix)
library(stringr)

# min cut off 600
Read_Cluster<-function(path, Run_Number) {
  Genomes_Cluster <- CreateSeuratObject(counts = path, min.cells = 5, min.features = 5)
  Genomes_Cluster[["percent.mt"]] <- PercentageFeatureSet(object = Genomes_Cluster, pattern = "ATMG")
  Genomes_Cluster[["percent.ch"]] <- PercentageFeatureSet(object = Genomes_Cluster, pattern = "ATCG")
  Genomes_Cluster$Run <- Run_Number
  Genomes_Cluster <- Genomes_Cluster[-grep("^ATCG", row.names(Genomes_Cluster)),]
  Genomes_Cluster <- Genomes_Cluster[-grep("^ATMG", row.names(Genomes_Cluster)),]
  Genomes_Cluster <- subset(x = Genomes_Cluster, nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 8 & percent.ch < 50)
  Genomes_Cluster <- NormalizeData(object = Genomes_Cluster, verbose = FALSE)
  Genomes_Cluster <- FindVariableFeatures(object = Genomes_Cluster, selection.method = "vst", nfeatures = 6000)
  Genomes_Cluster <- ScaleData(object = Genomes_Cluster,verbose = TRUE)
  Genomes_Cluster <- RunPCA(object = Genomes_Cluster, npcs = 30, verbose = FALSE)
  Genomes_Cluster}

print("start")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_1/CI_Run_1_full_matrix.RData")
Genomes_Cluster_Run_1<-Read_Cluster(full_matrix, "Run_1")
Genomes_Cluster_Run_1

print("run_1")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_3/CI_Run_3_full_matrix.RData")
Genomes_Cluster_Run_3<-Read_Cluster(full_matrix, "Run_3")
Genomes_Cluster_Run_3

print("run_3")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_4/CI_Run_4_full_matrix.RData")
Genomes_Cluster_Run_4<-Read_Cluster(full_matrix, "Run_4")
Genomes_Cluster_Run_4

print("run_4")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_5/CI_Run_5_full_matrix.RData")
Genomes_Cluster_Run_5<-Read_Cluster(full_matrix, "Run_5")
Genomes_Cluster_Run_5

print("run_5")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_6/CI_Run_6_full_matrix.RData")
Genomes_Cluster_Run_6<-Read_Cluster(full_matrix, "Run_6")
Genomes_Cluster_Run_6

print("run_6")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_7/CI_Run_7_full_matrix.RData")
Genomes_Cluster_Run_7<-Read_Cluster(full_matrix, "Run_7")
Genomes_Cluster_Run_7

print("run_7")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_8/CI_Run_8_full_matrix.RData")
Genomes_Cluster_Run_8<-Read_Cluster(full_matrix, "Run_8")
Genomes_Cluster_Run_8

print("run_8")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_9/CI_Run_9_full_matrix.RData")
Genomes_Cluster_Run_9<-Read_Cluster(full_matrix, "Run_9")
Genomes_Cluster_Run_9

print("run_9")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_1.2/CI_Run_1.2_full_matrix.RData")
Genomes_Cluster_Run_1.2<-Read_Cluster(full_matrix, "Run_1.2")
Genomes_Cluster_Run_1.2

print("run_1.2")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_2.2/CI_Run_2.2_full_matrix.RData")
Genomes_Cluster_Run_2.2<-Read_Cluster(full_matrix, "Run_2.2")
Genomes_Cluster_Run_2.2

print("run_2.2")
load("/gale/netapp/home/jswift/analysis/CI_1001/1001_combo_indexing_run_3.2/CI_Run_3.2_full_matrix.RData")
Genomes_Cluster_Run_3.2<-Read_Cluster(full_matrix, "Run_3.2")
Genomes_Cluster_Run_3.2

print("run3.2")
load("/gale/ddn/ddn_neomorph/tjain/full_matrix_files/CI_Run_4.2_full_matrix.RData")
Genomes_Cluster_Run_4.2<-Read_Cluster(full_matrix,"Run_4.2")
Genomes_Cluster_Run_4.2

load("/gale/ddn/ddn_neomorph/tjain/full_matrix_files/CI_Run_5.2_full_matrix.RData")
Genomes_Cluster_Run_5.2<-Read_Cluster(full_matrix,"Run_5.2")
Genomes_Cluster_Run_5.2


load("/gale/ddn/ddn_neomorph/tjain/full_matrix_files/CI_Run_6.2_full_matrix.RData")
Genomes_Cluster_Run_6.2<-Read_Cluster(full_matrix,"Run_6.2")
Genomes_Cluster_Run_6.2

load("/gale/ddn/ddn_neomorph/tjain/full_matrix_files/CI_Run_7.2_full_matrix.RData")
Genomes_Cluster_Run_7.2<-Read_Cluster(full_matrix,"Run_7.2")
Genomes_Cluster_Run_7.2

load("/gale/ddn/ddn_neomorph/tjain/full_matrix_files/CI_Run_8.2_full_matrix.RData")
Genomes_Cluster_Run_8.2<-Read_Cluster(full_matrix,"Run_8.2")
Genomes_Cluster_Run_8.2
print("done with load")
# Merge data
Genomes_Cluster_merged <- merge(Genomes_Cluster_Run_1, y = c(Genomes_Cluster_Run_3,Genomes_Cluster_Run_4,Genomes_Cluster_Run_5,Genomes_Cluster_Run_6,
                                                             Genomes_Cluster_Run_7,Genomes_Cluster_Run_8,Genomes_Cluster_Run_9,Genomes_Cluster_Run_1.2,
                                                             Genomes_Cluster_Run_2.2, Genomes_Cluster_Run_3.2,Genomes_Cluster_Run_4.2,Genomes_Cluster_Run_5.2,
							     Genomes_Cluster_Run_6.2,Genomes_Cluster_Run_7.2,Genomes_Cluster_Run_8.2))
print("checkpoint_1")
Genomes_Cluster_merged <- NormalizeData(object = Genomes_Cluster_merged, verbose = FALSE)
Genomes_Cluster_merged <- FindVariableFeatures(object = Genomes_Cluster_merged , selection.method = "vst", nfeatures = 2000)
save(Genomes_Cluster_merged, file = "/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster_merged.RData")

# Subset by variable features
Variable_Features<-VariableFeatures(Genomes_Cluster_merged)
length(Variable_Features)
Genomes_Cluster_Run_1<-Genomes_Cluster_Run_1[Variable_Features,]
Genomes_Cluster_Run_3<-Genomes_Cluster_Run_3[Variable_Features,]
Genomes_Cluster_Run_4<-Genomes_Cluster_Run_4[Variable_Features,]
Genomes_Cluster_Run_5<-Genomes_Cluster_Run_5[Variable_Features,]
Genomes_Cluster_Run_6<-Genomes_Cluster_Run_6[Variable_Features,]
Genomes_Cluster_Run_7<-Genomes_Cluster_Run_7[Variable_Features,]
Genomes_Cluster_Run_8<-Genomes_Cluster_Run_8[Variable_Features,]
Genomes_Cluster_Run_9<-Genomes_Cluster_Run_9[Variable_Features,]
Genomes_Cluster_Run_1.2<-Genomes_Cluster_Run_1.2[Variable_Features,]
Genomes_Cluster_Run_2.2<-Genomes_Cluster_Run_2.2[Variable_Features,]
Genomes_Cluster_Run_3.2<-Genomes_Cluster_Run_3.2[Variable_Features,]
Genomes_Cluster_Run_4.2<-Genomes_Cluster_Run_4.2[Variable_Features,]
Genomes_Cluster_Run_5.2<-Genomes_Cluster_Run_5.2[Variable_Features,]
Genomes_Cluster_Run_6.2<-Genomes_Cluster_Run_6.2[Variable_Features,]
Genomes_Cluster_Run_7.2<-Genomes_Cluster_Run_7.2[Variable_Features,]
Genomes_Cluster_Run_8.2<-Genomes_Cluster_Run_8.2[Variable_Features,]
# Integrate
plant_anchors <- FindIntegrationAnchors(object.list = list(Genomes_Cluster_Run_1, Genomes_Cluster_Run_3, Genomes_Cluster_Run_4,
                                                           Genomes_Cluster_Run_5, Genomes_Cluster_Run_6, Genomes_Cluster_Run_7,
                                                           Genomes_Cluster_Run_8, Genomes_Cluster_Run_9, Genomes_Cluster_Run_1.2,
                                                           Genomes_Cluster_Run_2.2, Genomes_Cluster_Run_3.2,Genomes_Cluster_Run_4.2,
							   Genomes_Cluster_Run_5.2,Genomes_Cluster_Run_6.2,Genomes_Cluster_Run_7.2,Genomes_Cluster_Run_8.2), reference= c(10), dims = 1:30, reduction = "rpca")

remove(Genomes_Cluster_Run_1, Genomes_Cluster_Run_3, Genomes_Cluster_Run_4, Genomes_Cluster_Run_5, Genomes_Cluster_Run_6,Genomes_Cluster_Run_7,Genomes_Cluster_Run_8, Genomes_Cluster_Run_9,Genomes_Cluster_Run_1.2,Genomes_Cluster_Run_2.2,Genomes_Cluster_Run_3.2,Genomes_Cluster_Run_4.2,Genomes_Cluster_Run_5.2,Genomes_Cluster_Run_6.2,Genomes_Cluster_Run_7.2,Genomes_Cluster_Run_8.2)

Genomes_Cluster_integrated <- IntegrateData(anchorset = plant_anchors, dims = 1:30)
remove(plant_anchors)

# Clustering
DefaultAssay(object = Genomes_Cluster_integrated) <- "integrated"
Genomes_Cluster_integrated <- ScaleData(object = Genomes_Cluster_integrated, vars.to.regress = c('nCount_RNA'), verbose = TRUE)
Genomes_Cluster_integrated <- RunPCA(object = Genomes_Cluster_integrated, npcs = 30, verbose = FALSE)
Genomes_Cluster_integrated <- RunUMAP(object = Genomes_Cluster_integrated, reduction = "pca", dims = 1:30)
Genomes_Cluster_integrated <- FindNeighbors(object = Genomes_Cluster_integrated, reduction = "pca", dims = 1:30)
Genomes_Cluster_integrated <- FindClusters(Genomes_Cluster_integrated, resolution = c(0,0.05,0.1,0.2,0.5,0.6,0.75,1.0))
# if change here, change below!
Genomes_Cluster_integrated <- RunTSNE(object = Genomes_Cluster_integrated, reduction = "pca", dims = 1:30)
save(Genomes_Cluster_integrated, file = "L0_Genomes_Cluster_integrated.RData")

# Transfering UMAP and Meta Data
Genomes_Cluster_merged[['umap']] <- Genomes_Cluster_integrated[['umap']]
Genomes_Cluster_merged[['tsne']] <- Genomes_Cluster_integrated[['tsne']]
Genomes_Cluster_merged[['integrated']] <- Genomes_Cluster_integrated[['integrated']]
Genomes_Cluster_merged@meta.data$seurat_clusters <- Genomes_Cluster_integrated@meta.data$seurat_clusters
Genomes_Cluster_merged@meta.data$integrated_snn_res.0 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.05 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.05
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.1 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.1
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.2 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.2
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.5 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.5
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.6 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.6
Genomes_Cluster_merged@meta.data$integrated_snn_res.0.75 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.0.75
Genomes_Cluster_merged@meta.data$integrated_snn_res.1 <- Genomes_Cluster_integrated@meta.data$integrated_snn_res.1
Idents(Genomes_Cluster_merged)<-Genomes_Cluster_merged@meta.data$integrated_snn_res.1
table(Idents(Genomes_Cluster_merged))
Genomes_Cluster_merged<-Genomes_Cluster_merged[,WhichCells(Genomes_Cluster_merged, idents = c(0:19))]  # remove tiny clusters
Genomes_Cluster<-Genomes_Cluster_merged

# Save
jpeg('/gale/ddn/ddn_neomorph/tjain/cluster_new.jpg')
Idents(Genomes_Cluster)<- Genomes_Cluster@meta.data$integrated_snn_res.1
plot(DimPlot(object = Genomes_Cluster, reduction = "tsne", label = TRUE))
dev.off()
save(Genomes_Cluster, file = "/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster.RData")
Genomes_Cluster_subsample <- Genomes_Cluster[, sample(colnames(Genomes_Cluster), replace=FALSE)]
save(Genomes_Cluster_subsample, file = "/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster_subsample.RData")
