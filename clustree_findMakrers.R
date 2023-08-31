library(Seurat)
library(Matrix)
library(stringr)
library(mrtree)

args <- commandArgs(trailingOnly = TRUE)
#load Seurat object
load("/gale/ddn/ddn_neomorph/tjain/L0_Genomes_Cluster.RData")

setwd("/gale/ddn/ddn_neomorph/tjain/")
#making clustree plots
out = mrtree(Genomes_Cluster, prefix = 'integrated_snn_res.', n.cores = 2, consensus=FALSE, augment.path=FALSE)
jpeg('clustree.jpg')
plot_clustree(labelmat= Genomes_Cluster[[]], prefix ='integrated_snn_res.', plot.ref = FALSE)
jpeg('MRTree.jpg')
plot_tree(labelmat=out$labelmat.mrtree, plot.piechart = TRUE,
          node.size = 0.2, tip.label.dist = 10, bottom.margin=30)

#FindAllMarkers()

#available_res <- rownames(Genomes_Cluster[[]])[grep(paste0("^", "integrated_snn_res."), rownames(Genomes_Cluster[[]]))]
#max(as.numeric(sub(paste0("^",  "integrated_snn_res."), "", available_res)))
Idents(Genomes_Cluster) <- args[1]
Genomes_Cluster.markers <- FindAllMarkers(Genomes_Cluster, assay = 'RNA',features = NULL, logfc.threshold = 0.25, test.use="wilcox",slot="data",min.pct=0.1,min.diff.pct=-Inf,node = NULL, verbose = TRUE, only.pos = FALSE, max.cells.per.ident = Inf,random.seed = 1, latent.vars = NULL, min.cells.feature = 3, min.cells.group = 3, mean.fxn = NULL,fc.name = NULL, base = 2, return.thresh = 0.01,densify = FALSE)
write.csv(Genomes_Cluster.markers, "FindAllMarkers.csv", row.names=FALSE) 
