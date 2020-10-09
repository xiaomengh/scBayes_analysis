library(dplyr)
library(Seurat)
library(cowplot)
library(sctransform)
require(gridExtra)

sample = "Pt4_AGG"
cll = Read10X("Pt4_AGG/outs/filtered_feature_bc_matrix")
cll=CreateSeuratObject(counts = cll, project="cll", min.cells = 3, min.features = 200)

cll[['percent.mt']] = PercentageFeatureSet(object=cll, pattern="^MT-")
cll = subset(x=cll, subset = nFeature_RNA >200 & nFeature_RNA <4000 & percent.mt<15)

cll = SCTransform(cll)

cll = RunPCA(object = cll, features = VariableFeatures(object=cll))
DimPlot(object = cll, reduction="pca")
DimHeatmap(object = cll, dims=1:15, cells=500, balanced = TRUE)
cll = JackStraw(object=cll, num.replicate = 100)
cll = ScoreJackStraw(object = cll, dims = 1:20)
JackStrawPlot(object = cll, dims=1:15)

####### Cluster the cells #############
cll = RunUMAP(cll, dims = 1:30)
cll = FindNeighbors(object = cll, dims = 1:30)
cll = FindClusters(object = cll, resolution = 0.5)

umapplot=paste(sample,"_umap.pdf",sep="")
pdf(umapplot, width=8, height = 6)
DimPlot(object = cll, reduction = "umap", label=TRUE) + NoLegend()
dev.off()

markerplot = paste(sample,"_markers.pdf", sep="")
pdf(markerplot, width = 8, height = 10)
FeaturePlot(object = cll, features = c("MS4A1","CD19","GNLY", "CD3E", "CD8A","CD14"))
dev.off()

savefile = paste(sample,"_cll.rds", sep="")
saveRDS(cll, file=savefile)



