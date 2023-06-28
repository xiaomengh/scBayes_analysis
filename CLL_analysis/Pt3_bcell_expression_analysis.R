library(Seurat)
require(gridExtra)
library(ggplot2)
library(ggrepel)

sample="Pt3_AGG"
bcell = readRDS("Pt3_Bcells.rds")

##### UMAP plot #####
pdf(paste(sample,"_assignment_byTP.pdf",sep=""), width = 18, height = 6)
print(DimPlot(object = bcell, reduction = "umap", label=F,split.by="tp",cols=c("#238443","#225ea8","#fe9929","#cccccc")))
dev.off()

##### Differential expressed Genes #####
Idents(bcell) = bcell$tp
pdf(paste(sample,"_genes_exp_tp.pdf",sep=""),width=12,height=3)
print(VlnPlot(bcell, features = c("TNFRSF13B","TXNIP","CD69"),pt.size=0.1))
dev.off()

Idents(bcell) = bcell$subclone
subclone=subset(bcell,idents=c("SC1","SC2"))
Idents(subclone) = paste(subclone$subclone,subclone$tp,sep="_")
reorder_levels <- c("SC1_T1","SC2_T1",
                    "SC1_T2","SC2_T2",
                    "SC1_T3","SC2_T3")
subclone@active.ident <- factor(x = subclone@active.ident, levels = reorder_levels)
pdf(paste(sample,"_genes_exp_subclone_tp.pdf",sep=""),width=12,height=3)
print(VlnPlot(subclone, features = c("TNFRSF13B","TXNIP","CD69"),pt.size=0.1))
dev.off()

##### Differential expressed Genes between subclones#####
subclone=subset(bcell,idents=c("SC1","SC2","normal"))
Idents(subclone) = subclone$subclone
method="negbinom"
diffgenes_sc_1=FindMarkers(subclone, ident.1="SC1",ident.2="SC2",test.use=method)
diffgenes_sc_2=FindMarkers(subclone, ident.1="SC1",ident.2="normal",test.use=method)
diffgenes_sc_3=FindMarkers(subclone, ident.1="SC2",ident.2="normal",test.use=method)

VolcanoPlot=function(data,color_value=c("blue","red")){
  data$gene=rownames(data)
  data$diffexpressed = "no"
  data$diffexpressed[data$avg_logFC>0 & data$p_val_adj<0.05] = "up"
  data$diffexpressed[data$avg_logFC< 0 & data$p_val_adj<0.05] = "down"
  data$delabel=NA
  data$delabel[data$diffexpressed !="no"] = data$gene[data$diffexpressed !="no"]
  p=ggplot(data, aes(x=avg_logFC,y=-log10(p_val_adj),col=diffexpressed,label=delabel))+
    geom_point()+
    theme_minimal()+
    scale_color_manual(values=color_value)+
    geom_text_repel()+
    scale_x_continuous(limits=c(-1,2))
  print(p)
}
pdf(paste(sample,"_DE_SC1vsSC2_volcanoplot.pdf",sep=""),width=6,height=6)
VolcanoPlot(diffgenes_sc_1)
dev.off()
pdf(paste(sample,"_DE_SC1vsNormal_volcanoplot.pdf",sep=""),width=6,height=6)
VolcanoPlot(diffgenes_sc_2,color_value = c("blue","black","red"))
dev.off()
pdf(paste(sample,"_DE_SC2vsNormal_volcanoplot.pdf",sep=""),width=6,height=6)
VolcanoPlot(diffgenes_sc_3,color_value = c("blue","black","red"))
dev.off()

reorder_levels <- c("SC1","SC2","normal")
subclone@active.ident <- factor(x = subclone@active.ident, levels = reorder_levels)
pdf(paste(sample,"_DEgenes_vlnplot.pdf",sep=""),width=8, height=8)
print(VlnPlot(subclone, features=c("MIR155HG","ID3","RAC2",
                             "FCER2","MS4A1","CD22"),pt.size=0.1))
dev.off()

