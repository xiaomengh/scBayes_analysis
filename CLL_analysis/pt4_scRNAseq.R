library(dplyr)
library(Seurat)
library(cowplot)
library(sctransform)
require(gridExtra)

####### pt4 for scBayes#######
sample="Pt4_AGG"
cll = readRDS("Pt4_AGG_cll.rds")
cll$cluster = Idents(cll)
bcell = subset(cll,idents=c("1", "2", "7","8"))

identfile = paste(sample, "_timepoints.csv",sep="")
tp = read.csv(identfile,header=F)
head(tp)
cll$tp = tp$V2

##### sample level expression ####
DimPlot(bcell, reduction="umap",split.by="tp")
bcell.tp = tp[tp$V1 %in% colnames(bcell),]
Idents(bcell) = bcell.tp$V2

diffgenes_tp_1=FindMarkers(bcell, ident.1="T1",ident.2="T2")
diffgenes_tp_2=FindMarkers(bcell, ident.1="T2",ident.2="T3")
diffgenes_tp_3=FindMarkers(bcell, ident.1="T1",ident.2="T3")

gene="CD69"
diffgenes_tp_1[gene,]
diffgenes_tp_2[gene,]
diffgenes_tp_3[gene,]
pdf(paste(sample,"_",gene,"_exp_tp.pdf",sep=""),width=4,height=3)
VlnPlot(bcell, features = gene,pt.size=0.1)
dev.off()

gene="TNFRSF13B"
diffgenes_tp_1[gene,]
diffgenes_tp_2[gene,]
diffgenes_tp_3[gene,]
pdf(paste(sample,"_",gene,"_exp_tp.pdf",sep=""),width=4,height=3)
VlnPlot(bcell, features = gene,pt.size=0.1)
dev.off()

gene="TXNIP"
diffgenes_tp_1[gene,]
diffgenes_tp_2[gene,]
diffgenes_tp_3[gene,]
pdf(paste(sample,"_",gene,"_exp_tp.pdf",sep=""),width=4,height=3)
VlnPlot(bcell, features = gene,pt.size=0.1)
dev.off()

assignment=read.csv("pt4_hg19_wes.assign.csv",header=F)
head(assignment)
rownames(assignment)=assignment$V1
bcell_assign= assignment[rownames(assignment) %in% colnames(bcell),]
bcell_tp = tp[tp$V1 %in% colnames(bcell),]
head(bcell_assign)
head(bcell_tp)
bcell$subclone=bcell_assign$V2
Idents(bcell)=bcell$subclone
pdf(paste(sample,"_assignment_byTP.pdf",sep=""), width = 18, height = 6)
DimPlot(object = bcell, reduction = "umap", label=F,split.by="tp",cols=c("#238443","#225ea8","#fe9929","#cccccc"))
dev.off()

############# subclone_tp phenotype analysis ###########
Idents(bcell)=paste(bcell_assign$V2,bcell_tp$V2,sep="_")
table(Idents(bcell))

subclone=subset(bcell,idents=c("SC1_T1","SC2_T1",
                               "SC1_T2","SC2_T2",
                               "SC1_T3","SC2_T3"))
reorder_levels <- c("SC1_T1","SC2_T1",
                    "SC1_T2","SC2_T2",
                    "SC1_T3","SC2_T3")
subclone@active.ident <- factor(x = subclone@active.ident, levels = reorder_levels)

diffgenes1=FindMarkers(subclone, ident.1="SC1_T1",ident.2="SC2_T1")
diffgenes2=FindMarkers(subclone, ident.1="SC1_T2",ident.2="SC2_T2")
diffgenes3=FindMarkers(subclone, ident.1="SC1_T3",ident.2="SC2_T3")
diffgenes4=FindMarkers(subclone, ident.1="SC1_T1",ident.2="SC1_T2")
diffgenes5=FindMarkers(subclone, ident.1="SC1_T2",ident.2="SC1_T3")
diffgenes6=FindMarkers(subclone, ident.1="SC1_T1",ident.2="SC1_T3")
diffgenes7=FindMarkers(subclone, ident.1="SC2_T1",ident.2="SC2_T2")
diffgenes8=FindMarkers(subclone, ident.1="SC2_T2",ident.2="SC2_T3")
diffgenes9=FindMarkers(subclone, ident.1="SC2_T1",ident.2="SC2_T3")

gene="TNFRSF13B"
diffgenes1[gene,]
diffgenes2[gene,]
diffgenes3[gene,]
diffgenes4[gene,]
diffgenes5[gene,]
diffgenes6[gene,]
diffgenes7[gene,]
diffgenes8[gene,]
diffgenes9[gene,]
pdf(paste(sample,"_",gene,"_subclone_exp_nonormal.pdf",sep=""),width=4,height=3)
VlnPlot(subclone, features = gene,pt.size=0.1)
dev.off()

gene="TXNIP"
diffgenes1[gene,]
diffgenes2[gene,]
diffgenes3[gene,]
diffgenes4[gene,]
diffgenes5[gene,]
diffgenes6[gene,]
diffgenes7[gene,]
diffgenes8[gene,]
diffgenes9[gene,]
pdf(paste(sample,"_",gene,"_subclone_exp_nonormal.pdf",sep=""),width=4,height=3)
VlnPlot(subclone, features = gene,pt.size=0.1)
dev.off()

gene="CD69"
diffgenes1[gene,]
diffgenes2[gene,]
diffgenes3[gene,]
diffgenes4[gene,]
diffgenes5[gene,]
diffgenes6[gene,]
diffgenes7[gene,]
diffgenes8[gene,]
diffgenes9[gene,]
pdf(paste(sample,"_",gene,"_subclone_exp_nonormal.pdf",sep=""),width=4,height=3)
VlnPlot(subclone, features = gene,pt.size=0.1)
dev.off()


