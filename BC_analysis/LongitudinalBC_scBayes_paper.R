library(RColorBrewer)
library(reshape2)
library(ggplot2)


#### ssGSEAdata analysis for subclones#####
ssGSEAdata = read.table("dataforssGSEA.txt", header=T,sep='\t')
ssGSEAdata = ssGSEAdata[,c(1,3:ncol(ssGSEAdata))]
num_samples = ncol(ssGSEAdata)

pheno=read.table("phenotable.txt",sep='\t',header= TRUE)
head(pheno)
colnames(pheno) = c("cells","assigned.cluster","treatment")
rownames(pheno) = pheno$cells

###### compare two clusters ######
tumor_cells = which(!pheno$assigned.cluster %in% "normal") +1

ssGSEAdata$tumorAve = apply(ssGSEAdata[,tumor_cells],1,mean)
ssGSEAdata$tumorSd = apply(ssGSEAdata[,tumor_cells],1,sd)
zscore = ssGSEAdata
zscore[,c(2:num_samples)]=(ssGSEAdata[,c(2:num_samples)]-ssGSEAdata$tumorAve)/ssGSEAdata$tumorSd

plotViolinForSubclones=function(pathway) {
  pathwayData=melt(zscore[zscore$Name==pathway,c(1:num_samples)])
  #print(pathwayData)
  pathwayData$cluster = as.factor(pheno$assigned.cluster)
  p = ggplot(pathwayData, aes(cluster, value,color = cluster)) +
    geom_violin()+
    scale_x_discrete(limits=c("SC2", "SC4","SC3","SC5")) + 
    geom_jitter(shape=16, position=position_jitter(0.2))+
    scale_color_manual(values=c("#FFFFFF","#FFFFFF","#56B4E9","#E94F0D","#B266FF","#E69F00","#FFFFFF")) +
    ylim(-2,3)
    #+ scale_color_brewer(palette="Dark2") 
  p+theme(legend.position = "none")
  print(p)
}


pdf("HALLMARK_EMT_subclonal.pdf",width=5.25, height=5.4)
plotViolinForSubclones(pathway = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
dev.off()

#### ssGSEAdata analysis for pre and post sample #####
plotViolin=function(pathway) {
  pathwayData=melt(zscore[zscore$Name==pathway,tumor_cells])
  #print(pathwayData)
  pathwayData$cluster = as.factor(pheno[tumor_cells-1,]$treatment)
  p = ggplot(pathwayData, aes(cluster, value,color = cluster)) +
    geom_violin()+
    scale_x_discrete(limits=c("pre","post")) + 
    geom_jitter(shape=16, position=position_jitter(0.2))+
    scale_color_manual(values=c("#E94F0D","#56B4E9")) +
    ylim(-2,3)
  #+ scale_color_brewer(palette="Dark2")
  p+theme(legend.position = "none")
  print(p)
}
pdf("HALLMARK_EMT_pre_post.pdf",width=2.7,height=5.4)
plotViolin(pathway = "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
dev.off()
