require(reshape)
require(ggplot2)


convertCellVariantMatrix = function(data,cells,variant_range=NA){
  if (!is.na(variant_range)) {
    matrix = data[variant_range,][,colnames(data) %in% cells]
    matrix = cbind(data[variant_range,]$variant,matrix)
  } else {
    matrix = data[,colnames(data) %in% cells]
    matrix = cbind(data$variant,matrix)
  }
  #print(head(matrix[,c(1:5)]))
  colnames(matrix) = c("variant", 1:(ncol(matrix)-1))
  matrix$variant = c(1:nrow(matrix))
  #print(head(matrix[,c(1:5)]))
  df_heatmap=melt(matrix,id.vars="variant")
  #print(head(df_heatmap))
  df_heatmap$variant = as.numeric(df_heatmap$variant)
  df_heatmap$variable = as.numeric(df_heatmap$variable)
  df_heatmap$value=as.numeric(df_heatmap$value)
  #print(head(df_heatmap))
  return(df_heatmap)
}

data=read.table("atac.scaf.tsv",sep="\t",header=T)
assignment=read.table("atac.forR.assign",header=T)
rownames(assignment) = assignment$Barcode
cells = rownames(assignment[assignment$ASIG == "tumor",])
df1 = convertCellVariantMatrix(data,cells)

cells = rownames(assignment[assignment$ASIG == "normal",])
df2 = convertCellVariantMatrix(data,cells)
df2$variable = df2$variable + max(df1$variable)

cells = rownames(assignment[assignment$ASIG == "UNASSIGN",])
df3 = convertCellVariantMatrix(data,cells)
df3$variable = df3$variable + max(df2$variable)


df_combined_1 = rbind(df1,df2,df3)
b = df_combined_1[df_combined_1$value==0,]
c = df_combined_1[df_combined_1$value>0,]
pdf(paste("atac_genotype.pdf",sep=""),width=8,height=4)
plot(b$variable, b$variant,pch=16, cex=0.5,col="#00FF0050",ylim=c(1,max(df_combined_1$variant)),xlim=c(1,max(df_combined_1$variable)),
     xlab="Cell index",ylab="Variant index", main="Genotype")
points(c$variable, c$variant,pch=16, cex=0.5,col="#FF000050",ylim=c(1,max(df_combined_1$variant)),xlim=c(1,max(df_combined_1$variable)))
abline(v=max(df1$variable) + 0.5)
abline(v=max(df2$variable) + 0.5)
dev.off()
