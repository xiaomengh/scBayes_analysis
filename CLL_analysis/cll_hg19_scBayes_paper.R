require(reshape);
require(ggplot2)

plotCellVariantMatrix = function(data,mainlabel,cells,variant_range=NA){
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
  b = df_heatmap[df_heatmap$value==0,]
  c = df_heatmap[df_heatmap$value>0,]
  plot(b$variable, b$variant,pch=16, cex=1,col="#00FF0020",ylim=c(1,nrow(matrix)),xlim=c(1,ncol(matrix)),
       xlab="Cell index",ylab="Variant index",main=mainlabel)
  points(c$variable, c$variant,pch=16, cex=1,col="#0000FF50",ylim=c(1,nrow(matrix)),xlim=c(1,ncol(matrix)))
}
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
positiveCellCount = function(data, cells) {
  matrix = data[,colnames(data) %in% cells]
  print("total cell number")
  print(sum(colnames(data) %in% cells))
  print("positive cell number")
  print(sum(apply(matrix, 2, max)>0))
}
maxScAf = function(data, cells) {
  matrix = data[,colnames(data) %in% cells]
  a = apply(matrix,2,max)
  print(a[a>0])
}
maxScAfGenotype = function(data,genotype, cells) {
  matrix = data[,colnames(data) %in% cells]
  matrix_geno = genotype[,colnames(genotype) %in% cells]
  a = apply(matrix,2,max)>0
  b = apply(matrix[,a],1,max)>0
  print(matrix_geno[,a][b,])
}

sample = "T1"
data = read.table(paste(sample,".scaf.tsv",sep=""),header=T)
assignment = read.table(paste(sample,".forR.assign",sep=""),header = T)
rownames(assignment) = assignment$Barcode
head(assignment)
idents_all = read.csv(paste(sample,".celltype.csv",sep=""),header=F)
rownames(idents_all)= idents_all$V1
cells = colnames(data)[c(2:ncol(data))]
plotCellVariantMatrix(data,mainlabel="genotype for all cells",cells)

## genotyping B cells ##
par(mfrow=c(1,2))
cells = rownames(idents_all[idents_all$V2 == "B cells",] )
head(cells)
plotCellVariantMatrix(data,mainlabel="genotype for B cells",cells)

## genotyping assigned cells ##
par(mfrow=c(2,4))
cells_1 = rownames(idents_all[idents_all$V2 == "B cells",])
cells_2 = rownames(assignment[assignment$ASIG == "SC1",])
cells = cells_1[cells_1 %in% cells_2]
head(cells)
plotCellVariantMatrix(data,mainlabel="genotype for SC1, B cells",cells)
df1 = convertCellVariantMatrix(data,cells)
max(df1$variable)

cells_2 = rownames(assignment[assignment$ASIG == "SC2",])
cells = cells_1[cells_1 %in% cells_2]
head(cells)
plotCellVariantMatrix(data,mainlabel="genotype for SC2, B cells",cells)
df2 = convertCellVariantMatrix(data,cells)
df2$variable = df2$variable + max(df1$variable)
max(df2$variable)

cells_2 = rownames(assignment[assignment$ASIG == "normal",])
cells = cells_1[cells_1 %in% cells_2]
head(cells)
plotCellVariantMatrix(data,mainlabel="genotype for normal, B cells",cells)
positiveCellCount(data, cells)
df3 = convertCellVariantMatrix(data,cells)
df3$variable = df3$variable + max(df2$variable)
max(df3$variable)


cells_2 = rownames(assignment[assignment$ASIG == "UNASSIGN",])
cells = cells_1[cells_1 %in% cells_2]
head(cells)
plotCellVariantMatrix(data,mainlabel="genotype for UNASSIGN, B cells",cells)
positiveCellCount(data, cells)
df4 = convertCellVariantMatrix(data,cells)
df4$variable = df4$variable + max(df3$variable)
max(df4$variable)


##### combine plot ######
par(mfrow=c(1,1))
df_combined_1 = rbind(df1,df2,df3,df4)
b = df_combined_1[df_combined_1$value==0,]
c = df_combined_1[df_combined_1$value>0,]
pdf(paste(sample,"_genotype_Bcell.pdf",sep=""),width=8,height=4)
plot(b$variable, b$variant,pch=16, cex=0.5,col="#00FF0050",ylim=c(1,max(df_combined_1$variant)),xlim=c(1,max(df_combined_1$variable)),
     xlab="Cell index",ylab="Variant index", main="Genotype of B cells")
points(c$variable, c$variant,pch=16, cex=0.5,col="#0000FF50",ylim=c(1,max(df_combined_1$variant)),xlim=c(1,max(df_combined_1$variable)))
abline(v=max(df1$variable) + 0.5)
abline(v=max(df2$variable) + 0.5)
abline(v=max(df3$variable) + 0.5)
abline(h = 4.5)
dev.off()