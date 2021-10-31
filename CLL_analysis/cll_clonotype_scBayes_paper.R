library(dplyr)
library(ggplot2)
library(gridExtra)

ct_sc = read.csv("pt3_CT_SC.csv",header=T)
ct_sc_no_unassign=ct_sc[!ct_sc$SC=="UNASSIGN",]
count=count(ct_sc_no_unassign,CT)

plotClonotype = function(sample){
  ct_sc_sample = ct_sc_no_unassign[ct_sc_no_unassign$sampleID==sample,]
  ct_sc1 = ct_sc_sample[ct_sc_sample$SC == "SC1",]
  ct_sc1_count = count(ct_sc1,CT)
  ct_sc2 = ct_sc_sample[ct_sc_sample$SC == "SC2",]
  ct_sc2_count = count(ct_sc2,CT)
  ct_normal = ct_sc_sample[ct_sc_sample$SC == "normal",]
  ct_normal_count = count(ct_normal,CT)
  ###### sample level clonotype #####
  ct_sc_sample_count = count
  for (i in 1:nrow(count)){
    x = count(ct_sc_sample,CT)
    if (ct_sc_sample_count[i,1] %in% x$CT) {
      ct_sc_sample_count[i,2] = x[x$CT == ct_sc_sample_count[i,1],2]
    } else 
    {
      ct_sc_sample_count[i,2] = 0
    }
  }
  matrix_sample = ct_sc_sample_count
  matrix_sample$percentage = matrix_sample$n/sum(matrix_sample$n)
  p0=ggplot(data=matrix_sample,aes(x=CT,y=percentage,fill="n")) +
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values="#636363") +
    ylim(0,1) +
    ylab("percentage") +
    theme_minimal() +
    theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4)
    )
  ###### subclonal clonotype#########
  ct_sc1_count = count
  for (i in 1:nrow(count)){
    x = count(ct_sc1,CT)
    if (ct_sc1_count[i,1] %in% x$CT) {
      ct_sc1_count[i,2] = x[x$CT == ct_sc1_count[i,1],2]
    } else 
    {
      ct_sc1_count[i,2] = 0
    }
  }
  ct_sc2_count = count
  for (i in 1:nrow(count)){
    x = count(ct_sc2,CT)
    if (ct_sc2_count[i,1] %in% x$CT) {
      ct_sc2_count[i,2] = x[x$CT == ct_sc2_count[i,1],2]
    } else 
    {
      ct_sc2_count[i,2] = 0
    }
  }
  ct_normal_count = count
  for (i in 1:nrow(count)){
    x = count(ct_normal,CT)
    if (ct_normal_count[i,1] %in% x$CT) {
      ct_normal_count[i,2] = x[x$CT == ct_normal_count[i,1],2]
    } else 
    {
      ct_normal_count[i,2] = 0
    }
  }
  matrix=as.data.frame(cbind(ct_sc1_count$n,ct_sc2_count$n,ct_normal_count$n))
  colnames(matrix) = c("SC1","SC2","normal")
  matrix$clonotype = count$CT
  #### percentage ####
  matrix$sc1_percent=matrix$SC1/sum(matrix$SC1)
  matrix$sc2_percent=matrix$SC2/sum(matrix$SC2)
  matrix$normal_percent=matrix$normal/sum(matrix$normal)
  p1 = ggplot(data=matrix,aes(x=clonotype,y=sc1_percent,fill="SC1")) +
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values="#2b8cbe") +
    ylim(0,1) +
    theme_minimal() +
    theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4)
    )
  p2 = ggplot(data=matrix,aes(x=clonotype,y=sc2_percent,fill="SC2")) +
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values="#fd8d3c") +
    ylim(0,1) +
    theme_minimal() +
    theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4)
    )
  p3 = ggplot(data=matrix,aes(x=clonotype,y=normal_percent,fill="Normal")) +
    geom_bar(stat="identity",position=position_dodge())+
    scale_fill_manual(values="#78c679") +
    ylim(0,1) +
    theme_minimal() +
    theme(
      legend.position="none",
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=4)
    )
  print(grid.arrange(p0,p1, p2, p3,nrow = 1))
}

pdf("T1_subclonal_clonotype_percentage.pdf",width=8,height=3)
plotClonotype(2)
dev.off()

pdf("T2_subclonal_clonotype_percentage.pdf",width=8,height=3)
plotClonotype(1)
dev.off()

pdf("T3_subclonal_clonotype_percentage.pdf",width=8,height=3)
plotClonotype(3)
dev.off()
















