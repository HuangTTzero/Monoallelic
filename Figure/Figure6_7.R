# Analysis the Epha7 gene
# the promoter track plot
# Figure 6a
for(i in c("Epha7")){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/high/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=8,ncol=2,heights = c(1,1.2,1.2,1.2,0.4,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}

# result of Fish
# Figure 6b-6e
# for multiomic seq
Fish_gene_dscore<-dscore_data%>%subset(gene %in% c("Adcy3","Epha7","Peg3") & allele_exp>0&allele_frag>0&(allele_exp+allele_frag>=4)) 
Fish_gene_dscore<-Fish_gene_dscore%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
Fish_gene_count<-Fish_gene_dscore%>%group_by(gene)%>%summarise(count=n())
Fish_gene_celltype_count<-Fish_gene_dscore%>%group_by(gene,merge_dscore_type)%>%summarise(cellnumber=n())
Fish_gene_celltype_ratio<-merge(Fish_gene_celltype_count,Fish_gene_count,by=c("gene"))%>%mutate(ratio=cellnumber/count)
# for single-cell RNA-seq
Fish_gene_dscore_rna<-dscore_data%>%subset(gene %in% c("Adcy3","Epha7","Peg3") & allele_exp>=4) 
Fish_gene_dscore_rna<-Fish_gene_dscore_rna%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))
Fish_gene_count_rna<-Fish_gene_dscore_rna%>%group_by(gene)%>%summarise(count=n())
Fish_gene_celltype_count_rna<-Fish_gene_dscore_rna%>%group_by(gene,rna_dscore_type)%>%summarise(cellnumber=n())
Fish_gene_celltype_ratio_rna<-merge(Fish_gene_celltype_count_rna,Fish_gene_count_rna,by=c("gene"))%>%mutate(ratio=cellnumber/count)


Epha7<-data.frame(allele_exp_type=c("P","P","P","M","M","M","B","B","B","P","M","B","P","M","B"),rep=c("rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep1","rep1","rep1","rep1","rep1"),
                ratio=c(29.5,32.6,41.7,39.3,34.9,27.1,31.1,32.6,31.2,11.3,7.7,81,13.5,9.6,76.9),
                method=c(rep("RNA_Fish",9),rep("Multiomics-seq",3),rep("Single-cell-RNA-seq",3)))
Epha7$method<-factor(Epha7$method,levels=c("Multiomics-seq","Single-cell-RNA-seq","RNA_Fish"))
Epha7_p<-ggplot(Epha7,aes(x=method,y=ratio,fill=allele_exp_type,color=allele_exp_type))+
    geom_bar(color="black",stat="summary",fun=mean,position="dodge",size=0)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1, shape = 21, stroke = 0.2,color="black") +
    stat_summary(fun.data='mean_sd',geom="errorbar",linewidth=0.15,size=1,width = 0.5,position = position_dodge(0.9))+
    scale_fill_manual(values=c("#dfe2d6","#c68989","#94b2c2"))+
    scale_color_manual(values=c("black","black","black"))+
    scale_y_continuous(limits=c(0,100))+
    labs(y="% cell",
       x="")+
    theme_classic()+
    theme(
        legend.title = element_blank(),
        text=element_text(size=6)
    )
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/Fish/experimental/Epha7_Fishresult.pdf",Epha7_p, device = "pdf", width = 6, height = 2, units = "cm", dpi = 60,bg="white")

Peg3<-data.frame(allele_exp_type=c("P","P","P","M","M","M","B","B","B","P","M","B","P","M","B"),rep=c("rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep1","rep1","rep1","rep1","rep1"),
                ratio=c(90,92,88.2,10,0,0,0,8,11.8,93.3,0,6.7,75,0,25),
                method=c(rep("RNA_Fish",9),rep("Multiomics-seq",3),rep("Single-cell-RNA-seq",3)))
Peg3$method<-factor(Peg3$method,levels=c("Multiomics-seq","Single-cell-RNA-seq","RNA_Fish"))
Peg3_p<-ggplot(Peg3,aes(x=method,y=ratio,fill=allele_exp_type,color=allele_exp_type))+
    geom_bar(color="black",stat="summary",fun=mean,position="dodge",size=0)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1, shape = 21, stroke = 0.2, color="black") +
    stat_summary(fun.data='mean_sd',geom="errorbar",linewidth=0.15,size=1,width = 0.5,position = position_dodge(0.9))+
    scale_fill_manual(values=c("#dfe2d6","#c68989","#94b2c2"))+
    scale_color_manual(values=c("black","black","black"))+
    scale_y_continuous(limits=c(0,100))+
    labs(y="% cell",
       x="")+
    theme_classic()+
    theme(
        legend.title = element_blank(),
        text=element_text(size=6)
    )
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/Fish/experimental/Peg3_Fishresult.pdf",Peg3_p, device = "pdf", width = 6, height = 2, units = "cm", dpi = 60,bg="white")


Adcy3<-data.frame(allele_exp_type=c("P","P","P","M","M","M","B","B","B","P","M","B","P","M","B"),rep=c("rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep1","rep1","rep1","rep1","rep1"),
                ratio=c(0,0,0,0,0,2,100,100,98,4.1,6.6,89.3,6.7,12.3,81),
                method=c(rep("RNA_Fish",9),rep("Multiomics-seq",3),rep("Single-cell-RNA-seq",3)))
Adcy3$method<-factor(Adcy3$method,levels=c("Multiomics-seq","Single-cell-RNA-seq","RNA_Fish"))
Adcy3_p<-ggplot(Adcy3,aes(x=method,y=ratio,fill=allele_exp_type,color=allele_exp_type))+
    geom_bar(color="black",stat="summary",fun=mean,position="dodge",size=0)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1, shape = 21, stroke = 0.2, color="black") +
    stat_summary(fun.data='mean_sd',geom="errorbar",linewidth=0.15,size=1,width = 0.5,position = position_dodge(0.9))+
    scale_fill_manual(values=c("#dfe2d6","#c68989","#94b2c2"))+
    scale_color_manual(values=c("black","black","black"))+
    scale_y_continuous(limits=c(0,100))+
    labs(y="% cell",
       x="")+
    theme_classic()+
    theme(
        legend.title = element_blank(),
        text=element_text(size=6)
    )
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/Fish/experimental/Adcy3_Fishresult.pdf",Adcy3_p, device = "pdf", width = 6, height = 2, units = "cm", dpi = 60,bg="white")


# Identification of expressed OR genes in each mOSN
# load OR gene list, and remove PS gene
ORgene<-read.xlsx("/data/R02/huangtt39/data/OSN/12864_2020_6583_MOESM2_ESM.xlsx",sheet = 2)
ORgene<-unique(ORgene$Gene.symbol)
load(file="~/ATAC-RNAseq/F1joint_mask/Combined_mOSN.RData")
OR_exp<-Combined_mOSN[["RNA"]]@data[which(rownames(Combined_mOSN[["RNA"]]@data) %in% ORgene),]
myORgene<-rownames(OR_exp)

# remove cells express multiple OR genes
calculate_ORNs_percentage <- function(Olfr_counts,cutoff){
  if(cutoff==0){
    each_cell_expressed_Olfr_N <- rowSums(Olfr_counts>cutoff)
  } else {
    each_cell_expressed_Olfr_N <- rowSums(Olfr_counts>=cutoff)
  }
  zero_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==0)*100
  one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N==1)*100
  more_than_one_Olfr_percentage <- mean(each_cell_expressed_Olfr_N>1)*100
  df <- data.frame(cutoff=cutoff,Olfr_N=c("0","1",">1"),percentage=c(zero_Olfr_percentage,one_Olfr_percentage,more_than_one_Olfr_percentage))
  return(df)
}
cutoffs <- c(0:8,16,32,64,128,256)
newpalette <- c(brewer.pal(8,"Set2")[8],brewer.pal(9,"Blues")[c(5,3)])
pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/multiome_scRNA_mature_ORNs_UMIs_cutoff_tunning.pdf",width=8)
i <- "_1"
cells <- grep(i,colnames(Combined_mOSN),value=TRUE)
Olfr_counts <- FetchData(object =Combined_mOSN,vars=myORgene,slot = "counts",cells=cells)
df <- c()
for (j in 1:length(cutoffs)){
  df <- calculate_ORNs_percentage(Olfr_counts,cutoffs[j]) %>% rbind(df, .)
}
df$cutoff <- factor(df$cutoff,levels=cutoffs)
df$Olfr_N <- factor(df$Olfr_N,levels=c("0","1",">1"))
p <- ggplot(data=df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of OSNs",color="# of ORs",title="joint202201")+
  geom_vline(xintercept=2,linetype="longdash")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.8),face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
print(p)
i <- 2
cells <- grep(i,colnames(Combined_mOSN),value=TRUE)
Olfr_counts <- FetchData(object =Combined_mOSN,vars=myORgene,slot = "counts",cells=cells)
df <- c()
for (j in 1:length(cutoffs)){
  df <- calculate_ORNs_percentage(Olfr_counts,cutoffs[j]) %>% rbind(df, .)
}
df$cutoff <- factor(df$cutoff,levels=cutoffs)
df$Olfr_N <- factor(df$Olfr_N,levels=c("0","1",">1"))
p <- ggplot(data=df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of OSNs",color="# of ORs",title="jointF1OSN")+
  geom_vline(xintercept=2,linetype="longdash")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.8),face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
print(p)
i <- 3
cells <- grep(i,colnames(Combined_mOSN),value=TRUE)
Olfr_counts <- FetchData(object =Combined_mOSN,vars=myORgene,slot = "counts",cells=cells)
df <- c()
for (j in 1:length(cutoffs)){
  df <- calculate_ORNs_percentage(Olfr_counts,cutoffs[j]) %>% rbind(df, .)
}
df$cutoff <- factor(df$cutoff,levels=cutoffs)
df$Olfr_N <- factor(df$Olfr_N,levels=c("0","1",">1"))
p <- ggplot(data=df,aes(x=cutoff,y=percentage,color=Olfr_N,group=Olfr_N))+
  geom_point(size=2) +
  geom_line() +
  scale_color_manual(values=newpalette)+
  labs(x="Threshold(# of UMIs)",y="Percent of OSNs",color="# of ORs",title="jointF18w")+
  geom_vline(xintercept=2,linetype="longdash")+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5,size=rel(1.8),face="bold"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.title=element_text(size=rel(1.6)),axis.text.y = element_text(color="black",size=rel(1.5)), axis.text.x = element_text(color="black",size=rel(1.5)),axis.line = element_line(colour="black",size = 1))
print(p)
dev.off()

# dataset 1 : an OR was considered expressed in any given OSN if at least 1 UMIs were detected
# dataset 2 : an OR was considered expressed in any given OSN if at least 2 UMIs were detected
# dataset 3 : an OR was considered expressed in any given OSN if at least 3 UMIs were detected
# any OSN expressing either zero or multiple OR genes was not considered further for any downstream analyses

cutoffs <- c(1,2,3)
meta_df <- c()
samples<-c("_1","_2","_3")
for (i in 1:length(samples)){
  cells <- grep(samples[i],colnames(Combined_mOSN),value=TRUE)
  Olfr_counts <- FetchData(object = Combined_mOSN,vars=myORgene,slot = "counts",cells=cells)
  kept_cells <- rownames(Olfr_counts)[which(rowSums(Olfr_counts>=cutoffs[i])==1)]
  filtered_Olfr_counts <- Olfr_counts[kept_cells,]
  expressed_OR_idx <- apply(filtered_Olfr_counts,1,which.max)
  filtered_Olfr_counts$expressed_OR <- colnames(filtered_Olfr_counts)[expressed_OR_idx]
  meta_df <- filtered_Olfr_counts[,"expressed_OR",drop=FALSE] %>% rbind(meta_df,.)
}

filtered_mature_OSNs <- subset(Combined_mOSN,cells=rownames(meta_df))
filtered_mature_OSNs <- AddMetaData(filtered_mature_OSNs,meta_df)
DefaultAssay(filtered_mature_OSNs) <- "RNA"
filtered_mature_OSNs <- ScaleData(filtered_mature_OSNs)
save(filtered_mature_OSNs,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/filtered_mature_OSNs.RData")


# get OR subgroups with cells >=30
length(which(table(filtered_mature_OSNs$expressed_OR)>=30)) # 191 groups
OR_group_up<-filtered_mature_OSNs@meta.data%>%group_by(expressed_OR)%>%summarise(count=n())%>%subset(count>=30)
# the median and mean epression of Epha7 in different OR subgroups 
filtered_mature_OSNs_up<-subset(filtered_mature_OSNs,cells=rownames(filtered_mature_OSNs@meta.data%>%subset(expressed_OR%in%OR_group_up$expressed_OR)))
Epha7_median<-c()
for(OR in OR_group_up$expressed_OR){
      temp_barcode<-rownames(filtered_mature_OSNs_up@meta.data%>%subset(expressed_OR==OR))
      temp_median<-median(filtered_mature_OSNs_up[["RNA"]]@data["Epha7",temp_barcode])
      Epha7_median<-c(Epha7_median,temp_median)
}
Epha7_median<-data.frame(median=Epha7_median,expressed_OR=OR_group_up$expressed_OR)%>%arrange(median)
Epha7_OR_up_exp<-AverageExpression(
  filtered_mature_OSNs,
  assay="RNA",
  layer="data",features=c("Epha7"),
  group.by = "expressed_OR"
)
Epha7_OR_up_exp<-as.data.frame(Epha7_OR_up_exp$RNA)%>%pivot_longer(cols = everything(),names_to="expressed_OR",values_to="avg_exp")%>%arrange(avg_exp)

OR_group_up<-merge(Epha7_median,OR_group_up,by="expressed_OR")
OR_group_up<-merge(Epha7_OR_up_exp,OR_group_up,by="expressed_OR")%>%arrange(median)
OR_group_up<-OR_group_up%>%arrange(avg_exp,median)

filtered_mature_OSNs_up$expressed_OR<-factor(filtered_mature_OSNs_up$expressed_OR,levels=OR_group_up$expressed_OR)

# vlnplot of Epha7 expression of 30 sampled OR subgroups
# Figure 7a
set.seed(1223)
OR_sample<-sample(1:nrow(OR_group_up),30)
OR_group_up_30<-OR_group_up[OR_sample,]%>%arrange(avg_exp,median)
palette <- colorRampPalette(c("#DFDDE4","#8fb4be", "#d04c50"))
# generate 30 colors
colors <- palette(30)
filtered_mature_OSNs_30<-subset(filtered_mature_OSNs,cells=rownames(filtered_mature_OSNs@meta.data%>%subset(expressed_OR%in%OR_group_up_30$expressed_OR)))
filtered_mature_OSNs_30$expressed_OR<-factor(filtered_mature_OSNs_up$expressed_OR,levels=OR_group_up_30$expressed_OR)

Epha7_OR_sample30<-VlnPlot(
  object =filtered_mature_OSNs_30,
  features = c("Epha7"),
  group.by =  "expressed_OR",
  pt.size =0,slot = 'data')&geom_boxplot()&NoLegend()&scale_fill_manual(values=rep("#507b8a",30))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/Epha7_OR_sample30.pdf",Epha7_OR_sample30, device = "pdf", width = 40, height = 10, units = "cm", dpi = 300)
# barplot of Epha7 expression of OR subgroups
Epha7_OR_meanexp_p<-ggplot(OR_group_up, aes( x =  avg_exp))+
  geom_rect(aes(xmin=0,xmax=6,ymin=0,ymax=25),fill="#e6e8ed")+
  geom_rect(aes(xmin=6,xmax=Inf,ymin=0,ymax=25),fill="#a2a3b7")+
  geom_histogram(binwidth = 0.5, fill = "#799399", color = "white",linewidth=0.1,alpha=0.8)+  theme_classic()+
  labs(x="",y="# OR goups",
       fill="mean expression")+scale_y_continuous(limits=c(0,25))+geom_vline(xintercept= c(6),linewidth=0.25,linetype="dashed")+
  theme(
    text=element_text(size=8),
    plot.title = element_text(hjust = 0.5,vjust = 0.5), 
    legend.title=element_text(size=8), 
    legend.text=element_text(size=8)
  )
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/Epha7_OR_meanexp.pdf",Epha7_OR_meanexp_p, device = "pdf", width = 10, height = 5, units = "cm", dpi = 300,bg="white")
# the OR subgroups can be divided into two groups with average expression of Epha7
OR_group_up<-OR_group_up%>%mutate(avg_exp_type=case_when(avg_exp>=0&avg_exp<6~"low",avg_exp>=6~"high"))
OR_group_up$avg_exp_type<-factor(OR_group_up$avg_exp_type,levels=c("low","high"))
table(OR_group_up$avg_exp_type)
#   high    low 
#     37    154


# Import the coordinates of the projection of OR subgroups to the olfactory bulb. 
# Based on the coordinates of OR, determine the relationship between Epha7 and DV as well as AP.
# as x increases, it moves toward P, and as y increases, it moves toward D. The results of predict and slide-seq are combined.
OR_OB<-read.table(file="/data/R02/huangtt39/data/OSN/projection/map_654.txt")%>%rename_with(~c("expressed_OR","cell","type","x","y"),1:5)
OR_OB_slide<-OR_OB%>%subset(type=="slide-seq")
OR_OB <- OR_OB%>%subset(! expressed_OR%in%OR_OB_slide$expressed_OR)
OR_OB<-rbind(OR_OB,OR_OB_slide)
# merge
OR_group_up<-merge(OR_group_up,OR_OB%>%select(-cell),by="expressed_OR",all.x=TRUE)
OR_group_up$avg_exp_rank<-rank(OR_group_up$avg_exp)
OR_group_up$x_rank<-rank(OR_group_up$x)
OR_group_up$y_rank<-rank(OR_group_up$y)


# Projection map of OR in the olfactory bulb
# Figure 7b
library(ggExtra)
plot_features <- function(df, options = "viridis", directions = -1) {
  require(ggplot2)
  require(viridis)
  require(jpeg)
  require(ggpubr)
  img = readJPEG("/data/R02/huangtt39/data/OSN/projection/projection3.jpg")
  p1<-ggplot(df, aes_string("x", "y", colour = paste0("avg_exp")))+background_image(img)+geom_point()+
      xlim(c(0,2450))+ylim(c(-3130,300))+scale_color_gradientn(colors =c("#cdd2db","#cab1b7","#45486f","#2d3861"))+theme_classic()+labs(x="A-P",y="D-V")
  p2<-ggplot(df, aes_string("x", "y", colour = paste0("avg_exp_type")))+background_image(img)+geom_point()+
      xlim(c(0,2450))+ylim(c(-3130,300))+theme_classic()+labs(x="A-P",y="D-V")+
      scale_fill_manual(values =c("#cdd2db","#45486f"))+
      scale_color_manual(values =c("#cdd2db","#45486f"))
  p2<-ggMarginal(p2,groupColour=TRUE,groupFill=TRUE,type = "density")
  return(p1+p2)
}
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/Epha7_OR_projection.pdf",plot_features(OR_group_up)+plot_layout(ncol=2,width=c(1,1.5)), device = "pdf", width = 25, height = 10, units = "cm", dpi = 300)


# DV and AP scores of genes with two expression patterns at different expression levels
# Figure 7c
OR_group_up_with_OB$scale_x<-(OR_group_up_with_OB$x - min(OR_group_up_with_OB$x)) / (max(OR_group_up_with_OB$x) - min(OR_group_up_with_OB$x))
OR_group_up_with_OB$scale_y<-(OR_group_up_with_OB$y - min(OR_group_up_with_OB$y)) / (max(OR_group_up_with_OB$y) - min(OR_group_up_with_OB$y))
Epha7_DV_p<-ggviolin(OR_group_up_with_OB, x = "avg_exp_type", y = "scale_y", fill = "avg_exp_type", 
              palette=c("#cdd2db","#45486f"),
              ylim=c(0, 1),
              alpha = 1,width = 0.5,
              legend = "none",#去掉legend
              xlab="", ylab="D-V score",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
              stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("high", "low")),label.y = c(0.95))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
              
Epha7_AP_p<-ggviolin(OR_group_up_with_OB, x = "avg_exp_type", y = "scale_x", fill = "avg_exp_type", 
              palette=c("#cdd2db","#45486f"), 
              ylim=c(0, 1),
              alpha = 1,width = 0.5,
              legend = "none",#去掉legend
              xlab="", ylab="A-P score",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
              stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("high", "low")),label.y = c(0.95))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/Epha7_OR_projection_vln1.pdf",Epha7_DV_p+Epha7_AP_p, device = "pdf", width = 18, height = 10, units = "cm", dpi = 300)



# Epha7 alleic expression pattern of each OR subgroups
Epha7_dscore<-final_gene_dscore%>%subset(gene=="Epha7")%>%select(gene,barcode,merge_dscore_type)
Epha7_dscore<-merge(Epha7_dscore,filtered_mature_OSNs_up@meta.data%>%select(barcode,expressed_OR),by="barcode")
Epha7_dscore_in_ORgroup<-Epha7_dscore%>%group_by(expressed_OR,merge_dscore_type)%>%summarise(cell=n())
temp<-Epha7_dscore%>%group_by(expressed_OR)%>%summarise(total=n())
Epha7_dscore_in_ORgroup<-merge(Epha7_dscore_in_ORgroup,temp,by="expressed_OR",all.x=TRUE)
Epha7_dscore_in_ORgroup<-Epha7_dscore_in_ORgroup%>%pivot_wider(names_from="merge_dscore_type",values_from="cell",id_cols=c(expressed_OR,total),values_fill=0)
Epha7_dscore_in_ORgroup<-merge(OR_group_up,Epha7_dscore_in_ORgroup,by="expressed_OR",all.x=TRUE)

Epha7_dscore_in_ORgroup$total[is.na(Epha7_dscore_in_ORgroup$total)] <- 0
Epha7_dscore_in_ORgroup$Biallelic[is.na(Epha7_dscore_in_ORgroup$Biallelic)] <- 0
Epha7_dscore_in_ORgroup$Paternal_specific[is.na(Epha7_dscore_in_ORgroup$Paternal_specific)] <- 0
Epha7_dscore_in_ORgroup$Maternal_specific[is.na(Epha7_dscore_in_ORgroup$Maternal_specific)] <- 0

Epha7_dscore_in_ORgroup$Bi_rate<-Epha7_dscore_in_ORgroup$Biallelic/Epha7_dscore_in_ORgroup$count
Epha7_dscore_in_ORgroup$Mono_rate<-(Epha7_dscore_in_ORgroup$Paternal_specific+Epha7_dscore_in_ORgroup$Maternal_specific)/Epha7_dscore_in_ORgroup$count


Epha7_dscore_in_ORgroup$avg_exp_type<-factor(Epha7_dscore_in_ORgroup$avg_exp_type,levels=c("low","high"))
# Biallelic / Monoalleic foldchange
Epha7_dscore_in_ORgroup$fc<-(Epha7_dscore_in_ORgroup$Biallelic+0.1)/(Epha7_dscore_in_ORgroup$Paternal_specific+Epha7_dscore_in_ORgroup$Maternal_specific+0.1)
Epha7_dscore_in_ORgroup$fc[Epha7_dscore_in_ORgroup$total==0]<-NA
Epha7_dscore_in_ORgroup$fc1<-Epha7_dscore_in_ORgroup$fc
Epha7_dscore_in_ORgroup$fc[Epha7_dscore_in_ORgroup$fc>=10]<-10 # max=10


# Epha7 expression level VS Epha7 biallelic expression cell rate
Epha7_expVSBicell_p<-ggplot(Epha7_dscore_in_ORgroup,aes(x =avg_exp , y =Bi_rate)) +
  geom_point(size=1,alpha=0.8,color="#799399") +
  # scale_y_continuous(limit=c(-0.2,0.5))+
  # scale_x_continuous(limit=c(-0.2,0.5))+
  geom_smooth(method="lm",span=1, size=0.5,  se = TRUE)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
             formula = y~x, parse = TRUE, size = 2,color="red")+
  labs(y=str_c(" % Epha7 biallelic expression cell in OR group"),
       x="Epha7 mean expression in OR group")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=6))

# Epha7 expression level VS Epha7 mono expression cell rate
Epha7_expVSMonocell_p<-ggplot(Epha7_dscore_in_ORgroup,aes(x =avg_exp , y =Mono_rate)) +
  geom_point(size=1,alpha=0.8,color="#799399") +
  geom_smooth(method="lm",span=1, size=0.5,  se = TRUE)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
             formula = y~x, parse = TRUE, size = 2,color="red")+
  labs(y=str_c(" % Epha7 monoallelic expression cell in OR group"),
       x="Epha7 mean expression in OR group")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=6))
ggsave(filename =  "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/Epha7_expVScell.pdf",Epha7_expVSBicell_p+Epha7_expVSMonocell_p, device = "pdf", width =15, height = 8, units = "cm", dpi = 300,bg="white")

# Figure 4d,e
# selected OR subgroups with more than three Epha7 allelic informative cells
# compare the foldchange between high and low OR subgroups
ORgroup_Epha7_fc<-ggplot(Epha7_dscore_in_ORgroup%>%subset(total>=3),aes(x = avg_exp_type, y = fc,fill=avg_exp_type))+geom_boxplot()+
            stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("low","high")))+
            scale_fill_manual(values=c("#cdd2db","#45486f"))+
            NoLegend()+labs(x=NULL,y=NULL)+theme_classic()+
            theme(text=element_text(size=6))

# compare the biallelic expression rate between high and low OR subgroups
ORgroup_Epha7_Bicell<-ggviolin(Epha7_dscore_in_ORgroup, x = "avg_exp_type", y = "Bi_rate", fill = "avg_exp_type", 
              palette=c("#cdd2db","#45486f"), 
              alpha = 1,width = 0.5,
              trim=TRUE,
              legend = "none",#去掉legend
              xlab="", ylab="% Epha7 biallelic expression cell in each OR group",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
                  stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("high", "low")))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
# compare the monoallelic expression rate between high and low OR subgroups
ORgroup_Epha7_Monocell<-ggviolin(Epha7_dscore_in_ORgroup, x = "avg_exp_type", y = "Mono_rate", fill = "avg_exp_type", 
              palette=c("#cdd2db","#45486f"), 
              alpha = 1,width = 0.5,
              trim=TRUE,
              legend = "none",#去掉legend
              xlab="", ylab="% Epha7 monoallelic expression cell in each OR group",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
              stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("high", "low")))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure4/ORgroup_Epha7_allele.pdf",ORgroup_Epha7_Bicell+ORgroup_Epha7_Monocell+ORgroup_Epha7_fc, device = "pdf", width = 24, height = 7, units = "cm", dpi = 300)

