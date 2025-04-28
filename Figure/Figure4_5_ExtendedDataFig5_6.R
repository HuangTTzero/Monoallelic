# Figure 4a, the experimental result is by GraphPad Prism 10

# Select information of OR genes
OR_dscore<-dscore_data%>%subset(gene %in% myORgene & allele_exp>0&(allele_exp+allele_frag)>=2) 
OR_dscore<-OR_dscore%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
# Here, the data from single-cell nuclei are removed to better calculate the distance between OR groups, because fewer genes are detected in single-cell nuclei
OR_cell<-filtered_mature_OSNs@meta.data%>%subset(!orig.ident=="dataset1")
OR_cell<-data.frame(barcode=rownames(OR_cell),gene=OR_cell$expressed_OR)
OR_dscore<-merge(OR_cell,OR_dscore,by=c("gene","barcode"))
save(OR_dscore,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/OR_dscore.RData")

# add OR dscore to single cell object
add_dscore_df<-data.frame(barcode=rownames(filtered_mature_OSNs@meta.data),gene=filtered_mature_OSNs$expressed_OR)
add_dscore_df<-merge(add_dscore_df,OR_dscore%>%dplyr::select(mergescore,merge_dscore_type,barcode,gene),by=c("barcode","gene"),all.x=TRUE)
rownames(add_dscore_df)<-add_dscore_df$barcode
filtered_mature_OSNs$expressed_OR_dscore<-add_dscore_df[rownames(filtered_mature_OSNs@meta.data),3]
filtered_mature_OSNs$expressed_OR_dscore_type<-add_dscore_df[rownames(filtered_mature_OSNs@meta.data),4]

OR_celltype_count<-OR_dscore%>%group_by(gene,merge_dscore_type)%>%summarise(cellnumber=n())
temp_OR_count<-OR_dscore%>%group_by(gene)%>%summarise(count=n())
OR_celltype_count<-merge(OR_celltype_count,temp_OR_count,by="gene")
OR_count_ratio<-OR_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(gene,count),values_fill=0)
OR_count<-OR_celltype_count%>%mutate(ratio=cellnumber)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(gene,count),values_fill=0)

# expression pattern of OR genes
# Extended Data Fig.5a
ORgene_point_p <- ggplot(OR_count) +
      geom_point(aes(x = count, y = log2((Paternal_specific+0.1)/(Maternal_specific+0.1)), color =1-Biallelic/count),size=1) +  # alpha 调整点的透明度；shape 调整点的形状
      theme_classic() + scale_color_gradientn(colors =c("#fdfdfd","#f0eef6","#dadbe6","#bdbbdb","#3b037f"))+
      scale_y_continuous(breaks=seq(-10,10,2))+scale_x_continuous(expand=c(0.05,1),breaks=seq(0,150,25))+
      theme(legend.position = "right",legend.key.size = unit(8, "pt"),legend.text=element_text(size=6),text = element_text(size = 6)) + 
      labs(x = "allelic informatic cell number", y = "log2FC (paternal / maternal specific expressed cell)",color="monoallelic\nexpressed\ncell ratio")
ORgene_cellnumber_p <-ggplot(data=OR_count,aes(x=count))+
      geom_histogram(bins = 20,
                 color="white",
                 fill="grey",alpha=0.5)+
      scale_x_continuous(expand=c(0.025,0),breaks = seq(0,150,25))+
      theme_classic()+theme(text=element_text(size=6))+
      labs(y="gene number",
          x="allelic informatic cell number")
ORgene_FC_p <-ggplot(data=OR_count,aes(x=log2((Paternal_specific+0.1)/(Maternal_specific+0.1))))+
      geom_histogram(bins = 20,
                 color="white",
                 fill="grey",alpha=0.5)+
      scale_x_continuous(expand=c(0.04,0),breaks = seq(-10,10,2))+
      theme_classic()+theme(text=element_text(size=6))+
      labs(y="gene number",
          x="log2FC (paternal / maternal specific expressed cell)",color="monoallelic\nexpressed\ncell ratio")+
      coord_flip()
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure3/OR_cutoff.pdf",ORgene_cellnumber_p+plot_spacer()+ORgene_point_p+ORgene_FC_p+plot_layout(nrow=2,ncol=2,height=c(1,2,1,2),width=c(2,1,2,1)), device = "pdf", width =12, height = 9, units = "cm", dpi = 60,bg="white")


# Downsample to determine the cell number at which transcriptional differences are not affected by the cell count.
library(lsa)
OR_count_30<-OR_count%>%subset(Maternal_specific>=15&Paternal_specific>=15)
cellnumber<-c(10,20,25,30)
set.seed(1234)
get_sample_distance<-function(OR,n){
  OR_df <- OR_dscore%>%subset(gene==OR)
  g1_cells <- OR_df%>%subset(merge_dscore_type=="Paternal_specific")
  g1_cells<-g1_cells$barcode
  g1_cells<-sample(g1_cells,n)
  g2_cells <- OR_df%>%subset(merge_dscore_type=="Maternal_specific")
  g2_cells<-g2_cells$barcode
  g2_cells<-sample(g2_cells,n)
  PCA_df <- filtered_mature_OSNs[["pca"]]@cell.embeddings
  PCA_df <- t(PCA_df)
  PCA_df <- PCA_df[1:20,c(g1_cells,g2_cells)]
  cosine_similarity_df <- cosine(PCA_df)
  cosine_distance_df <- 1-cosine_similarity_df
  within_g1_dis <- cosine_distance_df[g1_cells,g1_cells]
  within_g1_dis <- unlist(lapply(2:nrow(within_g1_dis),function(i){
    within_g1_dis[i,1:(i-1)]
  }))
  within_g2_dis <- cosine_distance_df[g2_cells,g2_cells]
  within_g2_dis <- unlist(sapply(2:nrow(within_g2_dis),function(i){
    within_g2_dis[i,1:(i-1)]
  }))
  within_dis <- c(within_g1_dis,within_g2_dis)
  between_dis <- as.vector(cosine_distance_df[g1_cells,g2_cells])
  return(data.frame(gene=OR,cellnumber=n*2,within_median=median(within_dis),between_median=median(between_dis)))
}
sample_distance<-data.frame()
for(i in 1:nrow(OR_count_30)){
  for(number in cellnumber ){
    distance<-get_sample_distance(OR_count_30$gene[i],number/2)
    sample_distance<-rbind(sample_distance,distance)
  }
}
# Calculate the distance of OR gene expression between different parental cells, and sample to the same number of cells.
get_distance<-function(OR){
  OR_df <- OR_dscore%>%subset(gene==OR)
  g1_cells <- OR_df%>%subset(merge_dscore_type=="Paternal_specific")
  g1_cells<-g1_cells$barcode
  g2_cells <- OR_df%>%subset(merge_dscore_type=="Maternal_specific")
  g2_cells<-g2_cells$barcode
  cell<-min(length(g1_cells),length(g2_cells))
  g1_cells<-sample(g1_cells,cell)
  g2_cells<-sample(g2_cells,cell)
  PCA_df <- filtered_mature_OSNs[["pca"]]@cell.embeddings
  PCA_df <- t(PCA_df)
  PCA_df <- PCA_df[1:20,c(g1_cells,g2_cells)]
  cosine_similarity_df <- cosine(PCA_df)
  cosine_distance_df <- 1-cosine_similarity_df
  within_g1_dis <- cosine_distance_df[g1_cells,g1_cells]
  within_g1_dis <- unlist(lapply(2:nrow(within_g1_dis),function(i){
    within_g1_dis[i,1:(i-1)]
  }))
  within_g2_dis <- cosine_distance_df[g2_cells,g2_cells]
  within_g2_dis <- unlist(sapply(2:nrow(within_g2_dis),function(i){
    within_g2_dis[i,1:(i-1)]
  }))
  within_dis <- c(within_g1_dis,within_g2_dis)
  between_dis <- as.vector(cosine_distance_df[g1_cells,g2_cells])
  test<-wilcox.test(within_dis,between_dis)
  return(data.frame(gene=OR,cellnumber=2*cell,within_median=median(within_dis),between_median=median(between_dis),wilcox_p=test$p.value))
}
sample_ori_distance<-data.frame()
for(i in 1:nrow(OR_count_30)){
  distance<-get_distance(OR_count_30$gene[i])
  sample_ori_distance<-rbind(sample_ori_distance,distance)
}


sample_distance<-merge(sample_distance,sample_ori_distance%>%select(gene,within_median,between_median),by="gene")
library("gg.gap")
# Downsampling vs. original cell number 
# # Extended Data Fig.5b
rm(sample_distance_point_p)
for(i in c(10,20,25,30)){
  temp_p<-sample_distance %>% subset(cellnumber==i)%>%ggplot(aes(x =between_median.y-within_median.y , y =between_median.x-within_median.x)) +
  geom_point(size=1,alpha=0.8,color="#799399") +
  scale_y_continuous(limit=c(-0.2,0.5))+
  scale_x_continuous(limit=c(-0.2,0.5))+
  geom_smooth(method="lm",span=1, size=0.5,  se = TRUE)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
             formula = y~x, parse = TRUE, size = 2,color="red")+
  labs(y=str_c("Expression difference between intra and inter group with ",i," cells"),
       x="Expression difference between intra and inter group")+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=6))
  if(exists("sample_distance_point_p")){
    sample_distance_point_p<-sample_distance_point_p+temp_p
  }else{
    sample_distance_point_p<-temp_p
  }
}
ggsave(filename =  "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure3/OR_sample_distance_point.pdf",sample_distance_point_p+plot_layout(ncol=4,guides = 'collect'), device = "pdf", width =32, height = 8, units = "cm", dpi = 300,bg="white")


# finally get the OR subgroups with more than 20 cells 
OR_final<-OR_count%>%subset(count>=20&Maternal_specific>=10&Paternal_specific>=10)### 75 genes
# calculate the difference
set.seed(500)
distance<-data.frame()
for(i in 1:nrow(OR_final)){
  temp_distance<-get_distance(OR_final$gene[i])
  distance<-rbind(distance,temp_distance)
}
distance$difference<-distance$between_median-distance$within_median
distance<-distance%>%mutate(wilcox.test=case_when(wilcox_p<0.05~"p.value<0.05",wilcox_p>=0.05~"p.value>=0.05"))

# load the snp information
# load pwk specific snp
mm10_129_PWK_SNP <- read.table("/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/pwk_specific.vcf",header=FALSE)[,c(2:6)]%>%rename_with(~c("chr","pos","id","ref","alt"),1:5)
mm10_129_PWK_SNP$chr<-paste0("chr",mm10_129_PWK_SNP$chr)
# load ccds of each gene 
# use consensus CDS(ccds),0-based
ccdsGene_df <- read.table("/public/home/shipy3/DB/mm10/annotation/mm10_ccdsGene.txt")[,1:9]
colnames(ccdsGene_df) <- c("bin","ccdsId","chrom","strand","txStart","txEnd","cdsStart","cdsEnd","exonCount")
ccdsKgMap_df <- read.table("/public/home/shipy3/DB/mm10/annotation/mm10_ccdsKgMap.txt", header = FALSE, sep = "\t", stringsAsFactors = FALSE)[,c(1:2,6)]
colnames(ccdsKgMap_df) <- c("ccdsId","geneId","cdsSimilarity")
ccdsKgMap_df <- ccdsKgMap_df[which(ccdsKgMap_df$cdsSimilarity==1),]
CCDS_df <- merge(ccdsGene_df,ccdsKgMap_df,by="ccdsId")
CCDS_df$geneId <- gsub("\\.\\d","",CCDS_df$geneId)
txid2gene <- read.table("/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/genes/txid_to_name.txt")
colnames(txid2gene )<-c("geneid","gene")
merged_CCDS_df <- merge(CCDS_df,txid2gene,by.x="geneId",by.y="geneid",all.x=TRUE)
avail_OR_CCDS <- merged_CCDS_df%>%subset(gene%in%myORgene)%>% filter(exonCount==1)
# the ccds may have more than one txid, save the first one
write.table(avail_OR_CCDS%>%group_by(gene)%>%dplyr::slice(1),file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/OR_CCDS_txid.txt",col.names = FALSE, row.names = FALSE,quote = FALSE,sep = "\t")

# select the ccds of each OR genes ,and remove the duplicates
avail_OR_CCDS <- unique(avail_OR_CCDS[,c("ccdsId","chrom","cdsStart","cdsEnd","gene")])
# the number of concensus CDS for each OR gene
OR_CCDS_num_df <- avail_OR_CCDS %>% subset(gene%in%OR_final$gene)%>%
  group_by(gene) %>% 
  summarise(counts=n())
# all genes have one ccds
# caculate the snp count
final_OR_CCDS<-avail_OR_CCDS %>% subset(gene%in%OR_final$gene)
SNP_counts <- sapply(1:nrow(final_OR_CCDS),function(i){
  counts <- 0
  chrom <- final_OR_CCDS$chrom[i]
  start <-  final_OR_CCDS$cdsStart[i]
  end <- final_OR_CCDS$cdsEnd[i]
  OR_mm10_129_PWK_SNP <- mm10_129_PWK_SNP%>%subset(chr==chrom&pos>=start+1&pos<=end)
  nrow(OR_mm10_129_PWK_SNP)
}) 
final_OR_CCDS$SNP_counts <- SNP_counts

OR_SNP <- lapply(1:nrow(final_OR_CCDS),function(i){
  counts <- 0
  chrom <- final_OR_CCDS$chrom[i]
  start <-  final_OR_CCDS$cdsStart[i]
  end <- final_OR_CCDS$cdsEnd[i]
  OR_mm10_129_PWK_SNP <- mm10_129_PWK_SNP%>%subset(chr==chrom&pos>=start+1&pos<=end)
  OR_mm10_129_PWK_SNP$gene<-rep(final_OR_CCDS$gene[i],nrow(OR_mm10_129_PWK_SNP))
  OR_mm10_129_PWK_SNP<-merge(final_OR_CCDS,OR_mm10_129_PWK_SNP,by="gene")
}) 
OR_SNP<-do.call(rbind,OR_SNP)

# load the Non-synonymous mutation.
unsyn_variant<-read.table("/data/R02/huangtt39/ATAC-RNAseq/snp/snp_eff/missense/pwk_specific_missense.txt")[,c(1,4)]
colnames(unsyn_variant)<-c("id","geneId")
unsyn_variant$geneId <- gsub("\\.\\d","",unsyn_variant$geneId)

# Determine whether the site is a synonymous or non-synonymous mutation
# Retrieve the txid for each ccdsID again.
OR_unsyn_variant<-merge(merged_CCDS_df%>%select(geneId,ccdsId),OR_SNP,by=c("ccdsId"))
# Link SNPs using txid.
OR_unsyn_variant<-merge(OR_unsyn_variant,unsyn_variant,by=c("geneId","id"))
OR_unsyn_variant_uniq<-OR_unsyn_variant%>%select(id,chr,pos,gene,ccdsId)%>%unique()%>%group_by(ccdsId)%>%summarise(unsyn_snp_counts=n())
final_OR_CCDS<-merge(final_OR_CCDS,OR_unsyn_variant_uniq,by="ccdsId",all.x=TRUE)
final_OR_CCDS[is.na(final_OR_CCDS)]<-0

# merge the distance and snp information
OR_distance<-merge(distance,final_OR_CCDS,by="gene",all.x=TRUE)
write.xlsx(OR_distance,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_distance_SNP.xlsx")
OR_distance<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_distance_SNP.xlsx")

# Figure 4a
# ranking by OR difference
OR_distance<-OR_distance%>%arrange(between_median-within_median)
OR_distance$gene<-factor(OR_distance$gene,level=OR_distance$gene)
OR_difference_ranking<-ggplot(OR_distance)+
    geom_point(aes(x=gene,y=between_median-within_median,size=cellnumber,color=wilcox.test),alpha=0.7)+
    scale_y_break(breaks=c(0.3,0.5),ticklabels=seq(0.5,1,0.1),scales=0.3)+
    scale_y_continuous(limit=c(-0.08,0.7))+
    geom_text(data=OR_distance%>%subset(wilcox.test=="p.value<0.05"), aes(x =gene , y =between_median-within_median-0.04,label = gene), vjust =-0.5, size = 2.5)+
    scale_color_manual(values=c("#587880","grey"))+
      theme_classic()+
    #scale_y_break(c(-1000, 1000))+theme_prism(palette = "pearl",
    #           base_size = 8,
    #           base_line_size = 0.2,axis_text_angle = 45)+
    # scale_size_continuous(range = c(1,3),breaks=seq(-0.5,0,0.25),labels=abs(seq(-0.5,0,0.25)))+
   labs(y="Transcriptome difference bewteen C57 and PWK group",x="Gene ranking by difference")+
    theme(text=element_text(size=10),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
      axis.text.y=element_text(size=10,color = "black"),
        axis.ticks.x = element_blank(), # 去掉 x 轴的刻度
        axis.text.x = element_blank()   # 去掉 x 轴的文字
        )    
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_difference_ranking.pdf",OR_difference_ranking, device = "pdf", width=16, height = 8, units = "cm", dpi = 300,bg="white")

# Plot a heatmap and boxplot based on genes with higher difference
# Figure 4c
library(ComplexHeatmap)
library(circlize)
distance_heatmap<-function(OR){
  col_fun <- colorRamp2(
    c(1.5,1,0.5,0), 
    # c("#ccd4df","#e3c5c4","#914350", "#561319")
    c("#ccd4df","#e3c5c4","#625f86","#1c2d58")
  )
  OR_df <- OR_dscore%>%subset(gene==OR&merge_dscore_type!="Biallelic")
  PCA_df <- filtered_mature_OSNs[["pca"]]@cell.embeddings
  PCA_df <- t(PCA_df)
  PCA_df <- PCA_df[1:20,OR_df$barcode]
  cosine_similarity_df <- cosine(PCA_df)
  cosine_distance_df <- 1-cosine_similarity_df
  ha2 <- HeatmapAnnotation(
    df = OR_df[,15,drop=FALSE]
    col = list(merge_dscore_type=c("Paternal_specific"="#94b2c2", "Biallelic"="#dfe2d6", "Maternal_specific"="#c68989"))
    )
  h <- Heatmap(cosine_distance_df,top_annotation=ha2,name="Transcriptome distance",col=col_fun,width = unit(11, "cm"),height = unit(11, "cm"),show_row_names=FALSE,show_column_names=FALSE,border_gp = gpar(col = "black"))
  draw(h,heatmap_legend_side = "right",annotation_legend_side = "right")
}
pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_distance_heatmap.pdf",height=6.5,width=8.5)
distance_heatmap("Olfr1209")
distance_heatmap("Olfr1484")
distance_heatmap("Olfr1423")
distance_heatmap("Olfr1206")
dev.off()

distance_violin<-function(OR){
  OR_df <- OR_dscore%>%subset(gene==OR)
  g1_cells <- OR_df%>%subset(merge_dscore_type=="Paternal_specific")
  g1_cells<-g1_cells$barcode
  g2_cells <- OR_df%>%subset(merge_dscore_type=="Maternal_specific")
  g2_cells<-g2_cells$barcode
  cell<-min(length(g1_cells),length(g2_cells))
  set.seed(1222)
  g1_cells<-sample(g1_cells,cell)
  g2_cells<-sample(g2_cells,cell)
  PCA_df <- filtered_mature_OSNs[["pca"]]@cell.embeddings
  PCA_df <- t(PCA_df)
  PCA_df <- PCA_df[1:20,c(g1_cells,g2_cells)]
  cosine_similarity_df <- cosine(PCA_df)
  cosine_distance_df <- 1-cosine_similarity_df
  within_g1_dis <- cosine_distance_df[g1_cells,g1_cells]
  within_g1_dis <- unlist(lapply(2:nrow(within_g1_dis),function(i){
    within_g1_dis[i,1:(i-1)]
  }))
  within_g2_dis <- cosine_distance_df[g2_cells,g2_cells]
  within_g2_dis <- unlist(sapply(2:nrow(within_g2_dis),function(i){
    within_g2_dis[i,1:(i-1)]
  }))
  within_dis <- c(within_g1_dis,within_g2_dis)
  between_dis <- as.vector(cosine_distance_df[g1_cells,g2_cells])
  ##permutation distance
  perm_cells<-sample(OR_df$barcode,cell)
  PCA_df <- filtered_mature_OSNs[["pca"]]@cell.embeddings
  PCA_df <- t(PCA_df)
  PCA_df <- PCA_df[1:20,perm_cells]
  cosine_similarity_df <- cosine(PCA_df)
  cosine_distance_df <- 1-cosine_similarity_df
  perm_dis <- unlist(lapply(2:nrow(cosine_distance_df),function(i){
    cosine_distance_df[i,1:(i-1)]
  }))
  plot_df <- data.frame(distance=c(within_dis,between_dis,perm_dis),group=c(rep("within",length(within_dis)),rep("between",length(between_dis)),rep("permutation",length(perm_dis))))
  boxplot_p<-ggviolin(plot_df, x = "group", y = "distance", fill = "#766E8E", 
              alpha = 1,width = 0.5,ylim=c(0, 2.5),
              legend = "none",#去掉legend
              xlab="", ylab="Transcriptome distance",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
              stat_compare_means(method = "wilcox.test", label = "p.signif",comparisons = list(c("within", "between"),c("within", "permutation"),c("between", "permutation")),label.y = c(1.9,2.1,2.3))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=6))
  print(boxplot_p)
}
pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_distance_violin.pdf",height=2.5,width=3.5)
distance_violin("Olfr1209")
distance_violin("Olfr1484")
distance_violin("Olfr1423")
distance_violin("Olfr1206")
dev.off()


myOR_far<-c("Olfr1209","Olfr1484","Olfr1423","Olfr1206")
# subset the cell of Olfr1209, Olfr1484 and Olfr1206, tp plot the umap
OR_high_distance_cell<-OR_dscore%>%subset(gene%in%myOR_far[1:3])
OR_high_distance <- subset(filtered_mature_OSNs, subset=expressed_OR%in%myOR_far[1:3])

#Gene expression data processing
DefaultAssay(OR_high_distance) <- "inteRNA"
OR_high_distance <- ScaleData(OR_high_distance)
OR_high_distance <- RunPCA(OR_high_distance)
OR_high_distance <- RunUMAP(OR_high_distance,dims = 1:5, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

OR_high_distance_umap1 <- DimPlot(OR_high_distance, reduction = "umap.rna",label = TRUE, label.size = 2, group.by="expressed_OR",repel = TRUE, pt.size = 0.3) + ggtitle("RNA")
OR_high_distance_umap2 <- DimPlot(OR_high_distance, reduction = "umap.rna",label = TRUE, label.size = 2, group.by="expressed_OR_dscore_type",repel = TRUE, pt.size = 0.3) + ggtitle("RNA")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_high_distance_umap.pdf",OR_high_distance_umap1+OR_high_distance_umap2 & NoLegend() & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width =20, height = 10, units = "cm", dpi = 300)


# Differential analysis on the parental groups of high transcriptional difference genes Olfr1209, Olfr1484, and Olfr1423, as well as the negative control Olfr1206.
OR_obj<-lapply(1:4,function(i){
  OR_cell<-OR_dscore%>%subset(gene==myOR_far[i]&merge_dscore_type!="Biallelic")
  temp_OR_obj<-subset(filtered_mature_OSNs,cells=OR_cell$barcode)
  temp_OR_obj$OR_allele_exptype<-OR_cell$merge_dscore_type
  temp_OR_obj
})

# plot the umap
neighbor<-c(3,4,3,3)
dim<-c(5,5,10,10)
umap_list<-list()
for(i in 1:4){
  obj_temp<-OR_obj[[i]]
  DefaultAssay(obj_temp) <- "inteRNA"
  obj_temp <- ScaleData(obj_temp)
  obj_temp <- RunPCA(obj_temp,npcs=20)
  obj_temp <- RunUMAP(obj_temp,dims = 1:dim[i], reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_',n.neighbors=neighbor[i])
  umap_list[[i]] <- DimPlot(obj_temp, reduction = "umap.rna",label = FALSE, pt.size = 2, group.by="OR_allele_exptype",repel = TRUE,cols=c("#C68989","#94B2C2")) + ggtitle("RNA")
}
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/obj_temp_umap.pdf",umap_list[[1]]+umap_list[[2]]+umap_list[[3]]+umap_list[[4]]+plot_layout(ncol=2) & NoLegend() & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width =16, height = 16, units = "cm", dpi = 300)


# Differential expression annalysis
OR_obj_marker<-lapply(1:3,function(i){
  temp_marker<-FindMarkers(OR_obj[[i]],
                        group.by = "OR_allele_exptype",
                        ident.1 ="Paternal_specific",
                        ident.2="Maternal_specific") 
  temp_marker<-temp_marker%>%mutate(Group=case_when(avg_log2FC>=1.5&p_val<0.05 ~ "Up",abs(avg_log2FC)<1.5|p_val>=0.05 ~ "NotSignificant",avg_log2FC<=(-1.5)&p_val<0.05 ~ "Down"))
  temp_marker$gene<-rownames(temp_marker)
  temp_marker
})

# save the differential expression gene
wb <- createWorkbook()
# Add each dataframe in OR_obj_marker sequentially as different sheets
for (i in 1:length(OR_obj_marker)) {
  addWorksheet(wb, myOR_far[i])  
  df <- OR_obj_marker[[i]]
  df <- df[, c(ncol(df), 1:(ncol(df)-1))]
  writeData(wb, sheet = myOR_far[i], df)
}
saveWorkbook(wb, "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_differential_gene.xlsx", overwrite = TRUE)


# Figure 5d and Extended Data Fig.6a
# marker gene volcano plot
library(ggrepel)
marker_Volcano_Plot<-function(df){
  df$gene<-rownames(df)
  temp_p<-ggplot(df, aes(x = avg_log2FC, y = -log10(p_val),colour=Group)) +geom_point(alpha=0.8, size=2) +
    geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.5) +scale_color_manual(values=c( "#8fb4be","grey","#D93f49"))+
    geom_vline(xintercept = c(-1.5,1.5),lty=4,col="black",lwd=0.5)+
    scale_x_continuous(limit=c(-5,5))+
    # 坐标轴
    labs(x="avg_log2FC (Paternal/Maternal)",
      y="-log10 (p-value)",title=str_c("Differentially expressed genes between Paternal and Maternal specific cells"))+
    theme_classic()+
    # 图例
    theme(plot.title = element_text(hjust = 0.5), 
      legend.position="right", 
      legend.title = element_blank(),
      text=element_text(size=6)
  )
  temp_p<-temp_p+geom_text_repel(data = df%>%subset(Group%in%c("Up","Down")&(-log10(p_val))>=2), aes(x = avg_log2FC, y = -log10(p_val), label = gene),
                  size = 2,box.padding = unit(0.1, "lines"),
                  point.padding = unit(0.3, "lines"), 
                  segment.color = "black", max.overlaps=30,
                  show.legend = FALSE)
  return(temp_p)
}
rm(marker_p)
for(i in 1:4){
  temp_p<-marker_Volcano_Plot(OR_obj_marker[[i]])
  if(exists("marker_p")){
    marker_p<-marker_p+temp_p
  }else{
    marker_p<-temp_p
  }
}
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_marker_Volcano_Plot.pdf",marker_p+plot_layout(ncol=4,guides="collect"), device = "pdf", width = 34, height = 8, units = "cm", dpi = 300,bg="white")

# marker gene heat map
library(ComplexHeatmap)
library(circlize)
col_fun <- colorRamp2(
  c(-2,0,1, 4), 
  c("#8fb4be","#d5e1e3","#EBBFC2","#D93f49")
)
get_dist<-function(x,y){
  return(1-cosine(x,y))
}
OR_DEG_heat<-function(i){
  temp_dscore<-OR_dscore%>%subset(gene==myOR_far[i]&merge_dscore_type!="Biallelic")
  marker<-OR_obj_marker[[i]]%>%subset(Group%in%c("Up","Down"))
  ha2 <- HeatmapAnnotation(
      df = temp_dscore[,15,drop=FALSE],
      col = list(merge_dscore_type=c("Paternal_specific"="#94b2c2", "Biallelic"="#dfe2d6", "Maternal_specific"="#c68989"))
  )
  h<- Heatmap(filtered_mature_OSNs[["RNA"]]@scale.data[rownames(marker),temp_dscore$barcode],clustering_distance_columns=get_dist,clustering_distance_rows=get_dist,top_annotation=ha2,col =col_fun,name="RNA z-score",show_row_dend=FALSE,show_row_names=FALSE,show_column_names=FALSE,border_gp = gpar(col = "black"))
  return(draw(h,heatmap_legend_side = "right",annotation_legend_side = "right"))
}
pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_DEG_heat.pdf",height=4.5,width=8)
for(i in 1:3){
  OR_DEG_heat(i)
}
dev.off()


# differential peak annalysis
OR_obj_Peakmarker<-lapply(1:3,function(i){
  DefaultAssay(OR_obj[[i]])<-"peak"
  temp_marker<-FindMarkers(OR_obj[[i]],
                        group.by = "OR_allele_exptype",
                        ident.1 ="Paternal_specific",
                        ident.2="Maternal_specific",
                        test.use = 'LR') 
  temp_marker<-temp_marker%>%mutate(Group=case_when(avg_log2FC>=1.5&p_val<0.05 ~ "Up",abs(avg_log2FC)<1.5|p_val>=0.05 ~ "NotSignificant",avg_log2FC<=(-1.5)&p_val<0.05 ~ "Down"))
})
wb <- createWorkbook()
# Add each dataframe in OR_obj_Peakmarker sequentially as different sheets
for (i in 1:length(OR_obj_Peakmarker)) {
  addWorksheet(wb, myOR_far[i])  
  df <- OR_obj_Peakmarker[[i]]
  df$peak<-rownames(OR_obj_Peakmarker[[i]])
  df <- df[, c(ncol(df), 1:(ncol(df)-1))]
  writeData(wb, sheet = myOR_far[i], df)
}
saveWorkbook(wb, "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_differential_peak.xlsx", overwrite = TRUE)


# differential peak heat map
DefaultAssay(filtered_mature_OSNs) <- "peak"
filtered_mature_OSNs <- ScaleData(filtered_mature_OSNs)
col_fun_peak <- colorRamp2(
  c(-1,0,2, 4), 
  c("#5c8b9b","#ccd5dc","#e0c7c3","#8d4550")
)
OR_DEP_heat<-function(i){
  temp_dscore<-OR_dscore%>%subset(gene==myOR_far[i])
  marker<-OR_obj_Peakmarker[[i]]%>%subset(Group%in%c("Up","Down"))
  ha2 <- HeatmapAnnotation(
      df = temp_dscore[,15,drop=FALSE],
      col = list(merge_dscore_type=c("Paternal_specific"="#94b2c2", "Biallelic"="#dfe2d6", "Maternal_specific"="#c68989"))
  )
  h<- Heatmap(filtered_mature_OSNs[["peak"]]@scale.data[rownames(marker),temp_dscore$barcode],clustering_distance_columns=get_dist,clustering_distance_rows=get_dist,top_annotation=ha2,col =col_fun_peak,name="ATAC z-score",show_row_dend=FALSE,show_row_names=FALSE,show_column_names=FALSE,border_gp = gpar(col = "black"))
  return(draw(h,heatmap_legend_side = "right",annotation_legend_side = "right"))
}

pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/OR_DEP_heat.pdf",height=4.5,width=6.5)
for(i in 1:3){
  OR_DEP_heat(i)
}
dev.off()



# peak annotation target gene is the nearest distance of tss
OR_obj_significantPeak_anno<-lapply(1:3,function(i){
  Significant<-OR_obj_Peakmarker[[i]]%>%subset(Group!="NotSignificant")
  peak_dist<-lapply(1:nrow(Significant),function(j){
    peak<-rownames(Significant)[j]
    peak<-unlist(strsplit(peak,":"))
    peak<-unlist(strsplit(peak,"-"))
    data<-tss_bed%>%subset(chr==peak[1])
    data$peak_start<-peak[2]
    data$peak_end<-peak[3]
    data<-data%>%mutate(dist=case_when(strand=="+"~start-as.numeric(peak_start),strand=="-"~as.numeric(peak_start)-start))%>%subset(abs(dist)<=10**5)%>%rename_with(~ c("geneid"), 4)
    data<-merge(data,id2gene,by="geneid")
    data$peak_Group<-rep(Significant$Group[j],nrow(data))
    data
  })
  peak_dist<-do.call(rbind,peak_dist)
})

OR_differential_expression_peak_combined<-lapply(1:3,function(i){
  temp<-merge(OR_obj_significantPeak_anno[[i]],OR_obj_marker[[i]]%>%subset(Group!="NotSignificant"),by="gene")
  temp_upstream<-temp%>%subset(dist>=0)
  temp_downstream<-temp%>%subset(dist<0)
  temp_upstream_min<-temp_upstream%>%group_by(gene)%>%summarise(dist=min(dist))
  temp_downstream_min<-temp_downstream%>%group_by(gene)%>%summarise(dist=max(dist))
  temp_upstream<-merge(temp_upstream,temp_upstream_min,by=c("gene","dist"))
  temp_downstream<-merge(temp_downstream,temp_downstream_min,by=c("gene","dist"))
  a<-rbind(temp_downstream,temp_upstream)
  a<-a%>%subset(peak_Group==Group)
})

# filtered the highly confidence differential expression-peak
peak_gene<-list()
peak_gene[[1]]<-c("Amfr","Kirrel3","Pak7","Rps24","Bsn","Hnrnpa1","Lsamp")
peak_gene[[2]]<-c("Pcdh7","Pcdh9","Calb2")
peak_gene[[3]]<-c("Kirrel","Kirrel3","Lrrc36","Cd55","Lrrc36")
OR_dep_highly<-list()
for(i in 1:3){ 
  OR_dep_highly[[i]]<-OR_differential_expression_peak_combined[[i]]%>%subset(gene%in%peak_gene[[i]])
  OR_dep_highly[[i]]$OR<-myOR_far[i]
}
OR_dep_highly<-do.call(rbind,OR_dep_highly)

# expression and peak track
# Figure 5b & Extended Data Fig.6b-d
DefaultAssay(filtered_mature_OSNs)<-"ATAC"
OR_dep_track<-function(i){
  peak<-OR_dep_highly[i,]
  dist<-peak[1,2]
  if(dist<0&peak[1,8]=="-"){
    total_bed<-str_c(peak[1,4],"-",as.numeric(peak[1,9])-5000,"-",as.numeric(peak[1,5]))
  }else if(dist<0&peak[1,8]=="+"){
    total_bed<-str_c(peak[1,4],"-",as.numeric(peak[1,5]),"-",as.numeric(peak[1,10])+5000)
  }else if(dist>0&peak[1,8]=="-"){
    total_bed<-str_c(peak[1,4],"-",as.numeric(peak[1,5])-10000,"-",as.numeric(peak[1,10])+1000)
  }else{
    total_bed<-str_c(peak[1,4],"-",as.numeric(peak[1,9])-1000,"-",as.numeric(peak[1,10])+10000)
  }
  ######################################
  OR_index<-which(myOR_far==peak[1,18])
  obj<-OR_obj[[OR_index]]
  OR_cell<-OR_dscore%>%subset(gene==peak[1,18]&merge_dscore_type!="Biallelic")
  obj<-subset(obj,cells=OR_cell$barcode)
  ##########补充annotation plot################
  DefaultAssay(obj)<-"ATAC"
  ####画atac allele图####
  ranges.show <- StringToGRanges(str_c(peak[1,4],"-",peak[1,9],"-",peak[1,10]))
  ranges.show$color <- "#584c7c"
  ap<-AnnotationPlot(object = obj, region = total_bed)+
      theme(text=element_text(size=8))
  p<-CoveragePlot(
        object = obj,
        region = total_bed,
        group.by="OR_allele_exptype",
        tile =FALSE,
        peaks=FALSE,   
        annotation=FALSE,     
        region.highlight = ranges.show)&scale_fill_manual(values=c("#b45a55","#505f84"))&
        theme(axis.title.x = element_blank(),  # 隐藏 y 轴标题
        axis.text.x=element_blank(),   # 隐藏 y 轴刻度标签
        axis.ticks.x= element_blank(),   # 隐藏 y 轴刻度线
        axis.line.x= element_blank(),
        text=element_text(size=8)
        )
  total_p<-CoveragePlot(
        object = filtered_mature_OSNs,
        region = total_bed,
        tile =FALSE,
        peaks=FALSE,  
        annotation=FALSE,     
        region.highlight = ranges.show)&scale_fill_manual(values=c("#c5cd94"))&
        theme(axis.title.x = element_blank(),  # 隐藏 y 轴标题
        axis.text.x=element_blank(),   # 隐藏 y 轴刻度标签
        axis.ticks.x= element_blank(),   # 隐藏 y 轴刻度线
        axis.line.x= element_blank(),
        text=element_text(size=8)
        )
  exp<-ExpressionPlot(obj, features = peak[1,1],group.by="OR_allele_exptype", assay = "RNA")+scale_fill_manual(values=c("#b45a55","#505f84"))+
      theme(text=element_text(size=8))
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/OR_gene_dep/",myOR_far[OR_index],"_",peak[1,1],"_",peak[1,2],".pdf")
  ggsave(filename = file,total_p+plot_spacer()+p+exp+ap+plot_spacer()+plot_layout(ncol=2,widths=c(1.5,0.5),heights=c(1,2,1)),device = "pdf", width =15, height =8, units = "cm", dpi = 300,bg="white")
  return()
}

for(i in 1:nrow(OR_dep_highly)){
  OR_dep_track(i)
}


# Plot using non-synonymous mutations
# Figure 5c
# The number of SNPs in each OR CDS vs length of CDS
SNP_VS_length_p<-ggplot(data=final_OR_CCDS,aes(x=CCDS_length,y=unsyn_snp_counts))+
  geom_point(color="#799399")+
  labs(x="Length of CCDS per OR gene",y="# Non-synonymous variants per OR CCDS")+
  theme_classic()+
  theme(
        axis.text = element_text(colour = 'black',size = 6),
        panel.grid = element_blank(),
        text = element_text(size = 8))
# Parental distance VS non-synonymous mutations count
distance_VS_snp_p<-ggplot(data=OR_distance, aes(x =unsyn_snp_counts , y =between_median-within_median)) +
  geom_point(aes(color=wilcox.test,size=cellnumber),alpha=0.8) +scale_y_break(breaks=c(0.3,0.5),ticklabels=seq(0.5,1,0.1),scales=0.3)+
  geom_text(data=OR_distance%>%subset(wilcox.test=="p.value<0.05"), aes(x =unsyn_snp_counts+0.2 , y =between_median-within_median-0.04,label = gene), vjust =-0.5, size = 2.5)+
  scale_color_manual(values=c("#587880","grey"))+geom_smooth(method="lm",span=1, size=0.5,  se = TRUE)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., after_stat(p.value.label), sep = '~`,`~')),
             formula = y~x, parse = TRUE, size = 2,color="red")+
  labs(y="transcriptome difference (between-within)",
       x="# Non-synonymous variants in CCDS for each ORN group",legend)+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=8)
) 
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/figure/Figure6/unsyn_snp_distance.pdf",SNP_VS_length_p+distance_VS_snp_p+plot_layout(width=c(1,1.5)), device = "pdf", width = 21, height = 8, units = "cm", dpi = 60,bg="white")


# Figure 5d
# get genes with unsyn_snp == 1 and high and low difference
OR_distance_snp1<-OR_distance%>%subset(unsyn_snp_counts==1&(difference>0.1|difference<0.02))
OR_final_ccds<-read.table(file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/OR_final_ccds.txt",header=TRUE)
OR_GPCR_ccds<-OR_final_ccds%>%subset(gene%in%OR_distance_snp1$gene)%>%group_by(gene)%>%dplyr::slice(1)%>%select(gene,geneId)

# select class II OR
class_OR<-read.xlsx("/data/R02/huangtt39/data/OSN/mmc3.xlsx",sheet = 2)
class_OR<-class_OR%>%mutate(class=case_when(is.Class.I.OR=="yes"~"class_I",is.Class.I.OR=="no"~"class_II"))
class_OR<-class_OR%>%subset(class %in% c("class_I","class_II"))%>%select(gene.name,class)

# get the txid of gene to mark the non-synonymous mutations
OR_GPCR_ccds<-merge(class_OR,OR_GPCR_ccds,by.x="gene.name",by.y="gene")%>%subset(class=="class_II")
write.table(OR_GPCR_ccds%>%select(-2),file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/GPCR/OR_GPCR_txid.txt",col.names = FALSE, row.names = FALSE,quote = FALSE,sep = "\t")

# Output all gene IDs of class II ORs for constructing the class II GPCR consensus sequence
write.table(merge(class_OR,avail_OR_CCDS,by.x="gene.name",by.y="gene")%>%subset(class=="class_II")%>%group_by(gene.name)%>%dplyr::slice(1)%>%select(gene.name,geneId)
,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/OR_classII_CCDS_txid.txt",col.names = FALSE, row.names = FALSE,quote = FALSE,sep = "\t")

# Import the constructed class II OR consensus sequence
classII_OR_consense<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/conservation.xlsx",sheet = 1)

# 使用ggplot绘制条形图
GPRC_p<-ggplot(classII_OR_consense, aes(x = 1:302, y = round(Max.ratio*100,2), fill =round(Max.ratio*100,2) )) +
  geom_bar(stat = "identity") +  # 绘制条形图
  scale_fill_gradient2(
    low = "#f6eae0",   # 浅蓝
    mid = "#dfdde4",   # 紫色（中间色）
    high = "#414887",  # 深红
    midpoint = 40     # 设置中间值的位置为0.5，适合数据范围从0到1
  ) +                            # 使用自定义颜色
  labs(title = "Barplot with GPCR color", 
       x = "Index", 
       y = "Value") +
  theme_classic() +                                  # 使用简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # x轴标签倾斜
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/OR_annalysis/GPCR/GPCR_color.pdf",GPRC_p, device = "pdf", width = 30, height = 5, units = "cm", dpi = 300)

gg_build <- ggplot_build(GPRC_p)

# 提取每个条形的颜色
bar_colors <- gg_build$data[[1]]$fill

classII_OR_consense$color<-bar_colors

write.xlsx(classII_OR_consense,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/GPCR/conservation2color.xlsx",sheet = 1)

distance_conversation<-classII_OR_consense%>%subset(distance%in%c("high","low"))
distance_conversation$distance<-factor(distance_conversation$distance,levels=c("low","high"))
#根据distance 高低画散点图
distance_conversation_p<-ggplot(distance_conversation, aes(x = distance, y = Max.ratio)) +
              geom_boxplot() +  
              geom_jitter(aes(color = distance), width = 0.2, size = 3, alpha = 0.6) +  # 添加散点
              scale_color_manual(values=c("#0E8A8F","#F44040"))+
              theme_classic()+ 
              labs(y="Animo acid conservation",x="Transcriptome difference between two alleles")+
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=6))

distance_conversation_p<-ggplot(classII_OR_consense%>%subset(distance%in%c("high","low")), aes(x = difference, y = Max.ratio, color = distance)) +
  geom_point(size = 4) +  # 使用geom_point绘制散点
  theme_classic() +  # 使用简洁主题
  scale_color_manual(values=c("#0E8A8F","#F44040"))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure6/GPCR/distance_conversation.pdf",distance_conversation_p, device = "pdf", width = 10, height = 5, units = "cm", dpi = 300)
