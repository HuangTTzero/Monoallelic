library(ggplot2)
library(dplyr)
library(dbplyr)
library(magrittr)
library(data.table)
library(openxlsx)
library(sctransform)
library(SummarizedExperiment)
library(stringr)
library(hrbrthemes)
library(readr)
library("parallel")
library(tidyr) 
library("ggunchained")
library(patchwork)
###showtext###
library("sysfonts")
library("showtextdb")
library("showtext")
library("bedtoolsr")
library(ggpmisc)
library(ggpubr)
library(purrr)
library(ggExtra)
showtext_auto()

# 突变频率数据
genome <- 2725521370
genome_exon <- 110722912

Cast_snp <- 20668274
Cast_snp_exon <- 0.01638 * Cast_snp
Cast_genome_ratio <- Cast_snp / genome 
Cast_exon_ratio <- Cast_snp_exon / genome_exon 

DBA_snp <- 5184367
DBA_snp_exon <- 0.01748 * DBA_snp
DBA_genome_ratio <- DBA_snp / genome
DBA_exon_ratio <- DBA_snp_exon / genome_exon 

# 创建数据框 - 按strain分组
rate_data <- data.frame(
  strain = c("Cast", "Cast", "DBA", "DBA"),
  region = c("whole genome", "exon", "whole genome", "exon"),
  frequency = c(Cast_genome_ratio, Cast_exon_ratio, DBA_genome_ratio, DBA_exon_ratio)
)

# 确保因子顺序
rate_data$strain <- factor(rate_data$strain, levels = c("Cast", "DBA"))
rate_data$region <- factor(rate_data$region, levels = c("whole genome", "exon"))

# 绘图 - 按strain分组，x轴为strain
rate_p <- ggplot(data = rate_data, 
                 aes(x = strain, y = frequency, fill = region)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%0.4f", round(frequency, digits = 4))),
            position = position_dodge(width = 0.7),
            size = 3, vjust = -0.5) +
  scale_fill_manual(values = c("whole genome" = "#2f4553", 
                               "exon" = "#e8cacb")) +
  labs(x = "Strain", y = "SNV density", 
       fill = "Region") +
  theme_classic() +
  theme(legend.key.size = unit(8, "pt"),
        text = element_text(size = 8),
        legend.position = "top")


# 读取LiMCA数据
LiMCA_reads <- read.table("~/F1_OSN/mapping/LiMCA/reads_distribution.txt", 
                          header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 根据Sample名称添加strain分类
LiMCA_reads$strain <- ifelse(grepl("bc", LiMCA_reads$Sample), "Cast", 
                             ifelse(grepl("bdf1", LiMCA_reads$Sample), "DBA", "Other"))

# 将百分比列转换为数值
LiMCA_reads$G1_Percentage_numeric <- as.numeric(sub("%", "", LiMCA_reads$G1_Percentage)) / 100
LiMCA_reads$G2_Percentage_numeric <- as.numeric(sub("%", "", LiMCA_reads$G2_Percentage)) / 100

# 转换数据为长格式
library(tidyr)
LiMCA_reads_long <- pivot_longer(
  LiMCA_reads,
  cols = c(G1_Percentage_numeric, G2_Percentage_numeric),
  names_to = "Group",
  values_to = "Percentage"
)

# 简化Group名称
LiMCA_reads_long$Group <- ifelse(LiMCA_reads_long$Group == "G1_Percentage_numeric", "G1", "G2")

# 确保因子顺序
LiMCA_reads_long$strain <- factor(LiMCA_reads_long$strain, levels = c("Cast", "DBA"))
LiMCA_reads_long$Group <- factor(LiMCA_reads_long$Group, levels = c("G1", "G2"))

# 绘图 - 按strain分组，x轴为strain
reads_vln_plot <- ggplot(LiMCA_reads_long, aes(x = strain, y = Percentage, fill = Group)) +
  geom_violin(position = position_dodge(width = 0.8), 
              alpha = 0.9, 
              width = 0.7,
              trim = FALSE) +
  geom_boxplot(width = 0.15, 
               position = position_dodge(width = 0.8),
               outlier.shape = NA,
               alpha = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.2, 
                                             dodge.width = 0.8),
             size = 1, 
             alpha = 0.6) +
  scale_fill_manual(values = c("G1" = "#526188", "G2" = "#bb5c56")) +
  labs(x = "Strain", y = "Percentage", fill = "Group") +
  theme_classic() +
  theme(
    text = element_text(size = 10),
    legend.position = "top",
    axis.text = element_text(color = "black")
  ) +
  scale_y_continuous(labels = scales::percent_format())

ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/LiMCA_reads.pdf",rate_p+reads_vln_plot, device = "pdf", width = 30, height =10 , units = "cm", dpi = 60,bg="white")


rfdir<-"~/F1_OSN/mapping/reference/mm10-mask/"
id2gene<-read.table(file=str_c(rfdir,"genes/id_to_name.txt"))
id2gene<-id2gene[!duplicated(id2gene$V2),]
colnames(id2gene)<-c("geneid","gene")
trans_bed<-read.table(file=str_c(rfdir,"regions/transcripts.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
tss_bed<-read.table(file=str_c(rfdir,"regions/tss.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
promoter_bed<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-1000-1,
                         strand=="-" ~end-1000),end=case_when(strand=="+" ~ end+1000-1,
                         strand=="-" ~end+1000))
promoter_bed<-merge(promoter_bed,id2gene,by.x="gene",by.y="geneid")
chrom<-paste0("chr",c(1:19,"X","Y"))
tss_bed_new<-tss_bed%>%subset(chr %in% chrom)
chr_gene<-id2gene%>%subset(geneid  %in% tss_bed_new$gene)
chrXgene<-unique(trans_bed$gene[! is.na(str_extract(trans_bed$chr,"chrX"))])
chrXgene<-id2gene%>%subset(geneid %in% chrXgene)
chrXgene<-chrXgene$gene
# load OR names
ORgene<-read.xlsx("~/F1_OSN/CODE/data/OSN/12864_2020_6583_MOESM2_ESM.xlsx",sheet = 2)
ORgene <- unique(ORgene$Gene.symbol)

G1_counts <- read.table(
  "~/F1_OSN/mapping/LiMCA/featurecounts/total_G1.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
G2_counts <- read.table(
  "~/F1_OSN/mapping/LiMCA/featurecounts/total_G2.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
ori_counts <- read.table(
  "~/F1_OSN/mapping/LiMCA/featurecounts/total_dedup.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

# 3. 转换函数 
convert_geneid_to_gene <- function(counts_data, gene_mapping) {
  counts_data$GeneID <- rownames(counts_data)
  
  merged <- merge(counts_data, gene_mapping, 
                  by.x = "GeneID", by.y = "geneid")
  
  rownames(merged) <- merged$gene
  
  merged <- merged[, !colnames(merged) %in% c("GeneID", "gene")]

  count_matrix <- as.matrix(merged)
  
  return(count_matrix)
}

G1_matrix <- convert_geneid_to_gene(G1_counts, id2gene)
G2_matrix <- convert_geneid_to_gene(G2_counts, id2gene)
Ori_matrix <- convert_geneid_to_gene(ori_counts, id2gene)

LiMCA_readList<-list()
LiMCA_readList[["G1_rna"]]<-G1_matrix
LiMCA_readList[["G2_rna"]]<-G2_matrix
LiMCA_readList[["rna"]]<-LiMCA_readList[["G1_rna"]]+LiMCA_readList[["G2_rna"]]
LiMCA_readList[["Orirna"]]<-Ori_matrix

LiMCA_merge_count <- cbind(
  # 第一列：rna大于0的计数
  as.data.frame(rowSums(LiMCA_readList[["rna"]] > 0)),
  
  # 第二列：rna大于0的均值（需要处理NA）
  as.data.frame({
    rna_means <- rowMeans(
      replace(LiMCA_readList[["rna"]], LiMCA_readList[["rna"]] <= 0, NA), 
      na.rm = TRUE
    )
    # 将NaN转换为NA
    rna_means[is.nan(rna_means)] <- NA
    rna_means
  }),
  
  # 第三列：Orirna大于0的计数
  as.data.frame(rowSums(LiMCA_readList[["Orirna"]] > 0))
) %>% 
  # 重命名列 - 注意顺序：第一列是rna计数，第二列是rna均值，第三列是Orirna计数
  rename_with(~ c("rna_count", "rna_mean", "Orirna_count"), 1:3) %>% 
  # 过滤：只保留Orirna_count > 0的行
  filter(Orirna_count > 0)
#27575
# add gene type
LiMCA_merge_count$genetype<-rep("Others",nrow(LiMCA_merge_count))
for (i in 1:nrow(LiMCA_merge_count)) { 
  gene <- rownames(LiMCA_merge_count)[i]
  if(gene %in% chrXgene ){
    type="chrX"
    LiMCA_merge_count[i,]$genetype=type
  }
  if(gene %in% ORgene ){
    type="OR"
    LiMCA_merge_count[i,]$genetype=type
  }
} 
LiMCA_merge_count$gene<-rownames(LiMCA_merge_count)
# get gene in chromsomes
LiMCA_merge_count<-LiMCA_merge_count%>%subset(gene %in%chr_gene$gene)
LiMCA_merge_count%>%group_by(genetype)%>%summarize(count=n())
# # A tibble: 3 × 2
#   genetype count
#   <chr>    <int>
# 1 OR        1048
# 2 Others   25696
# 3 chrX       787
save(LiMCA_merge_count,file="~/F1_OSN/analysis/merge/RData/LiMCA_merge_count.RData")

# get gene_cell_reads_dataframe
cl <- makeCluster(6)
clusterExport(cl, c("LiMCA_merge_count","LiMCA_readList"))
LiMCA_dscore_data <- parLapply(cl,1:nrow(LiMCA_merge_count),function(i){
  library(Seurat)
  library(dplyr)
  gene <- rownames(LiMCA_merge_count)[i]
  cells<-which(LiMCA_readList[["Orirna"]][gene,]>0)
  temp_origene<-LiMCA_readList[["Orirna"]][gene,]
  final_oriexp<-temp_origene[cells]
  temp_G1gene<-LiMCA_readList[["G1_rna"]][gene,]
  temp_G2gene<-LiMCA_readList[["G2_rna"]][gene,]
  final_G1gene<-temp_G1gene[cells]
  final_G2gene<-temp_G2gene[cells]
  temp_matrix<-cbind(as.data.frame(rep(gene, length(final_oriexp))),as.data.frame(rep(LiMCA_merge_count[i,]$genetype, length(final_oriexp))),as.data.frame(final_G1gene),as.data.frame(final_G2gene),as.data.frame(final_oriexp)) %>% rename_with(~ c("gene","type","G1rna","G2rna","Orirna"), 1:5)
  temp_matrix
})
stopCluster(cl)
LiMCA_dscore_data<-do.call(rbind,LiMCA_dscore_data)
LiMCA_dscore_data$barcode <- str_sub(rownames(LiMCA_dscore_data), start = 1,end=20)

LiMCA_dscore_data$rnascore<-LiMCA_dscore_data$G1rna/(LiMCA_dscore_data$G1rna+LiMCA_dscore_data$G2rna)-0.5
LiMCA_dscore_data$allele_exp<-LiMCA_dscore_data$G1rna+LiMCA_dscore_data$G2rna
save(LiMCA_dscore_data,file="~/F1_OSN/analysis/merge/RData/LiMCA_dscore_data.RData")
load(file="~/F1_OSN/analysis/merge/RData/LiMCA_dscore_data.RData")

# gene number with allelic reads
# Extended Data Fig. 2d
LiMCA_Orirna<-LiMCA_merge_count%>%subset(Orirna_count>0)%>%nrow()
LiMCA_allelerna<-LiMCA_merge_count%>%subset(rna_count>0)%>%nrow()
LiMCA_merge_count2plot<-data.frame(allele_type=c("total","allele"),count=c(LiMCA_Orirna,LiMCA_allelerna))
LiMCA_merge_count2plot$allele_type<-factor(LiMCA_merge_count2plot$allele_type,levels=c("total","allele"))

LiMCA_merge_count2plot_p<-ggplot(LiMCA_merge_count2plot, aes(x = allele_type, y = count, fill = allele_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("#dee2d6","#5c787c")) +  # 设置颜色
  scale_y_continuous(limits=c(0,30000))+
  labs(x = "", y = "# genes") +
  theme_classic() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/LiMCA_gene_count.pdf", LiMCA_merge_count2plot_p, device = "pdf", width = 8, height = 7, units = "cm", dpi = 300)


LiMCA_allele_count_meanRNA_p<-ggplot() +
  geom_point(data=LiMCA_merge_count%>%subset(rna_count>0), aes(x = rna_count, y =log10(rna_mean)),color="#567074",size=1,alpha=0.9) +
  labs(x="Number of cells with allelic reads",
       y="Mean allelic reads of cells (log10)")+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="None",text=element_text(size=10)
)  
LiMCA_allele_count_p <-ggplot(data=LiMCA_merge_count%>%subset(rna_count>0),aes(x=rna_count))+
      geom_histogram(bins = 20,
                 color="white",
                 fill="grey",alpha=0.5)+
      scale_y_log10() +
      theme_classic()+theme(text=element_text(size=6))+
      labs(y="gene number",
          x="allelic informatic cell number")
LiMCA_allele_meanRNA_p <-ggplot(data=LiMCA_merge_count%>%subset(rna_count>0),aes(x=log10(rna_mean)))+
      geom_histogram(bins = 20,
                 color="white",
                 fill="grey",alpha=0.5)+
      scale_y_log10() +
      theme_classic()+theme(text=element_text(size=6))+
      labs(y="gene number",
          x="Mean allelic reads of cells (log10)")+
      coord_flip()
ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/LiMCA_allele_count_meanRNA.pdf",LiMCA_allele_count_p+plot_spacer()+LiMCA_allele_count_meanRNA_p+LiMCA_allele_meanRNA_p+plot_layout(nrow=2,ncol=2,height=c(1,2,1,2),width=c(2,1,2,1)), device = "pdf", width = 10, height = 9, units = "cm", dpi = 300,bg="white")

load(file="~/F1_OSN/analysis/merge/RData/final_gene.RData")

color<-c("#9b9b9a","#d19bd1","#1c2d58")
cutoff<-function(LiMCA_dscore_data){
  final<-data.frame()
  ori<-LiMCA_dscore_data%>%group_by(gene,type)%>%summarise(count=sum(Orirna>0))%>%mutate(cutofftype="total UMIs (>=1)")
  cutoff_1<-LiMCA_dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>=10))%>%mutate(cutofftype="allele UMIs (>=10)")
  cutoff_2<-LiMCA_dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>=50))%>%mutate(cutofftype="allele UMIs (>=50)")
  final<-rbind(final,ori,cutoff_1,cutoff_2)
  final$cutofftype<-factor(final$cutofftype,levels=c("total UMIs (>=1)","allele UMIs (>=10)","allele UMIs (>=50)"))
  return(final)
} 
LiMCA_readscutoff<-cutoff(LiMCA_dscore_data%>%subset(gene%in%final_gene$gene))
LiMCA_readscutoff$type<-factor(LiMCA_readscutoff$type,levels=c("Others","chrX","OR"))
cutoff_p<-function(data){
  p<-ggplot(data, aes(x = log10(count+1), fill = cutofftype),alpha=0.6) +geom_histogram(alpha=0.6,position="identity",bins = 20) +
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text=element_text(size=6), axis.text = element_text(colour = 'black',size = 6),)+
  scale_fill_manual(values=color,name="Cell cutoff")+ylab("gene number")+xlab("log10(cell numbers+1)")
  return(p)
}
ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/LiMCA_cutoff_reads_total.pdf",cutoff_p(LiMCA_readscutoff)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =12, height = 6, units = "cm", dpi = 300,bg="white")


# 为了保存OR基因，这里先不筛选cell number >=10 的基因
cellcolor<-c("#597b83","#d8bf87","#d9c3c5")
LiMCA_final_gene<-LiMCA_readscutoff%>%subset(cutofftype=="allele UMIs (>=10)"&count>=10)
LiMCA_gene_dscore_10<-LiMCA_dscore_data%>%subset(allele_exp>=10) 
# Give the high reads of smart-seq2,here a gene was called monoallelically expressed in a cell if one allele contributed ≥98% of the genotype-informative bases.
LiMCA_gene_dscore_10<-LiMCA_gene_dscore_10%>%mutate(rna_dscore_type=case_when(rnascore<=-0.48 ~ "Maternal_specific",rnascore>(-0.48)&rnascore<0.48 ~ "Biallelic",rnascore>=0.48 ~ "Paternal_specific"))

# the expression pattern of each gene 
# add positive control, imprinting gene

LiMCA_gene_10_celltype_count<-LiMCA_gene_dscore_10%>%group_by(gene,rna_dscore_type)%>%summarise(cellnumber=n())
LiMCA_gene_10_celltype_count<-merge(LiMCA_gene_10_celltype_count,LiMCA_readscutoff%>%subset(cutofftype=="allele UMIs (>=10)")%>%dplyr::select(gene,count,type),by="gene")
LiMCA_10_count_ratio<-LiMCA_gene_10_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="rna_dscore_type",values_from="ratio",id_cols=c(type,gene,count),values_fill=0)
LiMCA_10_count<-LiMCA_gene_10_celltype_count%>%pivot_wider(names_from="rna_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count),values_fill=0)

LiMCA_10_count_ratio$type<-as.character(LiMCA_10_count_ratio$type)
LiMCA_10_count_ratio<- LiMCA_10_count_ratio %>%mutate(type = ifelse(gene %in% c("Peg3","Kcnq1ot1"), "imprinting", type))
LiMCA_10_count_ratio$type<-factor(LiMCA_10_count_ratio$type,levels=c("OR","chrX","Others","imprinting"))

LiMCA_10_count_ratio_OR<-LiMCA_10_count_ratio%>%subset(type=="OR")

LiMCA_final_count_ratio<-LiMCA_10_count_ratio%>%subset(gene%in% LiMCA_final_gene$gene)
save(LiMCA_final_count_ratio,file="~/F1_OSN/analysis/merge/RData/LiMCA_final_count_ratio.RData")

LiMCA_to_Positive<-rbind(LiMCA_10_count_ratio_OR,LiMCA_final_count_ratio%>%subset(type!="OR"))

# 首先统计每个位置点的数量
LiMCA_to_Positive_counts <- LiMCA_to_Positive %>%
  subset(type %in% c("OR", "chrX", "imprinting")) %>%
  group_by(Maternal_specific, Paternal_specific, type) %>%
  summarise(count = n(), .groups = "drop")

LiMCA_Positive_p <- ggplot() +
  # 使用统计后的数据，点的大小映射到count
  geom_point(data = LiMCA_to_Positive_counts, 
             aes(x = Maternal_specific, 
                 y = Paternal_specific, 
                 color = type, 
                 size = count), 
             alpha = 0.9) +
  
  # 设置颜色
  scale_color_manual(values = c(cellcolor[1:2], "#2323AF")) +
  
  # 设置点的大小范围（可根据需要调整）
  scale_size_continuous(range = c(1, 5), 
                        breaks = c(1, 5, 10, 20), 
                        name = "Number of points") +
  
  # 设置坐标轴范围
  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +
  
  # 坐标轴标签
  labs(x = "% Maternal specific cell",
       y = "% Paternal specific cell") +
  
  theme_classic() +
  
  # 图例设置
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right",
        text = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(colour = 'black', size = 8),
        legend.key.size = unit(10, "pt"))

# 保存图形
ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/Positve.pdf",
       LiMCA_Positive_p, 
       device = "pdf", 
       width = 11,  # 稍微增加宽度以容纳大小图例
       height = 7, 
       units = "cm", 
       dpi = 300,
       bg = "white")
# Add OR genes


high_RME<-read.xlsx("~/F1_OSN/new_analysis/high_RME.xlsx")
high_RME<-high_RME%>%arrange(desc(1-Biallelic), desc(Maternal_specific)) 

LiMCA_high_RME_ratio<-LiMCA_10_count_ratio%>%subset(gene%in%high_RME$gene) # 23 genes ,the Bamp6 has no allelic reads
LiMCA_high_RME_ratio$gene<-factor(LiMCA_high_RME_ratio$gene,levels=high_RME$gene)
LiMCA_high_RME_count<-LiMCA_10_count%>%subset(gene%in%high_RME$gene)

LiMCA_high_RME_count_p<- ggplot(LiMCA_high_RME_ratio, aes(x = gene, y = count)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_classic() +
  labs(x = "Gene", y = "Total cell") +
  theme(axis.text.x = element_blank())    

LiMCA_high_RME_ratio2p<-LiMCA_high_RME_ratio%>% 
            pivot_longer(cols = c(Paternal_specific, Maternal_specific, Biallelic), 
               names_to = "dscore_type", values_to = "ratio")
LiMCA_high_RME_ratio2p$dscore_type<-factor(LiMCA_high_RME_ratio2p$dscore_type,levels=c("Paternal_specific","Maternal_specific","Biallelic"))

LiMCA_high_RME_ratio_p<-ggplot(LiMCA_high_RME_ratio2p, aes(x = gene, y = ratio*100, fill = dscore_type)) +
  geom_bar(stat = "identity", position = "fill") +  # 100% 堆积柱状图
  scale_y_continuous(labels = scales::percent) +  # Y轴转换为百分比
  scale_fill_manual(values=c("#94b2c2","#c68989","#dfe2d6"))+
  theme_classic() +
  labs(x = "Gene", y = "Proportion", fill = "type")+
  theme(text=element_text(size=10),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
    axis.text.y=element_text(size=10,color = "black"),
    axis.text.x=element_text(size=10,  color = "black",angle = 45, hjust = 1))  
ggsave(filename = "~/F1_OSN/analysis/figure/LiMCA/LiMCA_high_RME_rank.pdf",LiMCA_high_RME_count_p+LiMCA_high_RME_ratio_p+plot_layout(heights=c(1,1.5)), device = "pdf", width = 24, height = 14, units = "cm", dpi = 300,bg="white")


LiMCA_high_RME_dscore<-LiMCA_gene_dscore_10%>%subset(gene%in%LiMCA_high_RME_ratio$gene)
LiMCA_high_RME_dscore$cell <- str_remove(LiMCA_high_RME_dscore$barcode, "_rna\\d+$")
LiMCA_high_RME_dscore$barcode <- str_replace(LiMCA_high_RME_dscore$barcode, "_rna\\d+$", "_rna")

LiMCA2Hic<-LiMCA_high_RME_count%>%subset(Maternal_specific>=5&Paternal_specific>=5)
LiMCA2Hic_dscore<-LiMCA_high_RME_dscore%>%subset(gene%in%LiMCA2Hic$gene)
# r$> LiMCA2Hic
# # A tibble: 11 × 6
#    type   gene    count Maternal_specific Paternal_specific Biallelic
#    <fct>  <chr>   <int>             <int>             <int>     <int>
#  1 Others Camkmt     14                 5                 9         0
#  2 Others Clec16a    28                11                10         7
#  3 Others Clip4      63                24                20        19
#  4 Others Dmxl2      47                20                13        14
#  5 Others Epha7      22                 6                10         6
#  6 Others Hecw1      28                 9                 7        12
#  7 Others Macrod2    41                19                16         6
#  8 Others Megf9      32                 6                15        11
#  9 Others Sorbs1     52                20                15        17
# 10 Others Wdr17      35                 5                20        10
# 11 Others Wwox       30                 5                16         9
write.table(LiMCA_high_RME_dscore, 
            file = "~/F1_OSN/LiMCA/HIC/LiMCA_high_RME_dscore.txt",
            quote = FALSE, 
            row.names = FALSE,sep="\t")
write.table(LiMCA2Hic_dscore, 
            file = "~/F1_OSN/LiMCA/HIC/LiMCA2Hic_dscore.txt",
            quote = FALSE, 
            row.names = FALSE,sep="\t")
LiMCA2Hic_dscore<-read.table( "~/F1_OSN/LiMCA/HIC/LiMCA2Hic_dscore.txt",header=TRUE)
load(file="~/F1_OSN/analysis/merge/RData/max_promoter_data.RData")

LiMCA2Hic_promoter<-max_promoter_data%>%subset(gene%in%LiMCA2Hic_dscore$gene)
write.table(LiMCA2Hic_promoter, 
            file = "~/F1_OSN/LiMCA/HIC/LiMCA2Hic_promoter.txt",
            quote = FALSE, 
            row.names = FALSE,sep="\t")
intra_ints <- read.table("~/F1_OSN/LiMCA/HIC/bulk/fithic_5k/out/intraOnly/FitHiC.spline_pass1.res5000.q0.05.txt",header=TRUE)
inter_ints <- read.table("~/F1_OSN/LiMCA/HIC/bulk/fithic_5k/out/interOnly/FitHiC.spline_pass1.res5000.q0.05.txt",header=TRUE)

find_promoter_interactions <- function(ints, promoters, bin) {
  # 初始化结果
  all_results <- data.frame()
  
  for(i in 1:nrow(promoters)) {
    # 当前promoter
    p <- promoters[i, ]
    
    # 计算anchor1和anchor2的区间
    ints$anchor1_start <- ints$fragmentMid1 - bin/2
    ints$anchor1_end <- ints$fragmentMid1 + bin/2
    ints$anchor2_start <- ints$fragmentMid2 - bin/2
    ints$anchor2_end <- ints$fragmentMid2 + bin/2
    
    # 筛选：区间有重叠（不是包含关系）
    matches <- subset(ints,
      (chr1 == p$chr & 
       !(anchor1_end < p$start | anchor1_start > p$end)) |  # anchor1与promoter重叠
      (chr2 == p$chr & 
       !(anchor2_end < p$start | anchor2_start > p$end))    # anchor2与promoter重叠
    )
    # matches <- subset(ints,
    #   (chr1 == p$chr & fragmentMid1 < p$end & fragmentMid1 > p$start) |  # fragment1在promoter中间
    #   (chr2 == p$chr & fragmentMid2 < p$end & fragmentMid2 > p$start)    # fragment2在promoter中间
    # )
    if(nrow(matches) > 0) {
      matches$gene <- p$gene
      matches$promoter_start<-p$start
      matches$promoter_end<-p$end
      all_results <- rbind(all_results, matches)
    }
  }
  
  return(all_results)
}

# 使用
intra_RME_ints <- find_promoter_interactions(intra_ints, LiMCA2Hic_promoter,5000)
inter_RME_ints <- find_promoter_interactions(inter_ints, LiMCA2Hic_promoter,5000)

library(tidyverse)

#
process_gene_simple <- function(gene, intra_RME_ints, base_dir) {
  promoter<-max_promoter_data%>%subset(gene==gene)
  gene_dir <- file.path(base_dir, gene)
  if (!dir.exists(gene_dir)) return(NULL)
  
  # 读取三种文件
  files <- c(
    maternal = "promoter_intra_Maternal_specific.impute.con.gz",
    paternal = "promoter_intra_Paternal_specific.impute.con.gz",
    biallelic = "promoter_intra_Biallelic.impute.con.gz"
  )
  
  # 读取数据
  data_list <- list()
  for (type in names(files)) {
    fpath <- file.path(gene_dir, files[[type]])
    if (file.exists(fpath)) {
      # 读取并解析 - 不创建allele_pair
      data <- read_delim(fpath, "\t", col_names = c("coord1", "coord2"))
      
      # 解析坐标，保持allele1和allele2分开
      parsed <- data %>%
        separate(coord1, c("chr1", "pos1", "allele1"), ",", convert = TRUE) %>%
        separate(coord2, c("chr2", "pos2", "allele2"), ",", convert = TRUE) %>%
        mutate(allelic_exp = type)
      
      data_list[[type]] <- parsed
    }
  }
  data_list<-bind_rows(data_list)

  # 筛选该基因的interactions
  gene_ints <- intra_RME_ints %>% subset(gene == gene)
  
  # 处理每个interaction
  results <- list()
  
  for (i in 1:nrow(gene_ints)) {
    int <- gene_ints[i, ]
    
    # 判断哪个anchor包含promoter
    # promoter在promoter_start和promoter_end之间
    promoter_in_anchor2 <- int$anchor1_start < int$promoter_end | int$anchor1_end > int$promoter_start
 
      # 在数据中查找匹配
    matched <- data_list %>%
        filter(
          pos1 >= int$anchor1_start & 
          pos1 <= int$anchor1_end &
          pos2 >= int$anchor2_start & 
          pos2 <= int$anchor2_end
        )

    adjusted <- matched
    if (!promoter_in_anchor2) {
        # 根据promoter位置调整allele顺序
          # promoter在anchor1，需要交换，让promoter的allele变成allele2
          adjusted <- matched %>%
            mutate(
              # 交换chr和pos
              temp_chr = chr1, temp_pos = pos1, temp_allele = allele1,
              chr1 = chr2, pos1 = pos2, allele1 = allele2,
              chr2 = temp_chr, pos2 = temp_pos, allele2 = temp_allele
            ) %>%
            select(-temp_chr, -temp_pos, -temp_allele)
    }
    adjusted$interaction_id<-i
    results[[i]] <- adjusted
  }
  
  # 合并所有结果
  if (length(results) > 0) {
    return(bind_rows(results))
  } else {
    return(NULL)
  }
}

# 批量处理
process_all_genes <- function(ints, base_dir) {
  # 读取interaction数据
  # 处理每个基因
  results <- list()
  genes <- unique(ints$gene)
  
  for (gene in genes) {
    cat("处理基因:", gene, "\n")
    res <- process_gene_simple(gene, ints, base_dir)
    res$gene<-gene
    if (!is.null(res)) results[[gene]] <- res
  }
  
  # 合并并保存结果
  final_results <- bind_rows(results)
  return(final_results)
}

# 使用
results <- process_all_genes(
  intra_RME_ints, 
  "~/F1_OSN/LiMCA/HIC/merge_impute"
)

# r$> results%>%group_by(gene,interaction_id)%>%summarise(count=n())
# `summarise()` has grouped output by 'gene'. You can override using the `.groups` argument.
# # A tibble: 13 × 3
# # Groups:   gene [4]
#    gene  interaction_id count
#    <chr>          <int> <int>
#  1 Hecw1              2     3
#  2 Hecw1              3     1
#  3 Megf9              2     1
#  4 Megf9              3     1
#  5 Megf9             11     1
#  6 Megf9             12     1
#  7 Megf9             13     1
#  8 Wdr17              2     1
#  9 Wdr17             13     1
# 10 Wdr17             19     1
# 11 Wwox               2     1
# 12 Wwox              16     2
# 13 Wwox              19     1


# TAD analysis
TAD<-read.table(file="~/F1_OSN/LiMCA/HIC/bulk/TAD/TAD_r5kb_contact_domains_list/5000_blocks.bedpe")
colnames(TAD) <- c("chr1", "x1", "x2", "chr2", "y1", "y2", "name", "score",
                    "strand1", "strand2", "color", "score", "uVarScore", 
                    "lVarScore", "upSign", "loSign")
TAD$chr1<-paste0("chr",TAD$chr1)
TAD$chr2<-paste0("chr",TAD$chr2)
promoter_TAD<- data.frame()

for(i in 1:nrow(LiMCA2Hic_promoter)) {
  gene_info <- LiMCA2Hic_promoter[i, ]
  
  # 直接匹配：染色体相同且promoter起始点在TAD区间内
  matching_tad <- TAD[TAD$chr1 == gene_info$chr & 
                      TAD$x1 <= gene_info$start & 
                      TAD$x2 >= gene_info$end, ]
  
  if(nrow(matching_tad) > 0) {
    for(j in 1:nrow(matching_tad)) {
      matching_tad$gene<-LiMCA2Hic_promoter[i,1]
      promoter_TAD <- rbind(promoter_TAD, matching_tad)
    }
  }
}

library(strawr)

# 处理单个基因 - 统计总接触数
process_gene_total_counts <- function(gene_row) {
  gene <- gene_row$gene
  chr <- gene_row$chr1
  tad_start <- gene_row$x1
  tad_end <- gene_row$x2
  
  # 去掉染色体名称中的"chr"前缀
  chr_no_prefix <- gsub("^chr", "", chr)
  
  # 基因目录
  gene_dir <- paste0("~/F1_OSN/LiMCA/HIC/merge_impute/", gene)
  
  # 创建结果行（包含TAD基本信息）
  result_row <- data.frame(
    gene = gene,
    chr = chr,
    TAD_start = tad_start,
    TAD_end = tad_end
  )
  
  # 处理所有样本类型和等位基因组合
  for (sample in c("Maternal_specific", "Paternal_specific")) {
    hic_file <- paste0(gene_dir, "/", sample, ".impute.hic")
    
    for (contact in c("MM", "PP", "MP", "PM")) {
      col_name <- paste0(sample, "_", contact)
      
      alleles <- switch(contact,
        "MM" = c("MAT", "MAT"),
        "PP" = c("PAT", "PAT"),
        "MP" = c("MAT", "PAT"),
        "PM" = c("PAT", "MAT")
      )
      
      if (file.exists(hic_file)) {
        # 提取接触矩阵
        contact_data <- straw(
          norm = "NONE",
          fname = hic_file,
          chr1loc = paste0(chr_no_prefix, "(", alleles[1], "):", tad_start, ":", tad_end),
          chr2loc = paste0(chr_no_prefix, "(", alleles[2], "):", tad_start, ":", tad_end),
          unit = "BP",
          binsize = 5000
        )
        
        # 统计总接触数
        result_row[[col_name]] <- sum(contact_data$counts)
      } else {
        result_row[[col_name]] <- NA
      }
    }
  }
  
  return(result_row)
}

# 批量处理所有基因
all_results <- list()

for (i in 1:nrow(promoter_TAD)) {
  cat("Processing gene", i, "/", nrow(promoter_TAD), "\n")
  gene_result <- process_gene_total_counts(promoter_TAD[i, ])
  all_results[[i]] <- gene_result
}

# 合并所有结果
final_df <- do.call(rbind, all_results)

# 查看结果
print(final_df)

# 保存结果
write.xlsx(final_df, "~/F1_OSN/LiMCA/HIC/merge_impute/TAD_total_contact_counts.xlsx", row.names = FALSE)

final_df$Maternal_specific_dscore<-final_df$Maternal_specific_PP/(final_df$Maternal_specific_MM+final_df$Maternal_specific_PP)-0.5
final_df$Paternal_specific_dscore<-final_df$Paternal_specific_PP/(final_df$Paternal_specific_MM+final_df$Paternal_specific_PP)-0.5
final_df<-final_df%>%mutate(total_counts=(Paternal_specific_PP+Paternal_specific_MM+Maternal_specific_PP+Maternal_specific_MM))
p <- ggplot(final_df, 
            aes(x = Paternal_specific_dscore, 
                y = Maternal_specific_dscore,
                size = total_counts)) +
  geom_rect(aes(xmin = 0.15, xmax = 0.5, ymin = -0.15, ymax = -0.5),
            fill = "lightblue", alpha = 0.1, inherit.aes = FALSE) +
  
  # 添加散点
  geom_point() +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  # 添加坐标轴交叉线（在0,0点）
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # 坐标轴标签
  labs(
    x = "Paternal Group TAD D-score",
    y = "Maternal Group TAD D-score",
    title = "TAD Interaction D-scores of RME"
  ) +
  
  # 使用theme_classic作为基础主题
  theme_classic() +
  
  # 调整主题设置
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    strip.background = element_rect(fill = "lightgray"),
    strip.text = element_text(size = 11),
    # 确保坐标轴交叉在0,0
    axis.line = element_line(color = "black", linewidth = 0.5),
    # 如果需要网格线，可以添加以下设置
    # panel.grid.major = element_line(color = "gray90", linewidth = 0.2),
    # panel.grid.minor = element_line(color = "gray95", linewidth = 0.1)
  )

ggsave(
  "~/F1_OSN/analysis/figure/LiMCA/TAD_dscore_scatterplot.pdf",
  p,
  width = 5.5,
  height = 4,
  dpi = 300
)

library(org.Mm.eg.db)
library("plotgardener")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plotgardener)
hicFile <- "~/F1_OSN/LiMCA/HIC/bulk/4DNFIVUGNDD7.hic"
library(plotgardener)
library(strawr)

# 主绘图函数
library(plotgardener)
library(strawr)

plot_gene_hic_comparison <- function(gene_row) {
  gene <- gene_row$gene
  chr <- gene_row$chr
  tad_start <- gene_row$TAD_start
  tad_end <- gene_row$TAD_end
  
  # 1. 创建PDF页面（增加高度）
  pageCreate(width = 7, height = 11, default.units = "inches", showGuides = FALSE)
  
  # 2. 扩展区间（向左右各扩展1/4）
  tad_length <- tad_end - tad_start
  extended_start <- max(0, tad_start - round(tad_length / 4))
  extended_end <- tad_end + round(tad_length / 4)
  
  # 设置参数
  params <- pgParams(
    chrom = chr,
    chromstart = extended_start,
    chromend = extended_end,
    assembly = "mm10",
    x = 0.5,
    width = 6,
    default.units = "inches"
  )
  
  # 3. 定义布局参数
  plot_width <- 6
  plot_x <- 0.5
  current_y <- 0.7  # 从更高的位置开始
  
  # 添加主标题
  plotText(
    label = paste("Gene:", gene, "(", chr, ":", 
                  format(tad_start, big.mark = ","), "-", 
                  format(tad_end, big.mark = ","), ")"),
    fontsize = 11,
    fontface = "bold",
    x = 3.5,  # 居中
    y = 0.35,
    default.units = "inches",
    just = "center"
  )
  
  # 4. 绘制Bulk Hi-C
  cat("Plotting bulk Hi-C for", gene, "...\n")
  
  # 子标题
  plotText(
    label = "Bulk Hi-C",
    fontsize = 10,
    fontface = "bold",
    x = plot_x,
    y = current_y - 0.2,
    default.units = "inches"
  )
  
  # 读取bulk Hi-C数据
  hicDataChrom <- readHic(
    file = hicFile,
    chrom = gsub("chr", "", chr),
    assembly = "mm10",
    resolution = 5000,
    res_scale = "BP",
    norm = "KR"
  )
  
  # 绘制三角形热图
  bulk_plot <- plotHicTriangle(
    data = hicDataChrom,
    params = params,
    zrange = c(0, 30),
    resolution = 5000,
    x = plot_x,
    y = current_y,
    width = plot_width,
    height = 2.0,
    default.units = "inches"
  )
  
  # 添加图例
  annoHeatmapLegend(
    plot = bulk_plot,
    fontsize = 8,
    x = plot_x + plot_width + 0.05,
    y = current_y,
    width = 0.08,
    height = 0.7,
    just = c("left", "top"),
    default.units = "inches"
  )
  
  current_y <- current_y + 2.2  # Bulk Hi-C后的间距
  
  # 添加分隔线
  plotSegments(
    x0 = plot_x, y0 = current_y,
    x1 = plot_x + plot_width, y1 = current_y,
    lty = "dashed", lwd = 0.8,
    default.units = "inches"
  )
  
  current_y <- current_y + 0.3  # 分隔线后的间距
  
  # 5. 绘制Allele-specific Hi-C
  gene_dir <- paste0("~/F1_OSN/LiMCA/HIC/merge_impute/", gene)
  
  # 样本类型和对应的高度
  sample_heights <- list(
    "Maternal_specific" = current_y,
    "Paternal_specific" = current_y + 3.2  # 为第二个样本留出更多空间
  )
  
  for (sample in c("Maternal_specific", "Paternal_specific")) {
    hic_file <- paste0(gene_dir, "/", sample, ".impute.hic")
    sample_y <- sample_heights[[sample]]
    
    if (file.exists(hic_file)) {
      cat("Processing", sample, "for", gene, "...\n")
      
      # 样本标题（更大间距）
      plotText(
        label = sample,
        fontsize = 10,
        fontface = "bold",
        x = plot_x,
        y = sample_y - 0.25,
        default.units = "inches"
      )
      
      # 提取MM和PP接触矩阵
      contact_types <- c("MM", "PP")
      plot_height <- 1.2  # 减小单个图高度
      vertical_spacing <- 0.3  # 增加图之间的垂直间距
      
      for (contact_idx in 1:2) {
        contact_type <- contact_types[contact_idx]
        alleles <- switch(contact_type,
          "MM" = c("MAT", "MAT"),
          "PP" = c("PAT", "PAT")
        )
        
        # 使用straw提取数据
        chr_no_prefix <- gsub("chr", "", chr)
        contact_data <- straw(
          norm = "NONE",
          fname = hic_file,
          chr1loc = paste0(chr_no_prefix, "(", alleles[1], "):", extended_start, ":", extended_end),
          chr2loc = paste0(chr_no_prefix, "(", alleles[2], "):", extended_start, ":", extended_end),
          unit = "BP",
          binsize = 5000
        )
        
        # 重命名列
        colnames(contact_data) <- c(paste0(chr_no_prefix,"_A"), paste0(chr_no_prefix,"_B"), "counts")
        
        # 计算当前接触类型的y位置（增加间距）
        plot_y <- sample_y + (contact_idx - 1) * (plot_height + vertical_spacing)
        
        # 绘制三角形热图
        as_plot <- plotHicTriangle(
          data = contact_data,
          params = params,
          zrange = c(0, 2),
          resolution = 5000,
          x = plot_x,
          y = plot_y,
          width = plot_width,
          height = plot_height,
          default.units = "inches"
        )
        
        # 添加图例
        annoHeatmapLegend(
          plot = as_plot,
          fontsize = 7,
          x = plot_x + plot_width + 0.05,
          y = plot_y,
          width = 0.07,
          height = 0.5,
          just = c("left", "top"),
          default.units = "inches"
        )
        
        # 添加接触类型标签（放在图左侧）
        plotText(
          label = contact_type,
          fontsize = 9,
          fontface = "bold",
          x = plot_x - 0.25,
          y = plot_y + plot_height/2,
          default.units = "inches",
          just = c("right", "center")
        )
      }
      
      # 样本内部的水平分隔线
      if (sample == "Maternal_specific") {
        plotSegments(
          x0 = plot_x, y0 = sample_y + 2*(plot_height + vertical_spacing) - 0.15,
          x1 = plot_x + plot_width, y1 = sample_y + 2*(plot_height + vertical_spacing) - 0.15,
          lty = "dotted", lwd = 0.5,
          default.units = "inches"
        )
      }
      
    } else {
      cat("File not found:", hic_file, "\n")
      # 如果没有文件，添加占位符
      plotText(
        label = paste("No", sample, "data"),
        fontsize = 9,
        fontcolor = "gray50",
        x = plot_x,
        y = sample_y,
        default.units = "inches"
      )
    }
  }
  
  # 6. 绘制基因注释（放在最下面，留出足够空间）
  genes_y <- current_y + 6.5  # 为allele-specific plots留出更多空间
  
  genes_plot <- plotGenes(
    params = params,
    stroke = 0.8,
    fontsize = 8,  # 增加字体大小
    strandLabels = FALSE,
    x = plot_x,
    y = genes_y,
    width = plot_width,
    height = 0.6,  # 增加基因注释高度
    default.units = "inches"
  )
  
  # 7. 添加基因组坐标标签
  annoGenomeLabel(
    plot = genes_plot,
    params = params,
    scale = "Kb",
    fontsize = 9,
    x = plot_x,
    y = genes_y + 0.65,
    width = plot_width,
    default.units = "inches"
  )
}

# 批量处理函数
plot_all_genes_hic <- function(final_df) {
  # 创建PDF文件（增加高度）
  pdf("~/F1_OSN/LiMCA/HIC/merge_impute/gene_hic_comparison.pdf", 
      width = 7, height = 11)
  
  # 遍历每个基因
  for(i in 1:nrow(final_df)) {
    cat("\n=== Processing gene", i, "of", nrow(final_df), ":",
        final_df$gene[i], "===\n")
    
    plot_gene_hic_comparison(final_df[i, ])
  }
  
  dev.off()
  cat("\nPDF saved with", nrow(final_df), "pages\n")
}
plot_all_genes_hic(final_df)


# 创建PDF文件
pdf("~/F1_OSN/LiMCA/HIC/merge_impute/TAD_plots.pdf", 
    width = 3.25, height = 3.5)
TAD_max<-c(15,15,30,35)
# 遍历final_df的每一行
for(i in 1:nrow(final_df)) {
  
  # 提取当前行的信息
  current_gene <- final_df$gene[i]
  current_chr <- final_df$chr[i]
  current_start <- final_df$TAD_start[i]
  current_end <- final_df$TAD_end[i]
  
  # 计算TAD长度
  tad_length <- current_end - current_start
  
  # 向左右各扩展三分之一的TAD长度
  extended_start <- max(0, current_start - round(tad_length / 4))
  extended_end <- current_end + round(tad_length / 4)
  
  cat(sprintf("Processing %s: %s:%s-%s (extended to %s-%s)\n", 
              current_gene, current_chr, 
              format(current_start, big.mark = ","), 
              format(current_end, big.mark = ","),
              format(extended_start, big.mark = ","),
              format(extended_end, big.mark = ",")))
  
  # 为每个TAD创建新页面
  pageCreate(width = 3.25, height = 3.5, default.units = "inches", 
             showGuides = FALSE)
  
  # 读取Hi-C数据
  chrom_num <- gsub("chr", "", current_chr)
  
  hicDataChrom <- readHic(file = hicFile,
      chrom = chrom_num, 
      assembly = "mm10",
      resolution = 5000, 
      res_scale = "BP", 
      norm = "KR"
  )
  
  # 设置参数 - 使用扩展后的区间
  params_a <- pgParams(
      chrom = current_chr,
      chromstart = extended_start,  # 使用扩展后的起始位置
      chromend = extended_end,      # 使用扩展后的结束位置
      assembly = "mm10",
      x = 0.25, 
      width = 2.5, 
      default.units = "inches"
  )
  
  # 创建三角形Hi-C图
  hicPlot_top <- plotHicTriangle(
      data = hicDataChrom, 
      params = params_a,
      zrange = c(0, TAD_max[i]),
      resolution = 5000,
      x = 0.25, 
      y = 0.5,
      width = 2.5,
      height = 2.0,
      default.units = "inches"
  )
  
  # 添加图例
  annoHeatmapLegend(
      plot = hicPlot_top, 
      fontsize = 6,
      x = 2.8,
      y = 0.5, 
      width = 0.06, 
      height = 0.4,
      just = c("right", "top"), 
      default.units = "inches"
  )
  
  # 绘制基因
  genes_gm <- plotGenes(
      params = params_a, 
      stroke = 0.5,
      fontsize = 5,
      strandLabels = FALSE,
      y = 2.6,
      height = 0.3
  )
  
  # 添加基因组标签
  annoGenomeLabel(
      plot = genes_gm, 
      params = params_a, 
      scale = "Kb", 
      fontsize = 6,
      y = 3.0
  )
  
  # 添加基因名称和位置信息
  plotText(
    label = paste0(current_gene, "\n", 
                   current_chr, ":", 
                   format(current_start, big.mark = ","), "-", 
                   format(current_end, big.mark = ","),
                   "\n(Extended: ", 
                   format(extended_start, big.mark = ","), "-",
                   format(extended_end, big.mark = ","), ")"),
    fontsize = 6,
    x = 0.25,
    y = 0.25,
    default.units = "inches",
    just = c("left", "top")
  )
}

# 关闭PDF设备
dev.off()

# LiMCA2hic_dscore_RNA vln plot

# 获取所有唯一的基因名
unique_genes <- unique(final_df$gene)

output_pdf <- "~/F1_OSN/LiMCA/HIC/merge_impute/all_genes_RNA_expression_log2.pdf"

# 打开PDF设备
pdf(output_pdf, width = 12/2.54, height = 9/2.54)  # 转换为英寸 (cm/2.54)

cat("Creating PDF file:", output_pdf, "\n")

# 为每个基因创建图并添加到PDF
for(gene_name in unique_genes) {
  cat("Processing gene:", gene_name, "\n")
  
  # 筛选当前基因的数据
  gene_data <- LiMCA2Hic_dscore %>%
    filter(gene == gene_name) %>%
    # 创建log2(x+1)转换后的列
    mutate(
      log2_G1rna = log2(G1rna + 1),  # log2(x+1)转换
      log2_G2rna = log2(G2rna + 1)   # log2(x+1)转换
    )
  
  # 将数据转换为长格式
  gene_data_long <- gene_data %>%
    pivot_longer(
      cols = c(log2_G1rna, log2_G2rna),
      names_to = "RNA_type",
      values_to = "log2_RNA_expression"
    ) %>%
    mutate(
      RNA_type = case_when(
        RNA_type == "log2_G1rna" ~ "G1rna",
        RNA_type == "log2_G2rna" ~ "G2rna",
        TRUE ~ RNA_type
      )
    )
  
  # 创建violin plot，使用指定的颜色
  plot <- ggviolin(
    gene_data_long, 
    x = "rna_dscore_type", 
    y = "log2_RNA_expression", 
    fill = "RNA_type",
    palette = c("#505f84", "#b45a55"),  # 使用指定的颜色
    alpha = 0.8,
    width = 0.7,
    position = position_dodge(0.8),
    legend = "top",
    xlab = "",
    ylab = "RNA Expression Level [log2(x+1)]",
    title = gene_name,
    font.title = c(14, "bold"),
    font.tickslab = c(11, "plain", "black"),
    add.params = list(
      fill = "white", 
      width = 0.12, 
      linetype = 1,
      position = position_dodge(0.8)
    )
  ) +
    theme_classic() + 
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11),
      axis.text.x = element_text(
        size = 11, 
        angle = 0, 
        hjust = 0.5,
        color = "black"
      ),
      axis.title.y = element_text(size = 12),
      axis.text.y = element_text(size = 10),
      panel.grid.major.y = element_line(color = "gray90", size = 0.2),
      panel.grid.minor.y = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.1))
    )
  
  # 将图添加到PDF
  print(plot)
}

# 关闭PDF设备
dev.off()

cat("\nPDF saved to:", output_pdf, "\n")
cat("Total genes processed:", length(unique_genes), "\n")
cat("Total pages in PDF:", length(unique_genes), "\n")


plot_gene_hic_comparison_with_bia <- function(gene_row) {
  gene <- gene_row$gene
  chr <- gene_row$chr
  tad_start <- gene_row$TAD_start
  tad_end <- gene_row$TAD_end
  
  # 1. 创建PDF页面（增加高度以容纳biallelic部分）
  pageCreate(width = 7, height = 14, default.units = "inches", showGuides = FALSE)
  
  # 2. 扩展区间（向左右各扩展1/5）
  tad_length <- tad_end - tad_start
  extended_start <- max(0, tad_start - round(tad_length / 5))
  extended_end <- tad_end + round(tad_length / 4)
  
  # 设置参数
  params <- pgParams(
    chrom = chr,
    chromstart = extended_start,
    chromend = extended_end,
    assembly = "mm10",
    x = 0.5,
    width = 6,
    default.units = "inches"
  )
  
  # 3. 定义布局参数
  plot_width <- 6
  plot_x <- 0.5
  current_y <- 0.7  # 从更高的位置开始
  
  # 添加主标题
  plotText(
    label = paste("Gene:", gene, "(", chr, ":", 
                  format(tad_start, big.mark = ","), "-", 
                  format(tad_end, big.mark = ","), ")"),
    fontsize = 11,
    fontface = "bold",
    x = 3.5,  # 居中
    y = 0.35,
    default.units = "inches",
    just = "center"
  )
  
  # 4. 绘制Bulk Hi-C
  cat("Plotting bulk Hi-C for", gene, "...\n")
  
  # 子标题
  plotText(
    label = "Bulk Hi-C",
    fontsize = 10,
    fontface = "bold",
    x = plot_x,
    y = current_y - 0.2,
    default.units = "inches"
  )
  
  # 读取bulk Hi-C数据
  hicDataChrom <- readHic(
    file = hicFile,
    chrom = gsub("chr", "", chr),
    assembly = "mm10",
    resolution = 5000,
    res_scale = "BP",
    norm = "KR"
  )
  
  # 绘制三角形热图
  bulk_plot <- plotHicTriangle(
    data = hicDataChrom,
    params = params,
    zrange = c(0, 30),
    resolution = 5000,
    x = plot_x,
    y = current_y,
    width = plot_width,
    height = 1.8,  # 稍微减小高度
    default.units = "inches"
  )
  
  # 添加图例
  annoHeatmapLegend(
    plot = bulk_plot,
    fontsize = 8,
    x = plot_x + plot_width + 0.05,
    y = current_y,
    width = 0.08,
    height = 0.6,
    just = c("left", "top"),
    default.units = "inches"
  )
  
  current_y <- current_y + 2.0  # Bulk Hi-C后的间距
  
  # 添加分隔线
  plotSegments(
    x0 = plot_x, y0 = current_y,
    x1 = plot_x + plot_width, y1 = current_y,
    lty = "dashed", lwd = 0.8,
    default.units = "inches"
  )
  
  current_y <- current_y + 0.3  # 分隔线后的间距
  
  # 5. 绘制Allele-specific和Biallelic Hi-C
  gene_dir <- paste0("~/F1_OSN/LiMCA/HIC/merge_impute/", gene)
  
  # 样本列表
  samples <- c("Maternal_specific", "Paternal_specific", "Biallelic")
  
  # 遍历所有样本
  for (sample_idx in 1:length(samples)) {
    sample <- samples[sample_idx]
    hic_file <- paste0(gene_dir, "/", sample, ".impute.hic")
    
    # 计算当前样本的y位置
    if (sample_idx == 1) {
      sample_y <- current_y
    } else {
      sample_y <- current_y + (sample_idx - 1) * 2.8
    }
    
    cat("Processing", sample, "for", gene, "...\n")
    
    if (file.exists(hic_file)) {
      # 样本标题
      plotText(
        label = sample,
        fontsize = 9,
        fontface = "bold",
        x = plot_x,
        y = sample_y - 0.2,
        default.units = "inches"
      )
      
      # 提取MM和PP接触矩阵（所有样本都一样）
      contact_types <- c("MM", "PP")
      plot_height <- 1.0  # 减小单个图高度
      vertical_spacing <- 0.25  # 减小间距
      
      for (contact_idx in 1:2) {
        contact_type <- contact_types[contact_idx]
        alleles <- switch(contact_type,
          "MM" = c("MAT", "MAT"),
          "PP" = c("PAT", "PAT")
        )
        
        # 使用straw提取数据
        chr_no_prefix <- gsub("chr", "", chr)
        contact_data <- straw(
          norm = "NONE",
          fname = hic_file,
          chr1loc = paste0(chr_no_prefix, "(", alleles[1], "):", extended_start, ":", extended_end),
          chr2loc = paste0(chr_no_prefix, "(", alleles[2], "):", extended_start, ":", extended_end),
          unit = "BP",
          binsize = 5000
        )
        
        # 重命名列
        colnames(contact_data) <- c(paste0(chr_no_prefix,"_A"), paste0(chr_no_prefix,"_B"), "counts")
        
        # 计算当前接触类型的y位置
        plot_y <- sample_y + (contact_idx - 1) * (plot_height + vertical_spacing)
        
        # 绘制三角形热图（所有样本使用相同的zrange）
        as_plot <- plotHicTriangle(
          data = contact_data,
          params = params,
          zrange = c(0, 2),  # 所有样本使用相同的颜色范围
          resolution = 5000,
          x = plot_x,
          y = plot_y,
          width = plot_width,
          height = plot_height,
          default.units = "inches"
        )
        
        # 添加图例
        annoHeatmapLegend(
          plot = as_plot,
          fontsize = 7,
          x = plot_x + plot_width + 0.05,
          y = plot_y,
          width = 0.06,
          height = 0.4,
          just = c("left", "top"),
          default.units = "inches"
        )
        
        # 添加接触类型标签
        plotText(
          label = contact_type,
          fontsize = 8,
          fontface = "bold",
          x = plot_x - 0.2,
          y = plot_y + plot_height/2,
          default.units = "inches",
          just = c("right", "center")
        )
      }
      
      # 样本内部的分隔线
      plotSegments(
        x0 = plot_x, y0 = sample_y + 2*(plot_height + vertical_spacing) - 0.1,
        x1 = plot_x + plot_width, y1 = sample_y + 2*(plot_height + vertical_spacing) - 0.1,
        lty = "dotted", lwd = 0.5,
        default.units = "inches"
      )
      
    } else {
      cat("File not found:", hic_file, "\n")
      # 如果没有文件，添加占位符
      plotText(
        label = paste("No", sample, "data"),
        fontsize = 8,
        fontcolor = "gray50",
        x = plot_x,
        y = sample_y,
        default.units = "inches"
      )
    }
  }
  
  # 6. 绘制基因注释（放在最下面，留出足够空间）
  genes_y <- current_y + 8.4  # 为所有三个样本留出空间
  
  genes_plot <- plotGenes(
    params = params,
    stroke = 0.8,
    fontsize = 8,
    strandLabels = TRUE,  # 显示正负链标签
    x = plot_x,
    y = genes_y,
    width = plot_width,
    height = 0.6,
    default.units = "inches"
  )
  
  # 7. 添加基因组坐标标签
  annoGenomeLabel(
    plot = genes_plot,
    params = params,
    scale = "Kb",
    fontsize = 9,
    x = plot_x,
    y = genes_y + 0.65,
    width = plot_width,
    default.units = "inches"
  )
}

# 批量处理函数（需要调整PDF高度）
plot_all_genes_hic_with_bia <- function(final_df) {
  # 创建PDF文件（增加高度以容纳biallelic部分）
  pdf("~/F1_OSN/LiMCA/HIC/merge_impute/gene_hic_comparison_with_bia.pdf", 
      width = 7, height = 14)  # 高度增加到14英寸
  
  # 遍历每个基因
  for(i in 1:nrow(final_df)) {
    cat("\n=== Processing gene", i, "of", nrow(final_df), ":",
        final_df$gene[i], "===\n")
    
    plot_gene_hic_comparison_with_bia(final_df[i, ])
  }
  
  dev.off()
  cat("\nPDF saved with", nrow(final_df), "pages\n")
}

# 使用函数
plot_all_genes_hic_with_bia(final_df[c(1,2,4),])

# 如果每个基因一个单独的pdf，使用以下函数
plot_gene_hic_comparison_separate <- function(gene_row) {
  gene <- gene_row$gene
  chr <- gene_row$chr
  tad_start <- gene_row$TAD_start
  tad_end <- gene_row$TAD_end
  
  # 扩展区间（向左右各扩展1/5）
  tad_length <- tad_end - tad_start
  extended_start <- max(0, tad_start - round(tad_length / 5))
  extended_end <- tad_end + round(tad_length / 5)
  
  # 设置公共参数
  params <- pgParams(
    chrom = chr,
    chromstart = extended_start,
    chromend = extended_end,
    assembly = "mm10",
    x = 0.5,
    width = 6,
    default.units = "inches"
  )
  
  plot_width <- 6
  plot_x <- 0.5
  
  # 1. 创建PDF文件，每个基因一个文件
  pdf_file <- paste0("~/F1_OSN/LiMCA/HIC/merge_impute/", 
                     gene, "_hic_comparison.pdf")
  pdf(pdf_file, width = 7, height = 11)
  
  # ========== 第1页: Bulk Hi-C ==========
  cat("Plotting bulk Hi-C for", gene, "...\n")
  pageCreate(width = 7, height = 11, default.units = "inches", showGuides = FALSE)
  
  # 主标题
  plotText(
    label = paste("Gene:", gene, "(", chr, ":", 
                  format(tad_start, big.mark = ","), "-", 
                  format(tad_end, big.mark = ","), ")"),
    fontsize = 11,
    fontface = "bold",
    x = 3.5,
    y = 0.35,
    default.units = "inches",
    just = "center"
  )
  
  # 子标题
  plotText(
    label = "Bulk Hi-C",
    fontsize = 12,
    fontface = "bold",
    x = 3.5,
    y = 0.7,
    default.units = "inches",
    just = "center"
  )
  
  current_y <- 1.0
  
  # 读取bulk Hi-C数据
  hicDataChrom <- readHic(
    file = hicFile,
    chrom = gsub("chr", "", chr),
    assembly = "mm10",
    resolution = 5000,
    res_scale = "BP",
    norm = "KR"
  )
  
  # 绘制三角形热图
  bulk_plot <- plotHicTriangle(
    data = hicDataChrom,
    params = params,
    zrange = c(0, 30),
    resolution = 5000,
    x = plot_x,
    y = current_y,
    width = plot_width,
    height = 3.0,  # 增加高度
    default.units = "inches"
  )
  
  # 添加图例
  annoHeatmapLegend(
    plot = bulk_plot,
    fontsize = 9,
    x = plot_x + plot_width + 0.05,
    y = current_y,
    width = 0.1,
    height = 1.0,
    just = c("left", "top"),
    default.units = "inches"
  )
  
  # 添加基因注释
  genes_y <- current_y + 3.2
  genes_plot <- plotGenes(
    params = params,
    stroke = 0.8,
    fontsize = 9,
    strandLabels = FALSE,
    x = plot_x,
    y = genes_y,
    width = plot_width,
    height = 0.7,
    default.units = "inches"
  )
  
  # 添加基因组坐标标签
  annoGenomeLabel(
    plot = genes_plot,
    params = params,
    scale = "Kb",
    fontsize = 10,
    x = plot_x,
    y = genes_y + 0.75,
    width = plot_width,
    default.units = "inches"
  )
  
  # ========== 第2页: Paternal-specific ==========
  pageCreate(width = 7, height = 11, default.units = "inches", showGuides = FALSE)
  
  # 主标题
  plotText(
    label = paste("Gene:", gene, "(", chr, ":", 
                  format(tad_start, big.mark = ","), "-", 
                  format(tad_end, big.mark = ","), ")"),
    fontsize = 11,
    fontface = "bold",
    x = 3.5,
    y = 0.35,
    default.units = "inches",
    just = "center"
  )
  
  # 子标题
  plotText(
    label = "Paternal-specific Hi-C",
    fontsize = 12,
    fontface = "bold",
    x = 3.5,
    y = 0.7,
    default.units = "inches",
    just = "center"
  )
  
  current_y <- 1.0
  
  # 绘制Paternal-specific Hi-C
  plot_allele_specific_hic(gene, chr, extended_start, extended_end, 
                           "Paternal_specific", current_y, plot_x, plot_width)
  
  # ========== 第3页: Maternal-specific ==========
  pageCreate(width = 7, height = 11, default.units = "inches", showGuides = FALSE)
  
  # 主标题
  plotText(
    label = paste("Gene:", gene, "(", chr, ":", 
                  format(tad_start, big.mark = ","), "-", 
                  format(tad_end, big.mark = ","), ")"),
    fontsize = 11,
    fontface = "bold",
    x = 3.5,
    y = 0.35,
    default.units = "inches",
    just = "center"
  )
  
  # 子标题
  plotText(
    label = "Maternal-specific Hi-C",
    fontsize = 12,
    fontface = "bold",
    x = 3.5,
    y = 0.7,
    default.units = "inches",
    just = "center"
  )
  
  current_y <- 1.0
  
  # 绘制Maternal-specific Hi-C
  plot_allele_specific_hic(gene, chr, extended_start, extended_end, 
                           "Maternal_specific", current_y, plot_x, plot_width)
  
  # 关闭PDF设备
  dev.off()
  
  cat("Saved to:", pdf_file, "\n")
}

# 辅助函数：绘制allele-specific Hi-C
plot_allele_specific_hic <- function(gene, chr, extended_start, extended_end, 
                                     sample_type, start_y, plot_x, plot_width) {
  
  gene_dir <- paste0("~/F1_OSN/LiMCA/HIC/merge_impute/", gene)
  hic_file <- paste0(gene_dir, "/", sample_type, ".impute.hic")
  
  if (!file.exists(hic_file)) {
    plotText(
      label = paste("No", sample_type, "data available"),
      fontsize = 10,
      fontcolor = "gray50",
      x = 3.5,
      y = 3.5,
      default.units = "inches",
      just = "center"
    )
    return()
  }
  
  cat("Processing", sample_type, "for", gene, "...\n")
  
  current_y <- start_y
  plot_height <- 1.5  # 增加单个图高度
  vertical_spacing <- 0.4  # 增加间距
  
  # 提取MM和PP接触矩阵
  contact_types <- c("MM", "PP")
  
  for (contact_idx in 1:2) {
    contact_type <- contact_types[contact_idx]
    alleles <- switch(contact_type,
      "MM" = c("MAT", "MAT"),
      "PP" = c("PAT", "PAT")
    )
    
    # 使用straw提取数据
    chr_no_prefix <- gsub("chr", "", chr)
    contact_data <- straw(
      norm = "NONE",
      fname = hic_file,
      chr1loc = paste0(chr_no_prefix, "(", alleles[1], "):", extended_start, ":", extended_end),
      chr2loc = paste0(chr_no_prefix, "(", alleles[2], "):", extended_start, ":", extended_end),
      unit = "BP",
      binsize = 5000
    )
    
    # 重命名列
    colnames(contact_data) <- c(paste0(chr_no_prefix,"_A"), paste0(chr_no_prefix,"_B"), "counts")
    
    # 计算y位置
    plot_y <- current_y + (contact_idx - 1) * (plot_height + vertical_spacing)
    
    # 绘制三角形热图
    params <- pgParams(
      chrom = chr,
      chromstart = extended_start,
      chromend = extended_end,
      assembly = "mm10",
      x = plot_x,
      width = plot_width,
      default.units = "inches"
    )
    
    as_plot <- plotHicTriangle(
      data = contact_data,
      params = params,
      zrange = c(0, 2),
      resolution = 5000,
      x = plot_x,
      y = plot_y,
      width = plot_width,
      height = plot_height,
      default.units = "inches"
    )
    
    # 添加图例
    annoHeatmapLegend(
      plot = as_plot,
      fontsize = 8,
      x = plot_x + plot_width + 0.05,
      y = plot_y,
      width = 0.08,
      height = 0.7,
      just = c("left", "top"),
      default.units = "inches"
    )
    
    # 添加接触类型标签
    plotText(
      label = paste(sample_type, contact_type),
      fontsize = 10,
      fontface = "bold",
      x = plot_x - 0.3,
      y = plot_y + plot_height/2,
      default.units = "inches",
      just = c("right", "center")
    )
  }
  
  # 添加基因注释
  genes_y <- current_y + 2*(plot_height + vertical_spacing) + 0.3
  params <- pgParams(
    chrom = chr,
    chromstart = extended_start,
    chromend = extended_end,
    assembly = "mm10",
    x = plot_x,
    width = plot_width,
    default.units = "inches"
  )
  
  genes_plot <- plotGenes(
    params = params,
    stroke = 0.8,
    fontsize = 9,
    strandLabels = FALSE,
    x = plot_x,
    y = genes_y,
    width = plot_width,
    height = 0.7,
    default.units = "inches"
  )
  
  # 添加基因组坐标标签
  annoGenomeLabel(
    plot = genes_plot,
    params = params,
    scale = "Kb",
    fontsize = 10,
    x = plot_x,
    y = genes_y + 0.75,
    width = plot_width,
    default.units = "inches"
  )
}

# 批量处理函数
plot_all_genes_separate <- function(final_df) {
  # 遍历每个基因
  for(i in 1:nrow(final_df)) {
    cat("\n=== Processing gene", i, "of", nrow(final_df), ":",
        final_df$gene[i], "===\n")
    
    plot_gene_hic_comparison_separate(final_df[i, ])
  }
  
  cat("\nAll PDFs saved in:", "~/F1_OSN/LiMCA/HIC/merge_impute/\n")
}

# 使用函数
plot_all_genes_separate(final_df)



# 创建PDF和页面 - 一步完成
pdf("~/F1_OSN/LiMCA/HIC/merge_impute/TAD_bulk.pdf", 
    width = 3.25, height = 3.5)

# 必须首先创建plotgardener页面
pageCreate(width = 3.25, height = 3.5, default.units = "inches", 
           showGuides = FALSE)  # 直接隐藏辅助线


# 读取Hi-C数据
hicDataChrom <- readHic(file = hicFile,
    chrom = "18", 
    assembly = "mm10",
    resolution = 5000, 
    res_scale = "BP", 
    norm = "KR"
)

# 设置参数
params_a <- pgParams(
    chrom = "chr18", 
    chromstart = 37699000, 
    chromend = 37979000, 
    assembly = "mm10",
    x = 0.25, 
    width = 2.5, 
    default.units = "inches"
)

# 创建三角形Hi-C图
hicPlot_top <- plotHicTriangle(
    data = hicDataChrom, 
    params = params_a,
    zrange = c(0, 30), 
    resolution = 5000,
    x = 0.25, 
    y = 0.5,
    width = 2.5,
    height = 2.0,
    default.units = "inches"
)

# 添加图例
annoHeatmapLegend(
    plot = hicPlot_top, 
    fontsize = 6,
    x = 2.8,
    y = 0.5, 
    width = 0.06, 
    height = 0.4,
    just = c("right", "top"), 
    default.units = "inches"
)

# 绘制基因
genes_gm <- plotGenes(
    params = params_a, 
    stroke = 0.5,
    fontsize = 5,
    strandLabels = FALSE,
    y = 2.6,
    height = 0.3
)

# 添加基因组标签
annoGenomeLabel(
    plot = genes_gm, 
    params = params_a, 
    scale = "Kb", 
    fontsize = 6,
    y = 3.0
)

# 关闭PDF设备
dev.off()

library(strawr)
hicFile <- "~/F1_OSN/LiMCA/HIC/merge_impute/Camkmt/Paternal_specific.impute.hic"

data1 <- straw(
    norm = "NONE",
    fname = hicFile,
    chr1loc = "17(MAT):84755000:85110000",
    chr2loc = "17(MAT):84755000:85110000",
    unit = "BP",
    binsize = 5000
)

sum(data1$counts)




r$> merge_count%>%subset(gene=="Epha7")
      Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype  gene mean_expression meanexp_rank
Epha7    8712     5623 1925  11035   11052    4589   Others Epha7        4.006308        30874

r$> merge_count%>%subset(gene=="Kcnh7")
      Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype  gene mean_expression meanexp_rank
Kcnh7    6368     1237  461   8201    2901    1297   Others Kcnh7        3.005917        30728

r$> merge_count%>%subset(gene=="Ago3")
     Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype gene mean_expression meanexp_rank
Ago3    7832     8330 2607  10111   11790    4537   Others Ago3        2.652439        30668

r$> merge_count%>%subset(gene=="Syn2")
     Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype gene mean_expression meanexp_rank
Syn2    7625     2741  841  10487    6785    2599   Others Syn2        2.793072        30699

r$> merge_count%>%subset(gene=="Dmxl2")
      Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype  gene mean_expression meanexp_rank
Dmxl2    5748    11475 2544   9877   14723    5422   Others Dmxl2        2.460704        30615

r$> merge_count%>%subset(gene=="Clip4")
      Onlyrna Onlypeak Both Orirna Oripeak Oriboth genetype  gene mean_expression meanexp_rank
Clip4    5528     7098 1987   6877   10295    3372   Others Clip4        1.607743        30147







epha7_g1 <- G1_matrix["Epha7",]
epha7_g2 <- G2_matrix["Epha7", ]

# 找出G1为0或G2为0，但不同时为0的列
selected_cols <- colnames(G1_matrix)[(G1_matrix["Epha7",] == 0 | G2_matrix["Epha7", ] == 0) & !(G1_matrix["Epha7",] == 0 & G2_matrix["Epha7", ] == 0)]

G1_matrix["Epha7",selected_cols]

selected_cols <- colnames(G1_matrix)[(G1_matrix["Adcy3",] == 0 | G2_matrix["Adcy3", ] == 0) & !(G1_matrix["Adcy3",] == 0 & G2_matrix["Adcy3", ] == 0)]

a<-names(G2_matrix["Adcy3",G2_matrix["Adcy3",] > 0])

b<-names(G1_matrix["Adcy3",G1_matrix["Adcy3",] > 0])

trans_bed<-read.table(file=str_c(rfdir,"regions/transcripts.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
tss_bed<-read.table(file=str_c(rfdir,"regions/tss.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
promoter_bed<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-1000-1,
                         strand=="-" ~end-1000),end=case_when(strand=="+" ~ end+1000-1,
                         strand=="-" ~end+1000))
chrom<-paste0("chr",c(1:19,"X","Y"))
tss_bed_new<-tss_bed%>%subset(chr %in% chrom)
chr_gene<-id2gene%>%subset(geneid  %in% tss_bed_new$gene)
chrXgene<-unique(trans_bed$gene[! is.na(str_extract(trans_bed$chr,"chrX"))])
chrXgene<-id2gene%>%subset(geneid %in% chrXgene)
chrXgene<-chrXgene$gene