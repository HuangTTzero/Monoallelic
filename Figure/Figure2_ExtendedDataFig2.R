# get allelic matrix from each sample
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
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
library("bedtoolsr")
library("parallel")
load("~/F1_OSN/F1joint_mask/Combined_mOSN.RData")
Idents(Combined_mOSN)<-Combined_mOSN@meta.data$celltype
mOSN_barcode<-WhichCells(Combined_mOSN, idents = "Mature OSNs")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"


# dataset1
data1_mOSN<-gsub('_1','',mOSN_barcode[grep('1$',mOSN_barcode)])
write.table(data.frame(barcode=data1_mOSN),'/data/R02/huangtt39/F1_OSN/analysis/joint202201/OSN/joint202201_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
#10294
# Load the data
dirs <- c("/data/R02/huangtt39/F1_OSN/mapping/joint202201/mapping_allele/joint202201_G1count/outs/","/data/R02/huangtt39/F1_OSN/mapping/joint202201/mapping_allele/joint202201_G2count/outs/","/data/R02/huangtt39/F1_OSN/mapping/joint202201/mapping_mask/joint202201_mask/outs/")
samples <- c("G1","G2","Ori")
objList <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_feature_bc_matrix.h5")) 
  fragpath <- str_c(dirs[i],"atac_fragments.tsv.gz")
  metadata <- read.csv(str_c(dirs[i],"per_barcode_metrics.csv"),row.names=1) %>% filter(is_cell==1) 
  obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    meta.data = metadata,
    project = samples[i]
  )
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  obj<-subset(x =obj, cells=data1_mOSN)
  obj
})
save(objList,file="/data/R02/huangtt39/F1_OSN/analysis/joint202201/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/F1_OSN/mapping/reference/mm10-mask/"
id2gene<-read.table(file=str_c(rfdir,"genes/id_to_name.txt"))
id2gene<-id2gene[!duplicated(id2gene$V2),]
colnames(id2gene)<-c("geneid","gene")
trans_bed<-read.table(file=str_c(rfdir,"regions/transcripts.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
tss_bed<-read.table(file=str_c(rfdir,"regions/tss.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
promoter_bed<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-1000-1,
                         strand=="-" ~end-1000),end=case_when(strand=="+" ~ end+1000-1,
                         strand=="-" ~end+1000))
 
fragmentobj<-lapply(1:3,function(i){
  obj<-bt.intersect(a = fragmentobj[[i]], b = promoter_bed,wa=TRUE,wb=TRUE)
})
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/joint202201/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/joint202201/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/joint202201/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data1_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/F1_OSN/analysis/joint202201/fragcount.RData")
###构建reads count objlist
library("parallel")
cl <- makeCluster(8)
clusterExport(cl, c("id2gene","data1_mOSN","fragcount","objList"))
readList <- parLapply(cl,c(1:3),function(i){
  objList<-objList[[i]]
  fragment<-fragcount[[i]]
  library(Seurat)
  library(dplyr)
  library(dbplyr)
  library(tidyr)
  obj <- list()
  barcode <- WhichCells(objList)
  fragment_gene<-unique(fragment$V10)
  obj[["gene"]] <- objList$RNA@counts[id2gene$V2,]
  null_rna_matrix<-matrix(0, nrow =length(id2gene$V2) , ncol =length(setdiff(data1_mOSN,barcode)) , byrow = TRUE, dimnames = list(id2gene$V2,setdiff(data1_mOSN,barcode)))
  new_rna_matrix<-cbind(obj[["gene"]],null_rna_matrix)
  obj[["gene"]]<-new_rna_matrix[id2gene$V2,data1_mOSN]
  obj[["peak"]]<-matrix(0, nrow =length(id2gene$V2) , ncol =length(data1_mOSN) , byrow = TRUE, dimnames = list(id2gene$V2, data1_mOSN))
  for(geneid in fragment_gene){
    genename<-id2gene%>%subset(V1==geneid)
    genename<-genename[1,2]
    obj[["peak"]][genename,subset(fragment,V10 ==geneid)$V4]<-subset(fragment,V10 ==geneid)$count
  }
  obj
})
stopCluster(cl)
save(readList,file="/data/R02/huangtt39/F1_OSN/analysis/joint202201/readList.RData")


# dataset2
data2_mOSN<-gsub('_2','',mOSN_barcode[grep('2$',mOSN_barcode)])
write.table(data.frame(barcode=data2_mOSN),'/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/OSN/jointF1OSN_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
# 9960
# Load the data
dirs <- c("/data/R02/huangtt39/F1_OSN/mapping/jointF1OSN/mapping_allele/jointF1OSN_G1count/outs/","/data/R02/huangtt39/F1_OSN/mapping/jointF1OSN/mapping_allele/jointF1OSN_G2count/outs/","/data/R02/huangtt39/F1_OSN/mapping/jointF1OSN/mapping_mask/jointF1OSN_mask/outs/")
samples <- c("G1","G2","Ori")
objList <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_feature_bc_matrix.h5")) 
  fragpath <- str_c(dirs[i],"atac_fragments.tsv.gz")
  metadata <- read.csv(str_c(dirs[i],"per_barcode_metrics.csv"),row.names=1) %>% filter(is_cell==1) 
  obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    meta.data = metadata,
    project = samples[i]
  )
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  obj<-subset(x =obj, cells=data2_mOSN)
  obj
})
save(objList,file="/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/F1_OSN/mapping/reference/mm10-mask/"
id2gene<-read.table(file=str_c(rfdir,"genes/id_to_name.txt"))
id2gene<-id2gene[!duplicated(id2gene$V2),]
trans_bed<-read.table(file=str_c(rfdir,"regions/transcripts.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
tss_bed<-read.table(file=str_c(rfdir,"regions/tss.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
promoter_bed<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-1000-1,
                         strand=="-" ~end-1000),end=case_when(strand=="+" ~ end+1000-1,
                         strand=="-" ~end+1000))
 
fragmentobj<-lapply(1:3,function(i){
  obj<-bt.intersect(a = fragmentobj[[i]], b = promoter_bed,wa=TRUE,wb=TRUE)
})
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data2_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/fragcount.RData")
###构建reads count objlist
library("parallel")
cl <- makeCluster(8)
clusterExport(cl, c("id2gene","data2_mOSN","fragcount","objList"))
readList <- parLapply(cl,c(1:3),function(i){
  objList<-objList[[i]]
  fragment<-fragcount[[i]]
  library(Seurat)
  library(dplyr)
  library(dbplyr)
  library(tidyr)
  obj <- list()
  barcode <- WhichCells(objList)
  fragment_gene<-unique(fragment$V10)
  obj[["gene"]] <- objList$RNA@counts[id2gene$V2,]
  null_rna_matrix<-matrix(0, nrow =length(id2gene$V2) , ncol =length(setdiff(data2_mOSN,barcode)) , byrow = TRUE, dimnames = list(id2gene$V2,setdiff(data2_mOSN,barcode)))
  new_rna_matrix<-cbind(obj[["gene"]],null_rna_matrix)
  obj[["gene"]]<-new_rna_matrix[id2gene$V2,data2_mOSN]
  obj[["peak"]]<-matrix(0, nrow =length(id2gene$V2) , ncol =length(data2_mOSN) , byrow = TRUE, dimnames = list(id2gene$V2, data2_mOSN))
  for(geneid in fragment_gene){
    genename<-id2gene%>%subset(V1==geneid)
    genename<-genename[1,2]
    obj[["peak"]][genename,subset(fragment,V10 ==geneid)$V4]<-subset(fragment,V10 ==geneid)$count
  }
  obj
})
stopCluster(cl)
save(readList,file="/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/readList.RData")


#dataset3
data3_mOSN<-gsub('_3','',mOSN_barcode[grep('3$',mOSN_barcode)])
write.table(data.frame(barcode=data3_mOSN),'/data/R02/huangtt39/F1_OSN/analysis/jointF18w/OSN/jointF18w_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
# 9647 cells
# Load the data
dirs <- c("/data/R02/huangtt39/F1_OSN/mapping/jointF18w/mapping_allele/jointF18w_G1count/outs/","/data/R02/huangtt39/F1_OSN/mapping/jointF18w/mapping_allele/jointF18w_G2count/outs/","/data/R02/huangtt39/F1_OSN/mapping/jointF18w/mapping_mask/jointF18w_mask/outs/")
samples <- c("G1","G2","Ori")
objList <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_feature_bc_matrix.h5")) 
  fragpath <- str_c(dirs[i],"atac_fragments.tsv.gz")
  metadata <- read.csv(str_c(dirs[i],"per_barcode_metrics.csv"),row.names=1) %>% filter(is_cell==1) 
  obj <- CreateSeuratObject(
    counts = counts$`Gene Expression`,
    assay = "RNA",
    meta.data = metadata,
    project = samples[i]
  )
  obj[["ATAC"]] <- CreateChromatinAssay(
    counts = counts$Peaks,
    sep = c(":", "-"),
    fragments = fragpath,
    annotation = annotation
  )
  obj<-subset(x =obj, cells=data3_mOSN)
  obj
})
save(objList,file="/data/R02/huangtt39/F1_OSN/analysis/jointF18w/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/F1_OSN/mapping/reference/mm10-mask/"
id2gene<-read.table(file=str_c(rfdir,"genes/id_to_name.txt"))
id2gene<-id2gene[!duplicated(id2gene$V2),]
trans_bed<-read.table(file=str_c(rfdir,"regions/transcripts.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
tss_bed<-read.table(file=str_c(rfdir,"regions/tss.bed"),header=FALSE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",gene="V4",score="V5",strand="V6"))
promoter_bed<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-1000-1,
                         strand=="-" ~end-1000),end=case_when(strand=="+" ~ end+1000-1,
                         strand=="-" ~end+1000))
 
fragmentobj<-lapply(1:3,function(i){
  obj<-bt.intersect(a = fragmentobj[[i]], b = promoter_bed,wa=TRUE,wb=TRUE)
})
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF18w/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF18w/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/F1_OSN/analysis/jointF18w/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data3_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/F1_OSN/analysis/jointF18w/fragcount.RData")
###构建reads count objlist
library("parallel")
cl <- makeCluster(8)
clusterExport(cl, c("id2gene","data3_mOSN","fragcount","objList"))
readList <- parLapply(cl,c(1:3),function(i){
  objList<-objList[[i]]
  fragment<-fragcount[[i]]
  library(Seurat)
  library(dplyr)
  library(dbplyr)
  library(tidyr)
  obj <- list()
  barcode <- WhichCells(objList)
  fragment_gene<-unique(fragment$V10)
  obj[["gene"]] <- objList$RNA@counts[id2gene$V2,]
  null_rna_matrix<-matrix(0, nrow =length(id2gene$V2) , ncol =length(setdiff(data3_mOSN,barcode)) , byrow = TRUE, dimnames = list(id2gene$V2,setdiff(data3_mOSN,barcode)))
  new_rna_matrix<-cbind(obj[["gene"]],null_rna_matrix)
  obj[["gene"]]<-new_rna_matrix[id2gene$V2,data3_mOSN]
  obj[["peak"]]<-matrix(0, nrow =length(id2gene$V2) , ncol =length(data3_mOSN) , byrow = TRUE, dimnames = list(id2gene$V2, data3_mOSN))
  for(geneid in fragment_gene){
    genename<-id2gene%>%subset(V1==geneid)
    genename<-genename[1,2]
    obj[["peak"]][genename,subset(fragment,V10 ==geneid)$V4]<-subset(fragment,V10 ==geneid)$count
  }
  obj
})
stopCluster(cl)
save(readList,file="/data/R02/huangtt39/F1_OSN/analysis/jointF18w/readList.RData")


# count allelic reads of two allele
# Extended Data Fig. 2a,b
genome<-2725521370
snp<-17686136
genome_exon<-110722912
genome_intron<-1122357321
genome_intergenic<-genome-genome_exon-genome_intron
snp_exon<-0.01626*snp
snp_intron<-0.5731*snp
snp_intergenic<-snp-snp_exon-snp_intron
# mutation frequence plot
exon_ratio<-snp_exon/genome_exon 
intron_ratio<-snp_intron/genome_intron 
intergenic_ratio<-snp_intergenic/genome_intergenic
x <- c('whole genome','exon','intron','intergenic')
y <- c(snp/genome,exon_ratio,intron_ratio,intergenic_ratio)
rate_data <- data.frame(x = x, y = y)
rate_data$x <- factor(rate_data$x, levels=x)
rate_p<-ggplot(data =rate_data, aes(x = x, y = y,fill=x)) + geom_bar(stat = 'identity')+
    geom_text(aes(label=sprintf("%0.4f", round(y, digits = 4))),size=3,vjust=-0.5)+
    scale_fill_manual(values=c("#2f4553","#e8cacb","#b7abd6","#d8bd86"))+labs(x="Rigion",y="SNV density",fill="")+
    theme_classic()+
    theme(legend.key.size = unit(8, "pt"),text=element_text(size=8))
write.xlsx(rate_data, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2a.xlsx")

# reads distribution in different regions
method<-c("RNA","ATAC")
method1<-c("rna","atac")
sample<-c("joint202201","jointF1OSN","jointF18w")
readtype<-c("G1","G2","ori")
dir<-"/data/R02/huangtt39/F1_OSN/analysis/"
reads<-lapply(1:2,function(i){
    final<-lapply(readtype,function(type){
      temp<-lapply(sample,function(smp){
        file<-str_c(dir,smp,"/OSN/",method[i],"/",type,"/",smp,"_OSN_",method1[i],"_",type,".txt")
        data<-read.table(file,skip=4,nrow=4,header=TRUE)
        total<-read.table(file,skip=1,nrow=1,header=FALSE)
        final_data<-data.frame(position=c("Exon","Intron","Intergenic"),values=c(sum(data$Tag_count[1:3]),data$Tag_count[4],total$V3-sum(data$Tag_count)),sample=rep(type,3))
        final_data
      })
      total_exon<-temp[[1]]$values[1]+temp[[2]]$values[1]+temp[[3]]$values[1]
      total_intron<-temp[[1]]$values[2]+temp[[2]]$values[2]+temp[[3]]$values[2]
      total_intergenic<-temp[[1]]$values[3]+temp[[2]]$values[3]+temp[[3]]$values[3]
      final_temp<-data.frame(position=c("Exon","Intron","Intergenic"),values=c(total_exon,total_intron,total_intergenic))
      final_temp$type<-rep(type,3)
      final_temp
    })
    Unassignable_exon<-final[[3]]$values[1]-final[[1]]$values[1]-final[[2]]$values[1]
    Unassignable_intron<-final[[3]]$values[2]-final[[1]]$values[2]-final[[2]]$values[2]
    Unassignable_intergenic<-final[[3]]$values[3]-final[[1]]$values[3]-final[[2]]$values[3]
    Unassignable<-data.frame(position=c("Exon","Intron","Intergenic"),values=c(Unassignable_exon,Unassignable_intron,Unassignable_intergenic))
    Unassignable$type=rep("un",3)
    data<-rbind(final[[1]],final[[2]],Unassignable)
})
for(i in 1:2){
  reads[[i]]$ratio<-reads[[i]]$values/sum(reads[[i]]$values)
  temp_a<-reads[[i]]%>%group_by(position)%>%summarise(ratio=sum(ratio),values=sum(values))
  temp_a$x<-rep("a",3)
  reads[[i]]$x<-rep("b",nrow(reads[[i]]))
  reads[[i]]$position<-paste(reads[[i]]$position,reads[[i]]$type,sep="_")
  reads[[i]]<-rbind(reads[[i]]%>%dplyr::select(position,values,ratio,x),temp_a)
  reads[[i]]$position<-factor(reads[[i]]$position,levels=c("Exon","Intron","Intergenic","Exon_G1","Exon_G2","Exon_un","Intron_G1","Intron_G2","Intron_un","Intergenic_G1","Intergenic_G2","Intergenic_un"))
}
p1<-ggplot(reads[[1]]) +
  # 绘制柱状图
  geom_bar(aes(x, 
               ratio, 
               fill = position),
           stat = 'identity', width = 1.3,color="white") +
  # 设置Y轴刻度
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("#e8cacb","#b7abd6","#d8bd86","#526188","#bb5c56","#c0cacc","#526188","#bb5c56","#c0cacc","#526188","#bb5c56","#c0cacc"))+
  coord_polar(theta = "y") + # 转换坐标轴
  theme_void() # 隐藏图例
p2<-ggplot(reads[[2]]) +
  # 绘制柱状图
  geom_bar(aes(x, 
               ratio, 
               fill = position),
           stat = 'identity', width = 1.3,color="white") +
  # 设置Y轴刻度
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c("#e8cacb","#b7abd6","#d8bd86","#526188","#bb5c56","#c0cacc","#526188","#bb5c56","#c0cacc","#526188","#bb5c56","#c0cacc"))+
  coord_polar(theta = "y") + # 转换坐标轴
  theme_void() # 隐藏图例
reads[[1]]$dataset<-"RNA"
reads[[2]]$dataset<-"ATAC"
reads<-rbind(reads[[1]],reads[[2]])
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/reads.pdf",rate_p+p1+p2, device = "pdf", width = 30, height =10 , units = "cm", dpi = 60,bg="white")
write.xlsx(reads%>%select(position,values,ratio,dataset), "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2b.xlsx")

# Extended Data Fig. 2c
reads<-lapply(1:2,function(i){
    final<-lapply(readtype[1:2],function(type){
      temp<-lapply(sample,function(smp){
        file<-str_c(dir,smp,"/OSN/",method[i],"/",type,"/",smp,"_OSN_",method1[i],"_",type,".txt")
        file_total<-str_c(dir,smp,"/OSN/",method[i],"/ori/",smp,"_OSN_",method1[i],"_ori.txt")
        allele<-read.table(file,skip=0,nrow=1,header=FALSE)
        total<-read.table(file_total,skip=0,nrow=1,header=FALSE)
        data<-data.frame(Parental=type,total=total$V3,allele=allele$V3,type=method[i])
      })
      temp<-do.call(rbind,temp)
      temp<-temp%>%group_by(Parental,type)%>%summarise(total=sum(total),allele=sum(allele))
    })
    final<-do.call(rbind,final)
})
reads_dis<-do.call(rbind,reads)
reads_dis$type<-factor(reads_dis$type,levels=c("RNA","ATAC"))
reads_dis$Parental<-factor(reads_dis$Parental,levels=c("G1","G2"))
reads_dis_p<-ggplot(reads_dis,aes(type,round(allele/total, 4)*100,fill=Parental))+geom_col(position=position_dodge(width=0.9))+
  scale_fill_manual(values=c("#526188","#bb5c56"))+
  theme_classic()+
  scale_y_continuous(expand = c(0,0),limits=c(0,25))+
  geom_text(aes(label = round(allele/total, 4)*100),colour = "black",size = 2,vjust =(-0.2),position = position_dodge(.9))+
  labs(x=NULL,title = NULL,y="% of allelic reads")+
  theme(legend.title = element_blank(),
        legend.position=c(.8, .9),
        legend.text = element_text(colour = 'black',size =8),
        legend.key.size = unit(10, "pt"),
        axis.text = element_text(colour = 'black',size = 8),
        panel.grid = element_blank(),
        text=element_text(size=10)
        )
pdf("/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/reads_dis.pdf",width = 5, height =7, )
print(reads_dis_p)
dev.off()
reads_dis$ratio<-reads_dis$allele/reads_dis$total
write.xlsx(reads_dis, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2c.xlsx")

# load readList of each sample to get total allelic reads matrix
# load library
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
######合并数据#####
read_List<-list()
load(file="/data/R02/huangtt39/F1_OSN/analysis/joint202201/readList.RData")

read_List[[1]]<-readList
rm(readList)
load(file="/data/R02/huangtt39/F1_OSN/analysis/jointF1OSN/readList.RData")

read_List[[2]]<-readList
rm(readList)
load(file="/data/R02/huangtt39/F1_OSN/analysis/jointF18w/readList.RData")

read_List[[3]]<-readList
rm(readList)

for(i in 1:3){
  colnames(read_List[[i]][[1]][["gene"]])<-str_c(colnames(read_List[[i]][[1]][["gene"]]),"_",i) 
  colnames(read_List[[i]][[1]][["peak"]])<-str_c(colnames(read_List[[i]][[1]][["peak"]]),"_",i) 
  colnames(read_List[[i]][[2]][["gene"]])<-str_c(colnames(read_List[[i]][[2]][["gene"]]),"_",i) 
  colnames(read_List[[i]][[2]][["peak"]])<-str_c(colnames(read_List[[i]][[2]][["peak"]]),"_",i) 
  colnames(read_List[[i]][[3]][["gene"]])<-str_c(colnames(read_List[[i]][[3]][["gene"]]),"_",i) 
  colnames(read_List[[i]][[3]][["peak"]])<-str_c(colnames(read_List[[i]][[3]][["peak"]]),"_",i) 
}
readList<-list()
readList[["G1_rna"]]<-cbind(read_List[[1]][[1]][["gene"]],read_List[[2]][[1]][["gene"]],read_List[[3]][[1]][["gene"]])
readList[["G2_rna"]]<-cbind(read_List[[1]][[2]][["gene"]],read_List[[2]][[2]][["gene"]],read_List[[3]][[2]][["gene"]])
readList[["Orirna"]]<-cbind(read_List[[1]][[3]][["gene"]],read_List[[2]][[3]][["gene"]],read_List[[3]][[3]][["gene"]])
readList[["G1_peak"]]<-cbind(read_List[[1]][[1]][["peak"]],read_List[[2]][[1]][["peak"]],read_List[[3]][[1]][["peak"]])
readList[["G2_peak"]]<-cbind(read_List[[1]][[2]][["peak"]],read_List[[2]][[2]][["peak"]],read_List[[3]][[2]][["peak"]])
readList[["Oripeak"]]<-cbind(read_List[[1]][[3]][["peak"]],read_List[[2]][[3]][["peak"]],read_List[[3]][[3]][["peak"]])
rm(read_List)

readList[["rna"]]<-readList[["G1_rna"]]+readList[["G2_rna"]]
readList[["peak"]]<-readList[["G1_peak"]]+readList[["G2_peak"]]
readList[["G1"]]<-readList[["G1_rna"]]+readList[["G1_peak"]]
readList[["G2"]]<-readList[["G2_rna"]]+readList[["G2_peak"]]
readList[["merge"]]<-readList[["G1"]]+readList[["G2"]]

# genes were detected
merge_count<-cbind(as.data.frame(rowSums(readList[["rna"]] >0&readList[["Orirna"]]>0)),as.data.frame(rowSums(readList[["peak"]]>0&readList[["Oripeak"]]>0)),as.data.frame(rowSums(readList[["rna"]] >0&readList[["peak"]]>0)),as.data.frame(rowSums(readList[["Orirna"]]>0)),as.data.frame(rowSums(readList[["Oripeak"]]>0)),as.data.frame(rowSums(readList[["Oripeak"]]>0&readList[["Orirna"]]>0))) %>% rename_with(~ c("Onlyrna","Onlypeak","Both","Orirna","Oripeak","Oriboth"), 1:6) %>% filter(Orirna>0|Oripeak>0)
#31222 genes
# load tss, transcript bedfile to identify chrX gene
rfdir<-"/data/R02/huangtt39/F1_OSN/mapping/reference/mm10-mask/"
id2gene<-read.table(file=str_c(rfdir,"genes/id_to_name.txt"))
id2gene<-id2gene[!duplicated(id2gene$V2),]
colnames(id2gene)<-c("geneid","gene")
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
# load OR names
ORgene<-read.xlsx("/data/R02/huangtt39/data/OSN/12864_2020_6583_MOESM2_ESM.xlsx",sheet = 2)
ORgene <- unique(ORgene$Gene.symbol)
# add gene type
merge_count$genetype<-rep("Others",nrow(merge_count))
for (i in 1:nrow(merge_count)) { 
  gene <- rownames(merge_count)[i]
  if(gene %in% chrXgene ){
    type="chrX"
    merge_count[i,]$genetype=type
  }
  if(gene %in% ORgene ){
    type="OR"
    merge_count[i,]$genetype=type
  }
} 
merge_count$gene<-rownames(merge_count)
# get gene in chromsomes
merge_count<-merge_count%>%subset(gene %in%chr_gene$gene)
merge_count%>%group_by(genetype)%>%summarize(count=n())
# # A tibble: 3 × 2
#   genetype count
#   <chr>    <int>
# 1 OR        1130
# 2 Others   29019
# 3 chrX      1025
save(merge_count,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/merge_count.RData")

# get gene_cell_reads_dataframe
cl <- makeCluster(6)
clusterExport(cl, c("merge_count","readList"))
dscore_data <- parLapply(cl,1:nrow(merge_count),function(i){
  library(Seurat)
  library(dplyr)
  gene <- rownames(merge_count)[i]
  cells<-which(readList[["Orirna"]][gene,]>0|readList[["Oripeak"]][gene,]>0)
  temp_oripeak<-readList[["Oripeak"]][gene,]
  temp_origene<-readList[["Orirna"]][gene,]
  final_oripeak<-temp_oripeak[cells]
  final_oriexp<-temp_origene[cells]
  temp_G1peak<-readList[["G1_peak"]][gene,]
  temp_G2peak<-readList[["G2_peak"]][gene,]
  temp_G1gene<-readList[["G1_rna"]][gene,]
  temp_G2gene<-readList[["G2_rna"]][gene,]
  final_G1peak<-temp_G1peak[cells]
  final_G2peak<-temp_G2peak[cells]
  final_G1gene<-temp_G1gene[cells]
  final_G2gene<-temp_G2gene[cells]
  temp_matrix<-cbind(as.data.frame(rep(gene, length(final_oripeak))),as.data.frame(rep(merge_count[i,]$genetype, length(final_oripeak))),as.data.frame(final_G1peak),as.data.frame(final_G2peak),as.data.frame(final_G1gene),as.data.frame(final_G2gene),as.data.frame(final_oripeak),as.data.frame(final_oriexp)) %>% rename_with(~ c("gene","type","G1peak","G2peak","G1rna","G2rna","Oripeak","Orirna"), 1:8)
  temp_matrix
})
stopCluster(cl)
dscore_data<-do.call(rbind,dscore_data)
dscore_data$barcode <- str_sub(rownames(dscore_data), start = 1,end=20)

dscore_data$rnascore<-dscore_data$G1rna/(dscore_data$G1rna+dscore_data$G2rna)-0.5
dscore_data$atacscore<-dscore_data$G1peak/(dscore_data$G1peak+dscore_data$G2peak)-0.5
dscore_data$mergescore<-(dscore_data$G1peak+dscore_data$G1rna)/(dscore_data$G1rna+dscore_data$G2rna+dscore_data$G1peak+dscore_data$G2peak)-0.5
dscore_data$allele_exp<-dscore_data$G1rna+dscore_data$G2rna
dscore_data$allele_frag<-dscore_data$G1peak+dscore_data$G2peak
save(dscore_data,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/dscore_data.RData")

# gene number with allelic reads
# Extended Data Fig. 2d
Orirna<-merge_count%>%subset(Orirna>0)%>%nrow()
Oripeak<-merge_count%>%subset(Oripeak>0)%>%nrow()
Oriboth<-merge_count%>%subset(Orirna>0&Oripeak>0)%>%nrow()
allelerna<-merge_count%>%subset(Onlyrna>0)%>%nrow()
allelepeak<-merge_count%>%subset(Onlypeak>0)%>%nrow()
alleleboth<-merge_count%>%subset(Onlyrna>0&Onlypeak>0)%>%nrow()
merge_count2plot<-data.frame(type=rep(c("Expression","Promoter","Combination"),2),allele_type=c("total","total","total","allele","allele","allele"),count=c(Orirna,Oripeak,Oriboth,allelerna,allelepeak,alleleboth))
#          type allele_type count
# 1  Expression      total 26277
# 2    Promoter      total 31069
# 3 Combination      total 26172
# 4  Expression      allele 22232
# 5    Promoter      allele 28139
# 6 Combination      allele 21599
write.xlsx(merge_count2plot, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2d.xlsx")
merge_count2plot_p<-ggplot(merge_count2plot, aes(x = type, y = count, fill = allele_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("lightgray", "darkcyan")) +  # 设置颜色
  labs(x = "", y = "# genes") +
  theme_minimal() +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/gene_count.pdf", merge_count2plot_p, device = "pdf", width = 12, height = 7, units = "cm", dpi = 300)

# get gene with both RNA and ATAC allelic reads
merge_count<-merge_count%>%subset(Onlyrna>0&Onlypeak>0)
# remove gene with low expression 
load(file="~/F1_OSN/F1joint_mask/Combined_mOSN.RData")
library(ggExtra)
mean_exp <-AverageExpression(Combined_mOSN,assay="RNA",layer="data",features=merge_count$gene)
mean_exp<-mean_exp$RNA
merge_count$mean_expression<-mean_exp[merge_count$gene,]
merge_count$meanexp_rank<-rank(merge_count$mean_expression)
gene_high<-merge_count%>%subset(log10(mean_expression)>=(-2)&Orirna>=100) # cutoff
gene_high$genetype<-factor(gene_high$genetype,levels=c("OR","chrX","Others"))
#     OR   chrX Others 
#    287    386  10986 


# Mountain Plot of the gene number with different thresholds 
# Extended Data Fig. 2e
# the fuction to build data
library(ggridges)
color<-c("#9b9b9a","#d19bd1","#1c2d58")
cutoff_rna<-function(dscore_data){
  final<-data.frame()
  ori<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(Orirna>0))%>%mutate(cutofftype="total UMIs (>=1)")
  cutoff_1<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>0))%>%mutate(cutofftype="allele UMIs (>=1)")
  cutoff_2<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>1))%>%mutate(cutofftype="allele UMIs (>=2)")
  final<-rbind(final,ori,cutoff_1,cutoff_2)
  final$cutofftype<-factor(final$cutofftype,levels=c("total UMIs (>=1)","allele UMIs (>=1)","allele UMIs (>=2)"))
  return(final)
} 
cutoff_atac<-function(dscore_data){
  final<-data.frame()
  ori<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(Oripeak>0))%>%mutate(cutofftype="total fragments (>=1)")
  cutoff_1<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_frag>0))%>%mutate(cutofftype="allele fragments(>=1)")
  cutoff_2<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_frag>1))%>%mutate(cutofftype="allele fragments(>=2)")
  final<-rbind(final,ori,cutoff_1,cutoff_2)
  final$cutofftype<-factor(final$cutofftype,levels=c("total fragments (>=1)","allele fragments(>=1)","allele fragments(>=2)"))
  return(final)
} 
cutoff_both<-function(dscore_data){
  final<-data.frame()
  ori<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(Orirna>0&Oripeak>0))%>%mutate(cutofftype="RNA total >=1 & ATAC total >=1")
  cutoff_1<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>0&allele_frag>0))%>%mutate(cutofftype="RNA allele >=1 & ATAC allele >=1")
  cutoff_2<-dscore_data%>%group_by(gene,type)%>%summarise(count=sum(allele_exp>0&allele_frag>0&(allele_exp+allele_frag>=4)))%>%mutate(cutofftype="RNA allele + ATAC allele >=4")
  final<-rbind(final,ori,cutoff_1,cutoff_2)
  final$cutofftype<-factor(final$cutofftype,levels=c("RNA total >=1 & ATAC total >=1","RNA allele >=1 & ATAC allele >=1","RNA allele + ATAC allele >=4"))
  return(final)
} 
readscutoff_rna<-cutoff_rna(dscore_data%>%subset(gene%in%gene_high$gene))
readscutoff_rna$type<-factor(readscutoff_rna$type,levels=c("Others","chrX","OR"))
readscutoff_atac<-cutoff_atac(dscore_data%>%subset(gene%in%gene_high$gene))
readscutoff_atac$type<-factor(readscutoff_atac$type,levels=c("Others","chrX","OR"))
readscutoff_both<-cutoff_both(dscore_data%>%subset(gene%in%gene_high$gene))
readscutoff_both$type<-factor(readscutoff_both$type,levels=c("Others","chrX","OR"))

cutoff_p<-function(data){
  p<-ggplot(data, aes(x = log10(count+1), fill = cutofftype),alpha=0.6) +geom_histogram(alpha=0.6,position="identity",bins = 20) +
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text=element_text(size=6), axis.text = element_text(colour = 'black',size = 6),)+
  scale_fill_manual(values=color,name="Cell cutoff")+ylab("gene number")+xlab("log10(cell numbers+1)")+
  scale_x_continuous(limit=c(-0.15,4.3),breaks = seq(0,4,1))
  return(p)
}

cutoff_genetype_p<-function(data){
  p<-ggplot(data, aes(x = log10(count+1), fill = cutofftype,y=type,color=cutofftype)) +geom_density_ridges(alpha=0.6)+theme_classic()+theme(plot.title = element_text(hjust = 0.5),legend.position = "none", text=element_text(size=6), axis.text = element_text(colour = 'black',size = 6))+scale_fill_manual(values=color)+scale_color_manual(values=color)+ylab("")+xlab("log10(cell numbers+1)")
  return(p)
}

ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/cutoff_reads_both_genetype.pdf",cutoff_p(readscutoff_both%>%subset(type=="Others"))+cutoff_p(readscutoff_both%>%subset(type=="chrX"))+cutoff_p(readscutoff_both%>%subset(type=="OR"))&geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =36, height = 6, units = "cm", dpi = 300,bg="white")
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/cutoff_reads_gene_type.pdf",cutoff_genetype_p(readscutoff_rna)+cutoff_genetype_p(readscutoff_atac)+cutoff_genetype_p(readscutoff_both)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width = 20, height = 5, units = "cm", dpi = 300,bg="white")
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/cutoff_reads_both_total.pdf",cutoff_p(readscutoff_both)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =12, height = 6, units = "cm", dpi = 300,bg="white")
horizontal_readscutoff_both <- readscutoff_both %>%
  select(gene, count, cutofftype) %>%
  pivot_wider(names_from = cutofftype, values_from = count)
write.xlsx(horizontal_readscutoff_both, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2e.xlsx")

# count the changes of ratio
readscutoff_rna_wide<-readscutoff_rna%>%pivot_wider(names_from="cutofftype",values_from="count",id_cols=c(type,gene),values_fill=0)
mean(readscutoff_rna_wide$`allele UMIs (>=1)`/readscutoff_rna_wide$`total UMIs (>=1)`)
# 0.4698551
mean(readscutoff_rna_wide$`allele UMIs (>=2)`/readscutoff_rna_wide$`total UMIs (>=1)`)
# 0.0459776

readscutoff_atac_wide<-readscutoff_atac%>%pivot_wider(names_from="cutofftype",values_from="count",id_cols=c(type,gene),values_fill=0)
mean(readscutoff_atac_wide$`allele fragments(>=1)`/readscutoff_atac_wide$`total fragments (>=1)`)
# 0.4010039
mean(readscutoff_atac_wide$`allele fragments(>=2)`/readscutoff_atac_wide$`total fragments (>=1)`)
# 0.04909761

readscutoff_both_wide<-readscutoff_both%>%pivot_wider(names_from="cutofftype",values_from="count",id_cols=c(type,gene),values_fill=0)
readscutoff_both_wide<-readscutoff_both_wide%>%subset(`RNA total >=1 & ATAC total >=1`>0)
mean(readscutoff_both_wide$`RNA allele >=1 & ATAC allele >=1`/readscutoff_both_wide$`RNA total >=1 & ATAC total >=1`)
# 0.2260857
mean(readscutoff_both_wide$`RNA allele + ATAC allele >=4`/readscutoff_both_wide$`RNA total >=1 & ATAC total >=1`)
# 0.0193597

# gene with sufficient allelic information
final_gene<-readscutoff_both%>%subset(cutofftype=="RNA allele + ATAC allele >=4"&count>=10)#2657
# Others   chrX     OR 
#   2591     37     29 

count_final_p<-ggplot(final_gene%>%group_by(type)%>%summarise(count=n()),aes(type,count,fill=type))+geom_col(position=position_dodge(width=0.4))+scale_fill_manual(values=cellcolor)+
  scale_y_break(breaks=c(50,2000),ticklabels=seq(2000,2800,200),scales=0.3)+
  theme_classic()+
  geom_text(aes(label = count), vjust =-0.5, size = 2)+
  labs(x=NULL,title = NULL,y="gene number",size=4)+
  theme(legend.title = element_blank(),
        legend.text = element_text(colour = 'black',size = 8),
        legend.key.size = unit(8, "pt"),
        axis.text = element_text(colour = 'black',size = 6),
        panel.grid = element_blank(),
        text = element_text(size = 8)
        )
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/final_gene.pdf",count_final_p,device = "pdf", width =10, height = 10, units = "cm", dpi = 300,bg="white")



####Smart-seq2 analysis
# ExtendedDataFig2f-k
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
  density = c(Cast_genome_ratio, Cast_exon_ratio, DBA_genome_ratio, DBA_exon_ratio)
)

# 确保因子顺序
rate_data$strain <- factor(rate_data$strain, levels = c("Cast", "DBA"))
rate_data$region <- factor(rate_data$region, levels = c("whole genome", "exon"))
write.xlsx(rate_data, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2f.xlsx")

# 绘图 - 按strain分组，x轴为strain
rate_p <- ggplot(data = rate_data, 
                 aes(x = strain, y = density, fill = region)) + 
  geom_bar(stat = 'identity', position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = sprintf("%0.4f", round(density, digits = 4))),
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
LiMCA_reads <- read.table("/data/R02/huangtt39/F1_OSN/mapping/LiMCA/reads_distribution.txt", 
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

ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/LiMCA_reads.pdf",rate_p+reads_vln_plot, device = "pdf", width = 30, height =10 , units = "cm", dpi = 60,bg="white")
write.xlsx(LiMCA_reads_long%>%select(Sample,Total_Fragments,strain,Group,Percentage), "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2g.xlsx")

## load G1 and G2 matrixs of Smart-seq2 dataset
G1_counts <- read.table(
  "/data/R02/huangtt39/F1_OSN/mapping/LiMCA/featurecounts/total_G1.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
G2_counts <- read.table(
  "/data/R02/huangtt39/F1_OSN/mapping/LiMCA/featurecounts/total_G2.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)
ori_counts <- read.table(
  "/data/R02/huangtt39/F1_OSN/mapping/LiMCA/featurecounts/total_dedup.counts_matrix.tsv",
  header = TRUE, sep = "\t", row.names = 1, check.names = FALSE
)

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

  as.data.frame(rowSums(LiMCA_readList[["rna"]] > 0)),

  as.data.frame({
    rna_means <- rowMeans(
      replace(LiMCA_readList[["rna"]], LiMCA_readList[["rna"]] <= 0, NA), 
      na.rm = TRUE
    )

    rna_means[is.nan(rna_means)] <- NA
    rna_means
  }),

  as.data.frame(rowSums(LiMCA_readList[["Orirna"]] > 0))
) %>% 

  rename_with(~ c("rna_count", "rna_mean", "Orirna_count"), 1:3) %>% 

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
save(LiMCA_merge_count,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/LiMCA_merge_count.RData")

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
save(LiMCA_dscore_data,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/LiMCA_dscore_data.RData")
load(file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/LiMCA_dscore_data.RData")

# gene number with allelic reads
LiMCA_Orirna<-LiMCA_merge_count%>%subset(Orirna_count>0)%>%nrow()
LiMCA_allelerna<-LiMCA_merge_count%>%subset(rna_count>0)%>%nrow()
LiMCA_merge_count2plot<-data.frame(allele_type=c("total","allele"),count=c(LiMCA_Orirna,LiMCA_allelerna))
LiMCA_merge_count2plot$allele_type<-factor(LiMCA_merge_count2plot$allele_type,levels=c("total","allele"))
write.xlsx(LiMCA_merge_count2plot, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2h.xlsx")
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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/LiMCA_gene_count.pdf", LiMCA_merge_count2plot_p, device = "pdf", width = 8, height = 7, units = "cm", dpi = 300)


LiMCA_allele_count_meanRNA_p<-ggplot() +
  geom_point(data=LiMCA_merge_count%>%subset(rna_count>0), aes(x = rna_count, y =log10(rna_mean)),color="#567074",size=1,alpha=0.9) +
  labs(x="Number of cells with allelic reads",
       y="Mean allelic reads of cells (log10)")+
  theme_classic()+

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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/LiMCA_allele_count_meanRNA.pdf",LiMCA_allele_count_p+plot_spacer()+LiMCA_allele_count_meanRNA_p+LiMCA_allele_meanRNA_p+plot_layout(nrow=2,ncol=2,height=c(1,2,1,2),width=c(2,1,2,1)), device = "pdf", width = 10, height = 9, units = "cm", dpi = 300,bg="white")
write.xlsx(LiMCA_merge_count, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2i.xlsx")

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
horizontal_LiMCA_readscutoff <- LiMCA_readscutoff %>%
  select(gene, count, cutofftype) %>%
  pivot_wider(names_from = cutofftype, values_from = count)
write.xlsx(horizontal_LiMCA_readscutoff, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2j.xlsx")

cutoff_p<-function(data){
  p<-ggplot(data, aes(x = log10(count+1), fill = cutofftype),alpha=0.6) +geom_histogram(alpha=0.6,position="identity",bins = 20) +
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text=element_text(size=6), axis.text = element_text(colour = 'black',size = 6),)+
  scale_fill_manual(values=color,name="Cell cutoff")+ylab("gene number")+xlab("log10(cell numbers+1)")
  return(p)
}
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/LiMCA_cutoff_reads_total.pdf",cutoff_p(LiMCA_readscutoff)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =12, height = 6, units = "cm", dpi = 300,bg="white")

# for OR genes , don't use the cutoff of cell number
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

LiMCA_to_Positive<-rbind(LiMCA_10_count_ratio_OR,LiMCA_final_count_ratio%>%subset(type!="OR"))
write.xlsx(LiMCA_to_Positive, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS2/FigureS2k.xlsx")


LiMCA_to_Positive_counts <- LiMCA_to_Positive %>%
  subset(type %in% c("OR", "chrX")) %>%
  group_by(Maternal_specific, Paternal_specific, type) %>%
  summarise(count = n(), .groups = "drop")

LiMCA_Positive_p <- ggplot() +
  geom_point(data = LiMCA_to_Positive_counts, 
             aes(x = Maternal_specific, 
                 y = Paternal_specific, 
                 color = type, 
                 size = count), 
             alpha = 0.9) +

  scale_color_manual(values = c(cellcolor[1:2], "#2323AF")) +

  scale_size_continuous(range = c(1, 5), 
                        breaks = c(1, 5, 10, 20), 
                        name = "Number of points") +

  scale_x_continuous(limits = c(0, 1)) +
  scale_y_continuous(limits = c(0, 1)) +

  labs(x = "% Maternal specific cell",
       y = "% Paternal specific cell") +
  
  theme_classic() +

  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right",
        text = element_text(size = 10),
        legend.title = element_text(size = 8),
        legend.text = element_text(colour = 'black', size = 8),
        legend.key.size = unit(10, "pt"))

# 保存图形
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/LiMCA/Positve.pdf",
       LiMCA_Positive_p, 
       device = "pdf", 
       width = 11,  
       height = 7, 
       units = "cm", 
       dpi = 300,
       bg = "white")

# define the allelic expression pattern of each cell
cellcolor<-c("#597b83","#d8bf87","#d9c3c5")

final_gene_dscore<-dscore_data%>%subset(gene %in% final_gene$gene & allele_exp>0&allele_frag>0&(allele_exp+allele_frag>=4)) 
final_gene_dscore<-final_gene_dscore%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
final_gene_dscore<-final_gene_dscore%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))

# the expression pattern of each gene 
# add positive control, imprinting gene
imprint<-read.xlsx("/data/R02/huangtt39/data/mouse/imprint_new.xlsx",sheet = 1)
imprint_Pa<-imprint%>%subset(Expressed.allele=="P")
imprint_Ma<-imprint%>%subset(Expressed.allele=="M")


final_gene_celltype_count<-final_gene_dscore%>%group_by(gene,merge_dscore_type)%>%summarise(cellnumber=n())
final_gene_celltype_count<-merge(final_gene_celltype_count,final_gene%>%dplyr::select(gene,count,type),by="gene")
final_count_ratio<-final_gene_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(type,gene,count),values_fill=0)

final_count_ratio$type<-as.character(final_count_ratio$type)
final_count_ratio<- final_count_ratio %>%mutate(type = ifelse(gene %in% c("Peg3","Kcnq1ot1"), "imprinting", type))
final_count_ratio$type<-factor(final_count_ratio$type,levels=c("OR","chrX","Others","imprinting"))

# Figure 2b
Positve_p<-ggplot() +
  geom_point(data=final_count_ratio%>%subset(type%in%c("chrX","OR")), aes(x = Maternal_specific, y =Paternal_specific ,color=type), size=1,alpha=0.9) +
  scale_color_manual(values=cellcolor)+scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
  # geom_hline(yintercept= c(0.15),linewidth=0.25,linetype="dashed")+geom_vline(xintercept= c(0.15),linewidth=0.25,linetype="dashed")+
  # 坐标轴
  labs(x="% Maternal specific cell",
       y="% Paternal specific cell")+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=10),
        legend.title = element_blank(),
        legend.text = element_text(colour = 'black',size = 8),
        legend.key.size = unit(10, "pt"),
) 
ggsave(filename ="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure1/Positve.pdf",Positve_p, device = "pdf", width = 10, height =7, units = "cm", dpi = 300,bg="white")


# permutation 100 times
test_cell_sample<-function(temp_dscore){
  cell_sample<-function(x){
    G1<-rbinom(100,size=x,p=0.5)
    G2<-x-G1
    return(cbind(G1,G2,1:100))
  }
  sample_rna<-map(temp_dscore$allele_exp,cell_sample)
  sample_atac<-map(temp_dscore$allele_frag,cell_sample)
  sample_rna<-do.call(rbind,sample_rna)%>%as.data.frame()
  colnames(sample_rna)<-c("G1rna","G2rna","round")
  sample_atac<-do.call(rbind,sample_atac)%>%as.data.frame()
  colnames(sample_atac)<-c("G1peak","G2peak","round")
  final<-cbind(sample_atac%>%dplyr::select(G1peak,G2peak),sample_rna)
  return(final)
}


get_sample<-function(genename){
  temp_dscore<-final_gene_dscore%>%subset(gene==genename)
  sample_100<-data.frame()
  sample_100<-test_cell_sample(temp_dscore)
  sample_100$gene<-rep(genename,nrow(sample_100))
  return(sample_100)
}

cl <- makeCluster(16)
clusterExport(cl, c("final_gene_dscore","get_sample","test_cell_sample"))
sample_100<-parLapply(cl,final_gene$gene,function(genename){
  library(dplyr)
  library(dbplyr)
  library(tidyr)
  library(purrr)
  get_sample(genename)
})
stopCluster(cl)
sample_100<-do.call(rbind,sample_100)
save(sample_100,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/sample_100.RData")

sample_100$mergescore<-(sample_100$G1rna+sample_100$G1peak)/(sample_100$G1rna+sample_100$G1peak+sample_100$G2rna+sample_100$G2peak)-0.5
sample_100$rnascore<-(sample_100$G1rna)/(sample_100$G1rna+sample_100$G2rna)-0.5

sample_100<-sample_100%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
sample_100<-sample_100%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))
sample_celltype_count<-sample_100%>%group_by(gene,merge_dscore_type,round)%>%summarise(cellnumber=n())
sample_celltype_count<-merge(sample_celltype_count,final_gene%>%dplyr::select(gene,type,count),by="gene")
sample_count_ratio<-sample_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(type,gene,count,round),values_fill=0)
sample_count<-sample_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count,round),values_fill=0)


# Compare the cell expression pattern proportions in RNA, RNA&ATAC and permutation RNA, permutation RNA&ATAC
# Only chose the reads==4，to control the reads number noise
# Figure 2c
rna<-final_gene_dscore%>%subset(G1rna+G2rna==4)%>%group_by(rna_dscore_type)%>%summarise(count=n())
rna$type<-rep("RNA",3)
rna$dscore_type<-rna$rna_dscore_type
rna$ratio<-rna$count/sum(rna$count)

total<-final_gene_dscore%>%subset(G1rna+G2rna+G1peak+G2peak==4)%>%group_by(merge_dscore_type)%>%summarise(count=n())
total$type<-rep("Total",3)
total$dscore_type<-total$merge_dscore_type
total$ratio<-total$count/sum(total$count)

permutation<-sample_100%>%subset(G1rna+G2rna+G1peak+G2peak==4)%>%group_by(merge_dscore_type)%>%summarise(count=n())
permutation$type<-rep("Permutation",3)
permutation$dscore_type<-permutation$merge_dscore_type
permutation$ratio<-permutation$count/sum(permutation$count)


final_data<-rbind(rna%>%select(type,dscore_type,ratio),permutation%>%select(type,dscore_type,ratio),total%>%select(type,dscore_type,ratio))
final_data$type<-factor(final_data$type,levels=c("Permutation","Total","RNA"))
final_data$dscore_type<-factor(final_data$dscore_type,levels=c("Biallelic","Paternal_specific","Maternal_specific"))

Total_ratio_p<- ggplot(final_data,aes(x=type,y=ratio,fill=dscore_type)) +geom_bar(position = "fill",stat= "identity")+theme_classic()+scale_fill_manual(values=c("#dfe2d6","#94b2c2","#c68989"))
write.xlsx(final_data,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2c.xlsx")


##### caculate monoallelic expression of final_gene in different datasets
# Figure2d
# for 10x RNA dataset cells with more than 2 RNA UMIs were selected
final_gene_dscore_RNA2<-dscore_data%>%subset(gene%in%final_gene$gene & allele_exp>=2)
# permutation 100 times for new cells
# permutation 100 times
test_cell_sample<-function(temp_dscore){
  cell_sample<-function(x){
    G1<-rbinom(100,size=x,p=0.5)
    G2<-x-G1
    return(cbind(G1,G2,1:100))
  }
  sample_rna<-map(temp_dscore$allele_exp,cell_sample)
  sample_atac<-map(temp_dscore$allele_frag,cell_sample)
  sample_rna<-do.call(rbind,sample_rna)%>%as.data.frame()
  colnames(sample_rna)<-c("G1rna","G2rna","round")
  sample_atac<-do.call(rbind,sample_atac)%>%as.data.frame()
  colnames(sample_atac)<-c("G1peak","G2peak","round")
  final<-cbind(sample_atac%>%dplyr::select(G1peak,G2peak),sample_rna)
  return(final)
}


get_sample<-function(genename){
  temp_dscore<-final_gene_dscore_RNA2%>%subset(gene==genename)
  if(nrow(temp_dscore)==0){
    return(data.frame())
  }
  sample_100<-data.frame()
  sample_100<-test_cell_sample(temp_dscore)
  sample_100$gene<-rep(genename,nrow(sample_100))
  return(sample_100)
}

cl <- makeCluster(16)
clusterExport(cl, c("final_gene_dscore_RNA2","get_sample","test_cell_sample"))
sample_100_RNA2<-parLapply(cl,final_gene$gene,function(genename){
  library(dplyr)
  library(purrr)
  get_sample(genename)
})
stopCluster(cl)
sample_100_RNA2<-do.call(rbind,sample_100_RNA2)
save(sample_100_RNA2,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/sample_100_RNA2.RData")
sample_100_RNA2$mergescore<-(sample_100_RNA2$G1rna+sample_100_RNA2$G1peak)/(sample_100_RNA2$G1rna+sample_100_RNA2$G1peak+sample_100_RNA2$G2rna+sample_100_RNA2$G2peak)-0.5
sample_100_RNA2$rnascore<-(sample_100_RNA2$G1rna)/(sample_100_RNA2$G1rna+sample_100_RNA2$G2rna)-0.5

sample_100_RNA2<-sample_100_RNA2%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
sample_100_RNA2<-sample_100_RNA2%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))

# Due to changes in the cell screening conditions, genes from 10 cells are still selected.
final_gene_RNA2<-final_gene_dscore_RNA2%>%subset(allele_frag>0)%>%group_by(gene)%>%summarise(count=n())%>%subset(count>=10)
# now select common genes in 10x genomic and smart-seq2
compare2common_gene<- intersect(final_gene_RNA2$gene, LiMCA_final_gene$gene) # 2199
sample_100_RNA2_common<-sample_100_RNA2%>%subset(gene%in%compare2common_gene)
sample_100_RNA2_common$allele_frag<-sample_100_RNA2_common$G1peak+sample_100_RNA2_common$G2peak
final_RNA2_common_dscore<-final_gene_dscore_RNA2%>%subset(gene%in%compare2common_gene)
final_RNA2_common_dscore<-final_RNA2_common_dscore%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
final_RNA2_common_dscore<-final_RNA2_common_dscore%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))

# calculate_ratios function to calculate cell ratios with different expression type
calculate_ratios <- function(data, type_col, count_info, filter_expr = NULL) {

  if (!is.null(filter_expr)) {
    data <- data %>% filter(eval(parse(text = filter_expr)))
  }
  
  result <- data %>%
    group_by(gene, !!sym(type_col)) %>%
    summarise(cellnumber = n(), .groups = "drop") %>%
    merge(count_info, by = "gene") %>%
    mutate(ratio = cellnumber / count) %>%
    pivot_wider(
      names_from = !!sym(type_col),
      values_from = ratio,
      id_cols = c(gene, type, count),
      values_fill = 0
    )
  return(result)
}

process_permutations<- function(perm_data, type_col, count_info, filter_expr = NULL) {
  if (!is.null(filter_expr)) {
    perm_data <- perm_data %>% filter(eval(parse(text = filter_expr)))
  }

  ratios <- perm_data %>%
    group_by(round, gene, !!sym(type_col)) %>%
    summarise(n = n(), .groups = "drop") %>%
    left_join(count_info, by = "gene") %>%
    mutate(ratio = n / count) %>%
    select(-n, -count)
  
  result <- ratios %>%
    pivot_wider(
      names_from = !!sym(type_col),
      values_from = ratio,
      values_fill = 0
    ) %>%

    {
      df <- .
      cols_needed <- c("Maternal_specific", "Biallelic", "Paternal_specific")
      for (col in cols_needed) {
        if (!col %in% names(df)) df[[col]] <- 0
      }
      df
    } %>%
    group_by(gene) %>%
    summarise(
      across(c("Maternal_specific", "Biallelic", "Paternal_specific"),
             ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    left_join(count_info, by = "gene")
  
  return(result)
}

#  calculate the difference from permutation and 10x multi-omic
calculate_difference<- function(original_ratio, perm_ratio) {
 
  common_genes <- intersect(original_ratio$gene, perm_ratio$gene)
  
  original_sub <- original_ratio %>% filter(gene %in% common_genes)
  perm_sub <- perm_ratio %>% filter(gene %in% common_genes)
  
  # < 0 refer to 0
  diff_result <- original_sub %>%
    select(gene, type, count, 
           monoallelic_original = monoallelic,
           Biallelic_original = Biallelic) %>%
    left_join(
      perm_sub %>% select(gene, 
                         monoallelic_perm = monoallelic,
                         Biallelic_perm = Biallelic),
      by = "gene"
    ) %>%
    mutate(
      monoallelic = pmax(monoallelic_original - monoallelic_perm, 0),
      Biallelic = 1-pmax(monoallelic_original - monoallelic_perm, 0),
    ) %>%
    select(gene, type, count, monoallelic, Biallelic)
  
  return(diff_result)
}
count_info_rna <- final_RNA2_common_dscore %>%
  group_by(gene) %>%
  summarise(
    count = n(),
    type = unique(type)
  )

# for 10x multi-omic data , allele_frag > 0
count_info_merge <- final_RNA2_common_dscore %>%
  filter(allele_frag > 0) %>%
  group_by(gene) %>%
  summarise(
    count = n(),
    type = unique(type)
  )

merge_ratio <- calculate_ratios(
  final_RNA2_common_dscore, 
  "merge_dscore_type", 
  count_info_merge,
  filter_expr = "allele_frag > 0"
)

# RNA ratio for total cells 
rna_ratio <- calculate_ratios(
  final_RNA2_common_dscore, 
  "rna_dscore_type", 
  count_info_rna
)

# for permutation, allele_frag > 0
perm_merge <- process_permutations(
  sample_100_RNA2_common, 
  "merge_dscore_type", 
  count_info_merge,
  filter_expr = "allele_frag > 0"
)

# for RNA permutation,total cells
perm_rna <- process_permutations(
  sample_100_RNA2_common, 
  "rna_dscore_type", 
  count_info_rna
)

merge_ratio$monoallelic<-merge_ratio$Paternal_specific+merge_ratio$Maternal_specific
rna_ratio$monoallelic<-rna_ratio$Paternal_specific+rna_ratio$Maternal_specific

perm_merge$monoallelic<-perm_merge$Paternal_specific+perm_merge$Maternal_specific
perm_rna$monoallelic<-perm_rna$Paternal_specific+perm_rna$Maternal_specific

merge_diff <- calculate_difference(merge_ratio, perm_merge) %>%
  mutate(dataset = "10x_merge_diff")

rna_diff <- calculate_difference(rna_ratio, perm_rna) %>%
  mutate(dataset = "10x_rna_diff")

LiMCA_count_ratio_common<-LiMCA_final_count_ratio%>%subset(gene%in%compare2common_gene)
LiMCA_count_ratio_common$monoallelic<-LiMCA_count_ratio_common$Paternal_specific+LiMCA_count_ratio_common$Maternal_specific
library(Seurat)
load(file="~/F1_OSN/F1joint_mask/Combined_mOSN.RData")
mean_exp <-AverageExpression(Combined_mOSN,assay="RNA",layer="data",features=compare2common_gene)
mean_exp<-mean_exp$RNA

prepare_for_plot <- function(ratio_data, name) {
  ratio_data %>%
    mutate(
      dataset = name,
      avg_expression = mean_exp[ratio_data$gene,],
      monoallelic_ratio = monoallelic
    ) %>%
    select(gene, avg_expression, monoallelic_ratio, dataset, type)
}

plot_data <- bind_rows(
  prepare_for_plot(merge_ratio, "10x_merge_original"),
  prepare_for_plot(rna_ratio, "10x_rna_original"),
  prepare_for_plot(perm_merge, "10x_merge_perm_mean"),
  prepare_for_plot(perm_rna, "10x_rna_perm_mean"),
  prepare_for_plot(merge_diff, "10x_merge_diff"),
  prepare_for_plot(rna_diff, "10x_rna_diff"),
  prepare_for_plot(LiMCA_count_ratio_common, "smart_seq2_rna")
)

library(patchwork)

# select gene type == "Others"
plot_data_others <- plot_data %>%subset(!gene%in%c("Peg3","Kcnq1ot1"))%>%
  filter(type == "Others") %>%
  filter(dataset %in% c("smart_seq2_rna", "10x_rna_original", "10x_rna_diff", 
                        "10x_merge_original", "10x_merge_diff")) %>%
  # rename
  mutate(
    dataset = factor(dataset,
      levels = c("smart_seq2_rna", "10x_rna_original", "10x_rna_diff", 
                 "10x_merge_original", "10x_merge_diff"),
      labels = c("Smart-seq2 RNA", "10X RNA Original", "10X RNA Diff", 
                 "10X Merge Original", "10X Merge Diff")
    )
  )
plot_data_others_toxlsx<-plot_data_others%>%subset(type!="10X RNA Diff") %>%
  pivot_wider(names_from = dataset, values_from = "monoallelic_ratio")
write.xlsx(plot_data_others_toxlsx,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2d.xlsx")

cat("数据分布检查:\n")
cat("总行数:", nrow(plot_data_others), "\n")
cat("基因数量:", n_distinct(plot_data_others$gene), "\n")
# 数据分布检查:
# 总行数: 10860 
# 基因数量: 2172 

# 2. 简洁的绘图函数
plot_scatter_lm <- function(data, ds_name, y_col) {
  
  # 提取该数据集的数据
  ds_data <- data %>% filter(dataset == ds_name)
  
  # 计算线性回归
  if (nrow(ds_data) > 1) {
    lm_fit <- lm(ds_data[[y_col]] ~ ds_data$avg_expression)
    r2 <- summary(lm_fit)$r.squared
    cor_val <- cor(ds_data$avg_expression, ds_data[[y_col]], use = "complete.obs")
    p_val <- summary(lm_fit)$coefficients[2, 4]
  }
  
  # 绘图
  p <- ggplot(ds_data, aes(x = log2(avg_expression), y = .data[[y_col]])) +
    geom_point(alpha = 0.6, size = 1.5, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", 
                fill = "pink", alpha = 0.3) +
    labs(
      title = paste(ds_name, "-", y_col),
      x = "Expression (log2)",
      y = y_col,
      subtitle = if (nrow(ds_data) > 1) 
        sprintf("R² = %.3f, cor = %.3f, p = %.3e", r2, cor_val, p_val)
    ) +
    scale_y_continuous(limits=c(0,1))+
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "darkred")
    )
  
  return(p)
}

plot_scatter_lm <- function(data, ds_name, y_col) {
  
  # 加载必要包
  if (!requireNamespace("hexbin", quietly = TRUE)) {
    install.packages("hexbin")
  }
  library(hexbin)
  
  # 提取该数据集的数据
  ds_data <- data %>% filter(dataset == ds_name)
  
  # 准备数据
  ds_data$x <- ds_data$avg_expression
  ds_data$y <- ds_data[[y_col]]
  
  # 计算线性回归
  if (nrow(ds_data) > 1) {
    lm_fit <- lm(y ~ x, data = ds_data)
    r2 <- summary(lm_fit)$r.squared
    cor_val <- cor(ds_data$x, ds_data$y, use = "complete.obs")
    p_val <- summary(lm_fit)$coefficients[2, 4]
  }
  
  # 创建ggplot版本
  p <- ggplot(ds_data, aes(x = log2(x), y = y)) +
    
    # 使用stat_binhex创建平滑散点图
    stat_binhex(
      bins = 50,  # 控制六边形数量
      aes(fill = after_stat(count)),  # 根据点数填充颜色
      alpha = 0.8
    ) +
    
    # 使用从grey到#d9c3c5的颜色渐变，设置最大值为10
    scale_fill_gradientn(
      colours = c("#5c8b9b","#ccd5dc","#e0c7c3","#8d4550"),
      name = "Density",
      limits = c(0, 10),  # 设置密度范围为0-10
      breaks = c(0, 2, 4, 6, 8, 10),  # 设置刻度
      labels = c("0", "2", "4", "6", "8", "10"),  # 标签，最后一个显示为10+
      guide = guide_colorbar(
        barwidth = 0.5,
        barheight = 8,
        title.position = "top",
        title.hjust = 0.5
      ),
      oob = scales::squish  # 超出范围的值压缩到边界
    ) +
    
    # 添加线性回归线
    geom_smooth(
      method = "lm",
      se = TRUE,
      color = "red",
      fill = rgb(1, 0.8, 0.8, 0.4),
      alpha = 0.3,
      linewidth = 1
    ) +
    
    labs(
      title = paste(ds_name, "-", y_col),
      x = "Expression (log2)",
      y = y_col
    ) +
    
    scale_y_continuous(limits = c(0, 1)) +
    
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, color = "darkred"),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8)
    )
  
  # 添加统计信息作为副标题
  if (nrow(ds_data) > 1) {
    p <- p + labs(subtitle = sprintf("R² = %.3f, cor = %.3f, p = %.3e", 
                                     r2, cor_val, p_val))
  }
  
  return(p)
}

# 3. 创建所有图形并组织布局
# 获取数据集列表
datasets <- c("Smart-seq2 RNA", "10X RNA Original", "10X RNA Diff", 
                 "10X Merge Original", "10X Merge Diff")

# # 创建Paternal图列表
# paternal_plots <- list()
# for (ds in datasets) {
#   paternal_plots[[ds]] <- plot_scatter_lm(plot_data_others, ds, "paternal_ratio")
# }

# # 创建Maternal图列表  
# maternal_plots <- list()
# for (ds in datasets) {
#   maternal_plots[[ds]] <- plot_scatter_lm(plot_data_others, ds, "maternal_ratio")
# }
monoallelic_plots <- list()
for (ds in datasets) {
  monoallelic_plots[[ds]] <- plot_scatter_lm(plot_data_others, ds, "monoallelic_ratio")
}
library(gridExtra)
# 4. 创建PDF文件
pdf("/data/R02/huangtt39/F1_OSN/analysis/figure/LiMCA/scatter_plots_mono.pdf", width = 23, height = 4)

grid.arrange(
  grobs = monoallelic_plots,
  nrow = 1,
  top = "Monoallelic Ratios of different datasets"
)
dev.off()

cat("PDF文件已生成: scatter_plots_combined.pdf\n")


# Get no_biallelic genes
final_count<-final_gene_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count),values_fill=0)
final_count$type<-as.character(final_count$type)
final_count<- final_count %>%mutate(type = ifelse(gene %in% c("Peg3","Kcnq1ot1"), "imprinting", type))

final_count$type<-factor(final_count$type,levels=c("OR","chrX","Others","imprinting"))

# 计算p_value
final_count <- as.data.frame(final_count)
final_count$monoallelic_count <- final_count$Maternal_specific + final_count$Paternal_specific
final_count$p_value <- pbinom(final_count$monoallelic_count - 1, 
                             size = final_count$count, 
                             prob = 0.125, 
                             lower.tail = FALSE)
final_count$q_value <- p.adjust(final_count$p_value, method = "BH")
final_count_ratio<-merge(final_count_ratio,final_count%>%select(gene,p_value,q_value))
save(final_count_ratio,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/final_count_ratio.RData")
write.xlsx(final_count_ratio,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2be.xlsx")

# 查看结果
head(final_count)
# no_bialleic genes
data_pvalue<-final_count%>%subset(p_value<0.01)#209,包括所有的OR gene和chrX gene
no_biallelic<-final_count%>%subset(q_value<0.05)#155,包括所有的OR gene和chrX gene
no_biallelic_ratio<-final_count_ratio%>%subset(gene%in%no_biallelic$gene)

# Figure 2d
final_count_ratio <- final_count_ratio %>%
  mutate(color_group = ifelse(q_value < 0.05, as.character(type), "Not significant"))
final_count_ratio$color_group<-factor(final_count_ratio$color_group,levels=c("OR","chrX","Others","imprinting","Not significant"))
count_ratio_cutoff_p <- ggplot() +
  geom_point(data = final_count_ratio, 
             aes(x = log10(count), 
                 y = Maternal_specific+Paternal_specific,  # 修正：这里应该是比例
                 color = color_group), 
             size = 1, alpha = 0.9) +
  scale_color_manual(
    values = c(
      cellcolor, "#2323AF", # 保留原有的type颜色
      "#e7e7e7"           # 为不显著基因设置灰色
    )
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "# cell of each gene (log10)",
       y = "% Parental specific cell",
       color = "Gene Type") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right",
        text = element_text(size = 10))


# 对permuation结果进行p_value 和q_value 计算，分组矫正
# 为sample_count计算p值和q值（按round分组校正）
sample_count <- sample_count %>%
  group_by(round) %>%
  mutate(
    monoallelic_count = Maternal_specific + Paternal_specific,
    p_value = pbinom(monoallelic_count - 1, size = count, prob = 0.125, lower.tail = FALSE),
    q_value = p.adjust(p_value, method = "BH")
  ) %>%
  ungroup()

sample_count%>%subset(q_value<0.05) #数据为0

ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/get_nobiallelic_gene.pdf",count_ratio_cutoff_p+Total_ratio_p+plot_layout(ncol=2,widths=c(1.5,1)), device = "pdf", width = 20, height = 7, units = "cm", dpi = 300,bg="white")

candicate_gene<-no_biallelic%>%subset(type=="Others"&gene!="Peg3"&gene!="Kcnq1ot1")#87,去除印记基因
candicate_gene_ratio<-final_count_ratio%>%subset(gene%in%candicate_gene$gene)%>%mutate(Parental_specific=Maternal_specific+Paternal_specific)
candicate_gene_ratio<-candicate_gene_ratio%>%mutate(PM_type=case_when(Paternal_specific==0~"Only_Ma",Maternal_specific==0~"Only_Pa",Paternal_specific>0&Maternal_specific>0~"BI_RME"))
table(candicate_gene_ratio$PM_type)
#  BI_RME Only_Ma Only_Pa 
#     72      6      9 

# Figure 2e
# M group & P group of no_biallelic genes
candicate_gene_p<-ggplot() +
  geom_point(data=candicate_gene_ratio, aes(x = Maternal_specific, y =Paternal_specific,color=PM_type,fill=PM_type),size=1,alpha=0.9) +
  scale_fill_manual(values=c("#828e79","#996277","#364f5e"))+
  scale_color_manual(values=c("#828e79","#996277","#364f5e"))+
  scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
  labs(x="% M group",
       y="% P group")+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="None",text=element_text(size=10)
)  
candicate_gene_p<-ggMarginal(candicate_gene_p,groupColour=TRUE,groupFill=TRUE,type = "density",alpha = 0.4)
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/nobiallelic_gene.pdf",candicate_gene_p, device = "pdf", width = 8, height = 7, units = "cm", dpi = 300,bg="white")


# final RME
RME<-candicate_gene_ratio%>%subset(PM_type=="BI_RME")

# Figure 2g cell with consistent alleic promoter accessibility and expression
gene_consistency <- final_gene_dscore %>%subset(!gene%in%c("Peg3","Kcnq1ot1")&type!="chrX")%>%
  group_by(gene, type) %>%
  summarise(
    total_cells = n(),
    consistent_cells = sum(rna_dscore_type == merge_dscore_type, na.rm = TRUE),
    consistency_ratio = consistent_cells / total_cells
  ) %>%
  ungroup()

gene_consistency_modified <- gene_consistency %>%
  mutate(type = ifelse(gene %in% RME$gene, "non-OR_RME", type))
write.xlsx(gene_consistency_modified,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2g.xlsx")

# Wilcoxon
stat_results <- rstatix::wilcox_test(
  gene_consistency_modified,
  consistency_ratio ~ type,
  comparisons = list(c("non-OR_RME", "OR"),
                    c("Others","OR"),
                    c("Others","non-OR_RME")),
  exact = FALSE,
  conf.level = 0.95,
  detailed = TRUE
)

effect_sizes <- gene_consistency_modified %>%
  rstatix::wilcox_effsize(consistency_ratio ~ type,
                 comparisons = list(c("non-OR_RME", "OR"),
                                   c("Others","OR"),
                                   c( "Others","non-OR_RME")))
final_results <- stat_results %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

# # A tibble: 3 × 10
#   group1     group2    n1    n2 statistic        p conf.low conf.high effsize magnitude
#   <chr>      <chr>  <int> <int>     <dbl>    <dbl>    <dbl>     <dbl>   <dbl> <ord>    
# 1 non-OR_RME OR        72    29       1.5 2.85e-15   -0.373    -0.305   0.786 large    
# 2 OR         Others    29  2517   72923   2.15e-20    0.440     0.578   0.183 small    
# 3 non-OR_RME Others    72  2517  129969   3.11e-10    0.112     0.205   0.124 small    

consistency_vlnplot <- ggviolin(gene_consistency_modified, 
                                x = "type", 
                                y = "consistency_ratio", 
                                ylim = c(0, 1.2),
                                alpha = 1,
                                width = 0.5,
                                legend = "none",
                                xlab = "", 
                                ylab = "Consistency ratio",
                                font.tickslab = c(15, "plain", "black"),
                                add = "boxplot", 
                                add.params = list(fill = "white", width = 0.1, linetype = 1)) +
  stat_compare_means(method = "wilcox.test", 
                      label = "p.format", 
                     comparisons = list(c("non-OR_RME", "OR"),c("OR", "Others"),c("non-OR_RME", "Others")), 
                     label.y = c(1.0,1.1,1.05)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position = "right", 
        legend.title = element_blank(),
        text = element_text(size = 8))

# 3. 如果需要保存图片
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/gene_consistency_violin.pdf", 
       consistency_vlnplot, 
       device = "pdf", 
       width = 10, 
       height = 6, 
       units = "cm", 
       dpi = 300)

# ranking the RME genes with the promoter signal of C57 and PWK group
# load the fragment file to identify the detail signal of promoter region
# Figure 2h
sample<-c("joint202201","jointF1OSN","jointF18w")
fragcount_dir<-"/data/R02/huangtt39/F1_OSN/analysis/"
fragcount <- lapply(1:3,function(i){
  frag<-list()
  filepath1<-str_c(fragcount_dir,sample[i],"/pwk_fragment_anno.txt")
  filepath2<-str_c(fragcount_dir,sample[i],"/c57_fragment_anno.txt")
  frag[[1]]<-read.table(file=filepath1,header=TRUE)
  frag[[2]]<-read.table(file=filepath2,header=TRUE)
  for(j in 1:2){
    frag[[j]]<-frag[[j]]%>%select(V1,V2,V3,V4,V10)
    names(frag[[j]])<-c("chr","start","end","barcode","gene")
    frag[[j]]$barcode<-paste(frag[[j]]$barcode,sep="_",i)
  }
  frag
})
fragment<-list()
fragment[[1]]<-rbind(fragcount[[1]][[1]],fragcount[[2]][[1]],fragcount[[3]][[1]])
fragment[[2]]<-rbind(fragcount[[1]][[2]],fragcount[[2]][[2]],fragcount[[3]][[2]])
rm(fragcount)
save(fragment,file="/data/R02/huangtt39/F1_OSN/analysis/merge/fragment.RData")

# the function to get promoter region with the most fragments when the gene has different transcription 
get_max_promoter<-function(gene_name){
  temp<-list()
  cell<-list()
  gene_id<-id2gene%>%subset(gene==gene_name)
  gene_id<-gene_id[1,1]
  cell[[1]]<-final_gene_dscore%>%subset(gene==gene_name&merge_dscore_type=="Paternal_specific")
  cell[[2]]<-final_gene_dscore%>%subset(gene==gene_name&merge_dscore_type=="Maternal_specific")
  for(i in 1:2){
    temp[[i]]<-fragment[[i]]%>%subset(gene==gene_id&barcode%in%cell[[i]]$barcode)
    temp[[i]]<-bt.sort(unique(temp[[i]]))
  }
  promoter_bed<-subset(promoter_bed,gene==gene_id)
  promoter_bed$name<-str_c(promoter_bed$chr,"-",promoter_bed$start,"-",promoter_bed$end)
  
  fragment_count<-bt.intersect(a = promoter_bed, b =rbind(temp[[1]],temp[[2]]))%>%group_by(V7)%>%summarise(count=n())
  max_count_index <- which.max(fragment_count$count)
  max_count_promoter <- fragment_count[max_count_index,1]
  max_promoter<-promoter_bed%>%subset(name==max_count_promoter$V7)
  return(data.frame(gene=c(gene_name),chr=c(max_promoter$chr),start=c(max_promoter$start),end=c(max_promoter$end)))
}
# the function to get the signal on max promoter region 
get_dis<-function(gene_name){
  temp<-list()
  cell<-list()
  gene_id<-id2gene%>%subset(gene==gene_name)
  gene_id<-gene_id[1,1]
  cell[[1]]<-final_gene_dscore%>%subset(gene==gene_name&merge_dscore_type=="Paternal_specific")
  cell[[1]]$sample<-gsub("^.{19}","",cell[[1]]$barcode)
  cell[[2]]<-final_gene_dscore%>%subset(gene==gene_name&merge_dscore_type=="Maternal_specific")
  cell[[2]]$sample<-gsub("^.{19}","",cell[[2]]$barcode)
  frag_total<-c(fragment[[1]]%>%subset(barcode%in%cell[[1]]$barcode)%>%nrow(),fragment[[2]]%>%subset(barcode%in%cell[[2]]$barcode)%>%nrow())
  for(i in 1:2){
    temp[[i]]<-fragment[[i]]%>%subset(gene==gene_id&barcode%in%cell[[i]]$barcode)
    temp[[i]]<-bt.sort(unique(temp[[i]]))
  }
  promoter_bed<-subset(promoter_bed,gene==gene_id)
  promoter_bed$name<-str_c(promoter_bed$chr,"-",promoter_bed$start,"-",promoter_bed$end)
  
  fragment_count<-bt.intersect(a = promoter_bed, b =rbind(temp[[1]],temp[[2]]))%>%group_by(V7)%>%summarise(count=n())
  # Narrow the promoter region to the most concentrated area.
  max_count_index <- which.max(fragment_count$count)
  max_count_promoter <- fragment_count[max_count_index,1]
  max_promoter<-promoter_bed%>%subset(name==max_count_promoter$V7)
  # divide the promoter into bins, to calculate the similarity of the signals. 
  promoter_bin<-bt.makewindows(b=max_promoter,w=50)
  final<-lapply(1:2,function(i){
    final_temp<-bt.intersect(a = promoter_bin, b =temp[[i]],c=TRUE,bed=TRUE)$V4/frag_total[i]*1000000
  })
  if(frag_total[1]==0){
    signal_max<-which.max(final[[2]])
    final[[1]]<-rep(0,40)
  }else if(frag_total[2]==0){
    signal_max<-which.max(final[[1]])
    final[[2]]<-rep(0,40)
  }else{
    signal_max<-which.max(final[[1]])
  }

  if(signal_max<=5){
    signal<-c(final[[1]][1:(signal_max+5)],final[[2]][1:(signal_max+5)])
  }else if(signal_max>=35){
    signal<-c(final[[1]][(signal_max-5):40],final[[2]][(signal_max-5):40])
  }else{
    signal<-c(final[[1]][(signal_max-5):(signal_max+5)],final[[2]][(signal_max-5):(signal_max+5)])
  }
  signal<-matrix(signal, nrow = 2, byrow = TRUE)
  rownames(signal)<-c("Pa","Ma")
  dist_df<-dist(signal,method = "euclidean")
  cor.df<-cor.test(signal[1,],signal[2,])
  dist_df<-as.matrix(dist_df)
  return(data.frame(gene=c(gene_name),dist=c(dist_df[1,2]),cor_p=c(cor.df$p.value),cor_r=c(cor.df$estimate),p_signal=c(max(signal[1,])),m_signal=c(max(signal[2,]))))
}

# use the function
max_promoter_data<-data.frame()
for(i in no_biallelic$gene){
  temp_promoter<-get_max_promoter(i)
  max_promoter_data<-rbind(max_promoter_data,temp_promoter)
}
dist_data<-data.frame()
for(i in no_biallelic$gene){
  temp_dis<-get_dis(i)
  dist_data<-rbind(dist_data,temp_dis)
}

save(dist_data,file="/data/R02/huangtt39/F1_OSN/new_analysis/dist.RData")

dist_data<-merge(dist_data,no_biallelic)
dist_data$difference_score<-abs(dist_data$p_signal/(dist_data$p_signal+dist_data$m_signal)-0.5)

no_biallelic_ratio<-merge(no_biallelic_ratio,dist_data%>%select(gene,p_value,q_value,difference_score))
RME<-no_biallelic_ratio%>%subset(gene%in%RME$gene)%>%select(-type)
write.xlsx(no_biallelic_ratio,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2f.xlsx")
write.xlsx(no_biallelic,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/no_biallelic.xlsx")
write.xlsx(RME,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2h.xlsx")


dist_data<-dist_data%>%arrange(desc(abs(p_signal/(p_signal+m_signal)-0.5)),monoallelic_count)
dist_data$gene<-factor(dist_data$gene,levels=dist_data$gene)
dist_p<-ggplot(dist_data%>%subset(gene%in%RME$gene))+
    geom_point(aes(x=gene,y=abs(p_signal/(p_signal+m_signal)-0.5),color=log10(monoallelic_count)))+
    theme_classic()+
    #scale_y_break(c(-1000, 1000))+theme_prism(palette = "pearl",
    #           base_size = 8,
    #           base_line_size = 0.2,axis_text_angle = 45)+
    # scale_size_continuous(range = c(1,3),breaks=seq(-0.5,0,0.25),labels=abs(seq(-0.5,0,0.25)))+
   scale_color_gradientn(colors =c("#f2f7e8","#D1Dbbd","#3E606F","#193441"))+
   labs(y="Difference score bewteen C57 and PWK group",x="Gene ranking by difference score")+
    theme(text=element_text(size=10),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
      axis.text.y=element_text(size=10,color = "black"),
        axis.ticks.x = element_blank(), # 去掉 x 轴的刻度
        axis.text.x = element_blank()   # 去掉 x 轴的文字
        )    
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/dist.pdf",dist_p, device = "pdf", width = 11, height = 3, units = "cm", dpi = 300,bg="white")
