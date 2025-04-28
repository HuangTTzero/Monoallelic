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
load("~/ATAC-RNAseq/F1joint_mask/combined.RData")
Idents(combined)<-combined@meta.data$celltype
mOSN_barcode<-WhichCells(combined, idents = "Mature OSNs")
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- "UCSC"


# dataset1
data1_mOSN<-gsub('_1','',mOSN_barcode[grep('1$',mOSN_barcode)])
write.table(data.frame(barcode=data1_mOSN),'/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/OSN/joint202201_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
#10294
# Load the data
dirs <- c("/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_allele/joint202201_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_allele/joint202201_G2count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_mask/joint202201_mask/outs/")
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
save(objList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/"
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
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data1_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/fragcount.RData")
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
save(readList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/readList.RData")


# dataset2
data2_mOSN<-gsub('_2','',mOSN_barcode[grep('2$',mOSN_barcode)])
write.table(data.frame(barcode=data2_mOSN),'/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/OSN/jointF1OSN_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
# 9960
# Load the data
dirs <- c("/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_allele/jointF1OSN_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_allele/jointF1OSN_G2count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_mask/jointF1OSN_mask/outs/")
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
save(objList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/"
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
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data2_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/fragcount.RData")
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
save(readList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/readList.RData")


#dataset3
data3_mOSN<-gsub('_3','',mOSN_barcode[grep('3$',mOSN_barcode)])
write.table(data.frame(barcode=data3_mOSN),'/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/OSN/jointF18w_OSN_barcode.txt',sep = '\t', quote = FALSE, row.names = FALSE,col.names=FALSE)
# 9647 cells
# Load the data
dirs <- c("/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_allele/jointF18w_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_allele/jointF18w_G2count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_mask/jointF18w_mask/outs/")
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
save(objList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/objList.RData")
# fragment annotation of promoter of each gene
fragmentobj <- lapply(1:3,function(i){
  obj<-read_tsv(file=str_c(dirs[i],"atac_fragments.tsv.gz"),comment = "#",col_names =c("chr","start","end","barcode","count"))
  obj$peak<-paste(obj$chr,obj$start,obj$end,sep="-")
  obj
})
rfdir<-"/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/"
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
write.table(fragmentobj[[1]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/pwk_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[2]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/c57_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
write.table(fragmentobj[[3]], gzfile('/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/ori_fragment_anno.txt.gz'),sep = '\t', quote = FALSE, row.names = FALSE)
fragcount<-lapply(1:3,function(i){
  fragmentobj[[i]]%>% group_by(V4,V10,V6)%>%summarise(count=n())%>%ungroup()%>%group_by(V4,V10)%>%summarise(count=n())%>%subset(V4 %in% data3_mOSN&V10 %in%id2gene$V1)
})
rm(fragmentobj)
save(fragcount,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/fragcount.RData")
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
save(readList,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/readList.RData")


# count allelic reads of two allele
# Extended Data Fig. 2a,b
genome<-2725521370
snp<-17686136
genome_exon<-110722912
genome_intron<-1122357321
genome_intergenic<-genome-genome_exon-genome_intron
snp_exon<-0.008*snp
snp_intron<-0.3833*snp
snp_intergenic<-snp-snp_exon-snp_intron
# mutation frequence plot
exon_ratio<-snp_exon/genome_exon #0.001277866
intron_ratio<-snp_intron/genome_intron #0.006040051
intergenic_ratio<-snp_intergenic/genome_intergenic #0.007213384
x <- c('whole genome','exon','intron','intergenic')
y <- c(snp/genome,exon_ratio,intron_ratio,intergenic_ratio)
rate_data <- data.frame(x = x, y = y)
rate_data$x <- factor(rate_data$x, levels=x)
rate_p<-ggplot(data =rate_data, aes(x = x, y = y,fill=x)) + geom_bar(stat = 'identity')+
    geom_text(aes(label=sprintf("%0.4f", round(y, digits = 4))),size=3,vjust=-0.5)+
    scale_fill_manual(values=c("#2f4553","#e8cacb","#b7abd6","#d8bd86"))+labs(x="Rigion",y="Mutation Frequnce",fill="")+
    theme_classic()+
    theme(legend.key.size = unit(8, "pt"),text=element_text(size=8))

# reads distribution in different regions
method<-c("RNA","ATAC")
method1<-c("rna","atac")
sample<-c("joint202201","jointF1OSN","jointF18w")
readtype<-c("G1","G2","ori")
dir<-"/data/R02/huangtt39/ATAC-RNAseq/analysis/"
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
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/reads.pdf",rate_p+p1+p2, device = "pdf", width = 30, height =10 , units = "cm", dpi = 60,bg="white")

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
pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/reads_dis.pdf",width = 5, height =7, )
print(reads_dis_p)
dev.off()


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
load(file="/data/R02/huangtt39/ATAC-RNAseq/analysis/joint202201/readList.RData")

read_List[[1]]<-readList
rm(readList)
load(file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF1OSN/readList.RData")

read_List[[2]]<-readList
rm(readList)
load(file="/data/R02/huangtt39/ATAC-RNAseq/analysis/jointF18w/readList.RData")

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
rfdir<-"/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/"
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
save(merge_count,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/merge_count.RData")

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
save(dscore_data,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/RData/dscore_data.RData")

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
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/gene_count.pdf", merge_count2plot_p, device = "pdf", width = 12, height = 7, units = "cm", dpi = 300)

# get gene with both RNA and ATAC allelic reads
merge_count<-merge_count%>%subset(Onlyrna>0&Onlypeak>0)
# remove gene with low expression 
load(file="~/ATAC-RNAseq/F1joint_mask/Combined_mOSN.RData")
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

ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/cutoff_reads_both_genetype.pdf",cutoff_p(readscutoff_both%>%subset(type=="Others"))+cutoff_p(readscutoff_both%>%subset(type=="chrX"))+cutoff_p(readscutoff_both%>%subset(type=="OR"))&geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =36, height = 6, units = "cm", dpi = 300,bg="white")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/cutoff_reads_gene_type.pdf",cutoff_genetype_p(readscutoff_rna)+cutoff_genetype_p(readscutoff_atac)+cutoff_genetype_p(readscutoff_both)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width = 20, height = 5, units = "cm", dpi = 300,bg="white")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/cutoff_reads_both_total.pdf",cutoff_p(readscutoff_both)+geom_vline(xintercept= c(log10(11)),linewidth=0.25,linetype="dashed"), device = "pdf", width =12, height = 6, units = "cm", dpi = 300,bg="white")

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
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/final_gene.pdf",count_final_p,device = "pdf", width =10, height = 10, units = "cm", dpi = 300,bg="white")


# define the allelic expression pattern of each cell
cellcolor<-c("#597b83","#d8bf87","#d9c3c5")

final_gene_dscore<-dscore_data%>%subset(gene %in% final_gene$gene & allele_exp>0&allele_frag>0&(allele_exp+allele_frag>=4)) 
final_gene_dscore<-final_gene_dscore%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
final_gene_dscore<-final_gene_dscore%>%mutate(rna_dscore_type=case_when(rnascore==-0.5 ~ "Maternal_specific",rnascore>(-0.5)&rnascore<0.5 ~ "Biallelic",rnascore==0.5 ~ "Paternal_specific"))

# the expression pattern of each gene 
# add positive control, imprinting gene
NPC_RME<-read.xlsx("/data/R02/huangtt39/data/mouse/ESC_NPC_RME.xlsx",sheet = 2,startRow=4)
imprint<-read.xlsx("/data/R02/huangtt39/data/mouse/imprint_new.xlsx",sheet = 1)
imprint_Pa<-imprint%>%subset(Expressed.allele=="P")
imprint_Ma<-imprint%>%subset(Expressed.allele=="M")


final_gene_celltype_count<-final_gene_dscore%>%group_by(gene,merge_dscore_type)%>%summarise(cellnumber=n())
                                                                                                                                                                                                                                                <-merge(final_gene_celltype_count,final_gene%>%dplyr::select(gene,count,type),by="gene")
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
ggsave(filename ="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure1/Positve.pdf",Positve_p, device = "pdf", width = 10, height =7, units = "cm", dpi = 300,bg="white")


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
save(sample_100,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/sample_100.RData")

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
# Figure 1c
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


# caculated the  threshold when cumulative probability across monoallelic cells =0.01   
simulation_count<-data.frame(count=10:10000)
get_simulation<-function(x){
  return(qbinom(0.01,x,0.125,lower.tail=FALSE))
}
simulation_count$cutoff<-do.call(rbind,map(simulation_count$count,get_simulation))


# the gene pass the cutoff
# Figure 1d
count_ratio_cutoff_p<-ggplot() +geom_area(data=simulation_count,aes(y=cutoff/count,x=log10(count)),fill="#DEDEDE",alpha=0.7)+
  geom_point(data=final_count_ratio, aes(x = log10(count), y =Maternal_specific+Paternal_specific,color=type),size=1,alpha=0.9) +
  scale_color_manual(values=c(cellcolor,"#2323AF"))+
  scale_y_continuous(limits=c(0,1))+
  geom_line(data=simulation_count,aes(x=log10(count),y=cutoff/count),color="grey")+
  labs(x="# cell of each gene (log10) ",
       y="% Parental specific cell",legend)+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=10)
)

# no_biallelic genes in permutation
get_permutation_fre<-function(df,num,type,color){
  p<-ggplot(data=df%>%group_by(round)%>%summarise(count=n()),aes(x=count))+
      geom_density()+
      geom_vline(xintercept = num,lty="dashed",color=color)+
      scale_x_continuous(limits=c(0,400))+
      annotate(geom = "text",x=num-50,y=0.05,label="Observed\n349",
            vjust=-2,size=3)+
      theme_classic()+theme(text=element_text(size=10))+
      labs(y="Density",
          x=paste0("# ",type," in permutations"))
      return(p)
}

permutation_fre<-get_permutation_fre(merge(sample_count,simulation_count,by="count")%>%subset((Maternal_specific+Paternal_specific)>=cutoff),349," RMEs","#851d34")
temp<-merge(sample_count,simulation_count,by="count")%>%subset((Maternal_specific+Paternal_specific)>=cutoff)
temp<-temp%>%group_by(round)%>%summarise(count=n())
mean(temp$count)
#[1] 22.36
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/get_nobiallelic_gene.pdf",count_ratio_cutoff_p+permutation_fre+Total_ratio_p+plot_layout(ncol=3,widths=c(1.5,1,1)), device = "pdf", width = 24, height = 7, units = "cm", dpi = 300,bg="white")


# Get no_biallelic genes
final_count<-final_gene_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count),values_fill=0)
final_count$type<-as.character(final_count$type)
final_count<- final_count %>%mutate(type = ifelse(gene %in% c("Peg3","Kcnq1ot1"), "imprinting", type))

final_count$type<-factor(final_count$type,levels=c("OR","chrX","Others","imprinting"))
temp_final_gene<-merge(final_count,simulation_count,by="count")

no_biallelic<-temp_final_gene%>%subset((Maternal_specific+Paternal_specific)>=cutoff)#349,包括所有的OR gene和chrX gene
no_biallelic_ratio<-final_count_ratio%>%subset(gene%in%no_biallelic$gene)
write.xlsx(no_biallelic_ratio,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/no_biallelic_gene.xlsx")
no_biallelic_ratio<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/no_biallelic_gene.xlsx")

candicate_gene<-no_biallelic%>%subset(type=="Others"&gene!="Peg3"&gene!="Kcnq1ot1")#281,去除印记基因
candicate_gene_ratio<-final_count_ratio%>%subset(gene%in%candicate_gene$gene)%>%mutate(Parental_specific=Maternal_specific+Paternal_specific)
candicate_gene_ratio<-candicate_gene_ratio%>%mutate(PM_type=case_when(Paternal_specific==0~"Only_Ma",Maternal_specific==0~"Only_Pa",Paternal_specific>0&Maternal_specific>0~"BI_RME"))
table(candicate_gene_ratio$PM_type)
#  BI_RME Only_Ma Only_Pa 
#     257      11      13 

# Figure 1e
# biallelic gene VS no_biallelic gene cellnumber 
temp_final_gene<-final_gene%>%subset(type=="Others")%>%mutate(type=case_when(gene%in%candicate_gene$gene~"P<0.01",!gene%in%candicate_gene$gene~"P>=0.01"))
no_biallelic_count_p<-temp_final_gene%>%ggboxplot(x="type", y="log10(count)",
              palette = "npg",
              xlab = F, #不显示x轴的标签
              bxp.errorbar=T,#显示误差条
              bxp.errorbar.width=0.5, #误差条大小
              size=1, #箱型图边线的粗细
              legend = "right")+scale_color_manual(values=c("#377ebd","#4eae48"))+theme(text = element_text(size = 10))+theme_classic()

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
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/nobiallelic_gene.pdf",no_biallelic_count_p+candicate_gene_p+plot_layout(ncol=2,widths=c(1,1.5)), device = "pdf", width = 14, height = 7, units = "cm", dpi = 300,bg="white")


# final RME
RME<-candicate_gene_ratio%>%subset(PM_type=="BI_RME")


# ranking the RME genes with the promoter signal of C57 and PWK group
# load the fragment file to identify the detail signal of promoter region
# Figure 1f
sample<-c("joint202201","jointF1OSN","jointF18w")
fragcount_dir<-"/data/R02/huangtt39/ATAC-RNAseq/analysis/"
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
save(fragment,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/fragment.RData")

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

no_biallelic_add<-no_biallelic%>%subset(!gene%in%dist_data$gene)

# use the function
max_promoter_data_add<-data.frame()
for(i in no_biallelic_add$gene){
  temp_promoter<-get_max_promoter(i)
  max_promoter_data_add<-rbind(max_promoter_data_add,temp_promoter)
}
dist_data_add<-data.frame()
for(i in no_biallelic_add$gene[6:91]){
  temp_dis<-get_dis(i)
  dist_data_add<-rbind(dist_data_add,temp_dis)
}

save(dist_data_add,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/dist_add.RData")
save(max_promoter_data_add,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/max_promoter_data_add.RData")

dist_data<-rbind(dist_data,dist_data_add)
max_promoter_data<-rbind(max_promoter_data,max_promoter_data_add)

save(dist_data,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/dist.RData")
save(max_promoter_data,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/max_promoter_data.RData")

dist_data<-dist_data%>%mutate(p_type=case_when(cor_p<0.05~"P<0.05",cor_p>0.05|is.na(cor_p)~"P>0.05|P=NA"))
dist_data<-merge(dist_data,no_biallelic)
dist_data$difference_score<-abs(dist_data$p_signal/(dist_data$p_signal+dist_data$m_signal)-0.5)

no_biallelic_ratio<-merge(no_biallelic_ratio,dist_data%>%select(gene,difference_score))
write.xlsx(no_biallelic_ratio,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/no_biallelic_gene.xlsx")


dist_data<-dist_data%>%arrange(desc(abs(p_signal/(p_signal+m_signal)-0.5)),count-Biallelic)
dist_data$gene<-factor(dist_data$gene,levels=dist_data$gene)
dist_p<-ggplot(dist_data%>%subset(gene%in%RME$gene))+
    geom_point(aes(x=gene,y=abs(p_signal/(p_signal+m_signal)-0.5),color=log10(count-Biallelic)))+
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
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/dist.pdf",dist_p, device = "pdf", width = 11, height = 6, units = "cm", dpi = 300,bg="white")
