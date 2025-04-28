# plot the track on alleles of different expression pattern group
# first make a big Obj that merge the cell from maternal obj and paternal obj
#####################################################################
mOSN_barcode<-WhichCells(Combined_mOSN, idents = "Mature OSNs")
annotation<-import('/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/gencode.vM23.primary_assembly.annotation.gtf')
annotation<-annotation[,c("type","transcript_id","gene_name","gene_id","gene_type")][annotation$type%in%c("exon","CDS","UTR"),]
annotation$tx_id<-annotation$transcript_id
annotation$gene_biotype<-annotation$gene_type
annotation$type<-as.character(annotation$type)
annotation$type<-replace(annotation$type,annotation$type=="CDS","cds")
annotation$type<-replace(annotation$type,annotation$type=="UTR","utr")
annotation$type<-factor(annotation$type,levels = c("utr","cds","exon"))
annotation<-annotation[,c("type","tx_id","gene_name","gene_id","gene_biotype")]

dirs<-c("/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_allele/joint202201_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_allele/joint202201_G2count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_allele/jointF1OSN_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_allele/jointF1OSN_G2count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_allele/jointF18w_G1count/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_allele/jointF18w_G2count/outs/") 
dirs<-matrix(dirs,ncol=3,nrow=2)
objList <- lapply(1:nrow(dirs),function(i){
  objsub<-lapply(1:ncol(dirs),function(j){
    counts <- Read10X_h5(str_c(dirs[i,j],"filtered_feature_bc_matrix.h5")) 
    fragpath <- str_c(dirs[i,j],"atac_fragments.tsv.gz")
    metadata <- read.csv(str_c(dirs[i,j],"per_barcode_metrics.csv"),row.names=1) %>% dplyr::filter(is_cell==1) 
    obj <- CreateSeuratObject(
      counts = counts$`Gene Expression`,
      assay = "RNA",
      meta.data = metadata,
    )
    obj[["ATAC"]] <- CreateChromatinAssay(
      counts = counts$Peaks,
      sep = c(":", "-"),
      fragments = fragpath,
      annotation = annotation
    )
    obj
  })
  combined.peaks <- UnifyPeaks(object.list = list(objsub[[1]][["ATAC"]], objsub[[2]][["ATAC"]],objsub[[3]][["ATAC"]]), mode = "reduce")
  data1.new_counts <- FeatureMatrix(
    fragments = Fragments(objsub[[1]][["ATAC"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(objsub[[1]][["ATAC"]])
  )
  data2.new_counts <- FeatureMatrix(
    fragments = Fragments(objsub[[2]][["ATAC"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(objsub[[2]][["ATAC"]])
  )
  data3.new_counts <- FeatureMatrix(
    fragments = Fragments(objsub[[3]][["ATAC"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(objsub[[3]][["ATAC"]])
  )
  objsub[[1]][['ATAC']] <- CreateChromatinAssay(counts =data1.new_counts,sep = c(":", "-"),fragments = Fragments(objsub[[1]][["ATAC"]]),annotation = annotation)
  objsub[[2]][['ATAC']] <- CreateChromatinAssay(counts =data2.new_counts,sep = c(":", "-"),fragments = Fragments(objsub[[2]][["ATAC"]]),annotation = annotation)
  objsub[[3]][['ATAC']] <- CreateChromatinAssay(counts =data3.new_counts,sep = c(":", "-"),fragments = Fragments(objsub[[3]][["ATAC"]]),annotation = annotation)
  obj.merge <- merge(x =objsub[[1]],y =list(objsub[[2]],objsub[[3]]),merge.data = TRUE,project = "OSN")
  obj.merge<-subset(obj.merge,cells=mOSN_barcode)
  obj.merge
})
combined.peaks <- UnifyPeaks(object.list = list(objList[[1]][["ATAC"]], objList[[2]][["ATAC"]]), mode = "reduce")
data1.new_counts <- FeatureMatrix(
    fragments = Fragments(objList[[1]][["ATAC"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(objList[[1]][["ATAC"]])
)
data2.new_counts <- FeatureMatrix(
    fragments = Fragments(objList[[2]][["ATAC"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(objList[[2]][["ATAC"]])
)
objList[[1]][['ATAC']] <- CreateChromatinAssay(counts =data1.new_counts,sep = c(":", "-"),fragments = Fragments(objList[[1]][["ATAC"]]),annotation = annotation)
objList[[2]][['ATAC']] <- CreateChromatinAssay(counts =data2.new_counts,sep = c(":", "-"),fragments = Fragments(objList[[2]][["ATAC"]]),annotation = annotation)
# add G1 and G2 to cell ids to distinguish the information from two alleles
Obj<-merge(x =objList[[1]],y =list(objList[[2]]),merge.data = TRUE,add.cell.ids =c("G1","G2"))
save(Obj,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/Obj.RData")

# function of promoter track
DefaultAssay(Obj)<-"ATAC"
heat_track<-function(gene_name){
  gene_id<-id2gene%>%subset(gene==gene_name)
  gene_id<-gene_id[1,1]
  promoter_bed<-max_promoter_data%>%subset(gene==gene_name)%>%select(chr,start,end)
  gene_bed<-bt.merge(subset(trans_bed,gene==gene_id))
  gene_bed<-data.frame(chr=c(gene_bed[1,1]),start=min(gene_bed[,2]),end=max(gene_bed[,3]))
  if(promoter_bed$start<gene_bed$start){
    total_bed<-str_c(gene_bed[1,1],"-",promoter_bed$start,"-",promoter_bed$end+10000)
  }else{
    total_bed<-str_c(gene_bed[1,1],"-",promoter_bed$start-10000,"-",promoter_bed$end)
  }
  
  heat_data<-final_gene_dscore%>%subset(gene==gene_name,select=c(gene,barcode,mergescore,merge_dscore_type))
  rownames(heat_data)<-heat_data$barcode
  # annotation plot
  ap<- AnnotationPlot(object = Obj,region = total_bed)+theme(text=element_text(size=8),
                axis.title.y = element_blank(),  
          axis.text.y = element_blank(),   
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  total_barcode<-c(str_c("G1_",heat_data$barcode),str_c("G2_",heat_data$barcode))
  Obj<-subset(x =Obj, cells=total_barcode)
  max_exp<-max(Obj[["RNA"]]@data[gene_name,])
  max_exp<-ceiling(max_exp)
  # track plot
  trackp <- lapply(c("Paternal_specific","Maternal_specific","Biallelic"),function(i){
    barcode<-heat_data%>%subset(merge_dscore_type==i)
    if(nrow(barcode)>0){
      barcode<-barcode$barcode
      barcode<-c(str_c("G1_",barcode),str_c("G2_",barcode))
      temp_obj<-subset(x =Obj, cells=barcode)
      label<-do.call(rbind,strsplit(WhichCells(temp_obj),"_"))
      label<-label[,1]
      temp_obj@meta.data$type<-str_c(i,"_",label)
      p<-CoveragePlot(
        object = temp_obj,
        region = total_bed,
        group.by="type",
        annotation = FALSE,
        tile = FALSE,
        peaks=FALSE)&scale_fill_manual(values=c("#505f84","#b45a55"))&theme(
          axis.title = element_blank(),  
          axis.text = element_blank(),  
          axis.ticks = element_blank(),
          axis.line = element_blank()
          )
      exp<-ExpressionPlot(temp_obj, features = gene_name,group.by="type", assay = "RNA")+scale_fill_manual(values=c("#505f84","#b45a55"))+
      theme(text=element_text(size=8))+scale_x_continuous(limits = c(0, max_exp))
      list(p,exp)
    }
    else{
      NULL                                                                                                                                                                                                                                                                                                                            
    }
  })
  # track of total reads and total cells
  p<-CoveragePlot(
    object = Combined_mOSN,
    region = total_bed,
    annotation = FALSE,
    tile = FALSE,
    peaks=FALSE)&scale_fill_manual(values=c("#c5cd94"))&theme(
          axis.title = element_blank(),  
          axis.text = element_blank(),  
          axis.line = element_blank(),
          axis.ticks = element_blank()
          )
  
  return(p+plot_spacer()+trackp[[1]][[1]]+trackp[[1]][[2]]+trackp[[2]][[1]]+trackp[[2]][[2]]+trackp[[3]][[1]]+trackp[[3]][[2]]+ap+plot_spacer())
}
# Circular stacked profile 
reads_bar_gene<-function(gene_name){
  data<-final_gene_dscore%>%subset(gene==gene_name)
  data$total<-data$allele_exp+data$allele_frag
  data$G1<-data$G1rna+data$G1peak
  data_temp<-data%>%arrange(total,desc(G1),desc(allele_frag))
  data<-data%>%select(gene,G1peak,G2peak,G1rna,G2rna,barcode)%>%pivot_longer(!c("gene","barcode"),names_to="Group",values_to="count")
  data<-data%>%mutate(Group=case_when(Group=="G1peak"~"ATAC_Paternal",Group=="G2peak"~"ATAC_Maternal",Group=="G1rna"~"RNA_Paternal",Group=="G2rna"~"RNA_Maternal"))
  data$Group<-factor(data$Group,levels=c("ATAC_Paternal","RNA_Paternal","ATAC_Maternal","RNA_Maternal"))
  data$barcode<-factor(data$barcode,levels=data_temp$barcode)
  P<-ggplot(data%>%subset(count>0),aes(x=barcode,y=count,fill=Group))+geom_bar(stat = "identity")+coord_polar()+
    scale_fill_manual(values=c("#526188","#91aebd","#8c0703","#c18686"))+
    labs(x="",y="")+
    theme(legend.title = element_blank(),
        legend.text = element_text(colour = 'black',size =8),
        legend.key.size = unit(10, "pt"),
        legend.position = "top",  # 将图例放置在底部
        legend.box = "horizontal",   # 将图例项并排放置
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.ticks = element_blank(),
        text=element_text(size=10)
        )
  return(P)
}

# Figure 3a and 3b
for(i in c("Lingo2","Hsbp1")){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/high/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}

# Extended Data Fig.3 a,b
# OR ranking
dist_data%>%subset(type=="OR")%>%arrange(abs(p_signal/(m_signal+p_signal)-0.5),desc(count))
# select OR (Olfr1178,Olfr309)
top_OR<-c("Olfr1178","Olfr309")
for(i in top_OR){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/OR/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}

for(i in high_RME$gene){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/high/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}

# Circular stacked profile 
# Extended Data Fig.3 c,d
# chrX ranking
top2_chrX<-dist_data%>%subset(type=="chrX"&count>=100)%>%arrange(abs(p_signal/(m_signal+p_signal)-0.5),desc(count))%>%head(n=2)
for(i in top2_chrX$gene){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/chrX/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}

for(i in c("Ophn1","Diaph2")){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/chrX/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}


# 与NPC gene之间的交集
NPC_RME<-read.xlsx("/data/R02/huangtt39/data/mouse/ESC_NPC_RME.xlsx",sheet = 2,startRow=4)
RME_sameNPC<-RME%>%subset(gene%in%NPC_RME$Gene.Name) #2 genes

# 与Adut DAR DAEE gene之间的交集
DRN_RME<-read.xlsx("/data/R02/huangtt39/data/mouse/DAEE_in_Mouse.xlsx",sheet = 1)
DRN_RME<-DRN_RME%>%subset(DAEEs==TRUE&X_Auto=="Autosomal"&`rab_value_95%CI_UpperLimit`<=0.5&ImprintedFDR05==FALSE) # 328 genes
RME_sameDRN<-RME%>%subset(gene%in%DRN_RME$external_gene_id) # 5 genes

# 与Adut ARN DAEE gene之间的交集
ARN_RME<-read.xlsx("/data/R02/huangtt39/data/mouse/DAEE_in_Mouse.xlsx",sheet = 4)
ARN_RME<-ARN_RME%>%subset(DAEEs==TRUE&X_Auto=="Autosomal"&`rab_value_95%CI_UpperLimit`<=0.5&ImprintedFDR05==FALSE)
RME_sameARN<-RME%>%subset(gene%in%ARN_RME$external_gene_id) #6 genes

# get the gene with similar promoter accessibility signal 
high_RME<-dist_data%>%subset(type=="Others"&count-Biallelic>=25&abs(p_signal/(p_signal+m_signal)-0.5)<=0.1)
write.xlsx(high_RME,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/high_RME.xlsx")
high_RME<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/high_RME.xlsx")
high_RME1<-no_biallelic_ratio%>%subset(gene%in%high_RME$gene)
write.xlsx(high_RME1,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/high_RME1.xlsx")

# high RME rank
# Figure 3c
library(ggprism)
high_RME_ratio<-no_biallelic_ratio%>%subset(gene%in%high_RME$gene)%>%arrange(desc(1-Biallelic), desc(Maternal_specific)) 
high_RME_ratio$gene<-factor(high_RME_ratio$gene,levels=high_RME_ratio$gene)

high_RME_count_p<- ggplot(high_RME_ratio, aes(x = gene, y = count)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_classic() +
  labs(x = "Gene", y = "Total cell") +
  theme(axis.text.x = element_blank())    

high_RME_ratio2p<-high_RME_ratio%>% 
            pivot_longer(cols = c(Paternal_specific, Maternal_specific, Biallelic), 
               names_to = "dscore_type", values_to = "ratio")
high_RME_ratio2p$dscore_type<-factor(high_RME_ratio2p$dscore_type,levels=c("Paternal_specific","Maternal_specific","Biallelic"))

high_RME_ratio_p<-ggplot(high_RME_ratio2p, aes(x = gene, y = ratio*100, fill = dscore_type)) +
  geom_bar(stat = "identity", position = "fill") +  # 100% 堆积柱状图
  scale_y_continuous(labels = scales::percent) +  # Y轴转换为百分比
  scale_fill_manual(values=c("#94b2c2","#c68989","#dfe2d6"))+
  theme_classic() +
  labs(x = "Gene", y = "Proportion", fill = "type")+
  theme(text=element_text(size=10),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
    axis.text.y=element_text(size=10,color = "black"),
    axis.text.x=element_text(size=10,  color = "black",angle = 45, hjust = 1))  
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/high_RME_rank.pdf",high_RME_count_p+high_RME_ratio_p+plot_layout(heights=c(1,1.5)), device = "pdf", width = 24, height = 14, units = "cm", dpi = 300,bg="white")

# Extended Data Fig. 4a,b
Canonical_RME_ratio<-no_biallelic_ratio%>%subset(type%in%c("chrX","OR"))%>%arrange(desc(1-Biallelic), desc(Maternal_specific)) 
Canonical_RME_ratio$gene<-factor(Canonical_RME_ratio$gene,levels=Canonical_RME_ratio$gene)

OR_count_p<- ggplot(Canonical_RME_ratio%>%subset(type=="OR"), aes(x = gene, y = count)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_classic() +
  labs(x = "Gene", y = "Total cell") +
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())   

chrX_count_p<- ggplot(Canonical_RME_ratio%>%subset(type=="chrX"), aes(x = gene, y = count)) +
  geom_bar(stat = "identity", fill = "gray") +
  theme_classic() +
  labs(x = "Gene", y = "Total cell") +
  theme(axis.text.x = element_blank(),axis.ticks.x=element_blank())    

Canonical_RME_ratio2p<-Canonical_RME_ratio%>% 
            pivot_longer(cols = c(Paternal_specific, Maternal_specific, Biallelic), 
               names_to = "dscore_type", values_to = "ratio")
Canonical_RME_ratio2p$dscore_type<-factor(Canonical_RME_ratio2p$dscore_type,levels=c("Paternal_specific","Maternal_specific","Biallelic"))

OR_ratio_p<-ggplot(Canonical_RME_ratio2p%>%subset(type=="OR"), aes(x = gene, y = ratio*100, fill = dscore_type)) +
  geom_bar(stat = "identity", position = "fill") +  # 100% 堆积柱状图
  scale_y_continuous(labels = scales::percent) +  # Y轴转换为百分比
  scale_fill_manual(values=c("#94b2c2","#c68989","#dfe2d6"))+
  theme_classic() +
  labs(x = "Gene", y = "Proportion", fill = "type")+
  theme(text=element_text(size=10),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
    axis.text.y=element_text(size=10,color = "black"),
    axis.text.x=element_blank(),axis.ticks.x=element_blank())  

chrX_ratio_p<-ggplot(Canonical_RME_ratio2p%>%subset(type=="chrX"), aes(x = gene, y = ratio*100, fill = dscore_type)) +
  geom_bar(stat = "identity", position = "fill") +  # 100% 堆积柱状图
  scale_y_continuous(labels = scales::percent) +  # Y轴转换为百分比
  scale_fill_manual(values=c("#94b2c2","#c68989","#dfe2d6"))+
  theme_classic() +
  labs(x = "Gene", y = "Proportion", fill = "type")+
  theme(text=element_text(size=10),
    plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
    axis.text.y=element_text(size=10,color = "black"),
    axis.text.x=element_blank(),axis.ticks.x=element_blank())  
 
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/OR_rank.pdf",OR_count_p+OR_ratio_p+plot_layout(heights=c(0.5,1.5)), device = "pdf", width = 16, height = 10, units = "cm", dpi = 300,bg="white")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/chrX_rank.pdf",chrX_count_p+chrX_ratio_p+plot_layout(heights=c(0.5,1.5)), device = "pdf", width = 21, height = 10, units = "cm", dpi = 300,bg="white")






