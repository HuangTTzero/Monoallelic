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
final_gene_celltype_count<-merge(final_gene_celltype_count,final_gene%>%dplyr::select(gene,count,type),by="gene")
final_count_ratio<-final_gene_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(type,gene,count),values_fill=0)

final_count_ratio$type<-factor(final_count_ratio$type,levels=c("OR","chrX","Others"))

# Figure S2A
Positve_p<-ggplot() +
  geom_point(data=final_count_ratio, aes(x = Maternal_specific, y =Paternal_specific ,color=type), size=1,alpha=0.9) +
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
ggsave(filename ="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/Positve.pdf",Positve_p, device = "pdf", width = 10, height =7, units = "cm", dpi = 300,bg="white")


# caculated the  threshold when cumulative probability across monoallelic cells =0.01   
simulation_count<-data.frame(count=10:10000)
get_simulation<-function(x){
  return(qbinom(0.01,x,0.125,lower.tail=FALSE))
}
simulation_count$cutoff<-do.call(rbind,map(simulation_count$count,get_simulation))


# the gene pass the cutoff
#Figure 2B
count_ratio_cutoff_p<-ggplot() +geom_area(data=simulation_count,aes(y=cutoff/count,x=log10(count)),fill="#DEDEDE",alpha=0.7)+
  geom_point(data=final_count_ratio, aes(x = log10(count), y =Maternal_specific+Paternal_specific,color=type),size=1,alpha=0.9) +
  scale_color_manual(values=cellcolor)+
  scale_y_continuous(limits=c(0,1))+
  geom_line(data=simulation_count,aes(x=log10(count),y=cutoff/count),color="grey")+
  labs(x="# cell of each gene (log10) ",
       y="% Parental specific cell",legend)+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right",text=element_text(size=10)
)


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
sample_100<-sample_100%>%mutate(merge_dscore_type=case_when(mergescore==-0.5 ~ "Maternal_specific",mergescore>(-0.5)&mergescore<0.5 ~ "Biallelic",mergescore==0.5 ~ "Paternal_specific"))
sample_celltype_count<-sample_100%>%group_by(gene,merge_dscore_type,round)%>%summarise(cellnumber=n())
sample_celltype_count<-merge(sample_celltype_count,final_gene%>%dplyr::select(gene,type,count),by="gene")
sample_count_ratio<-sample_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="ratio",id_cols=c(type,gene,count,round),values_fill=0)
sample_count<-sample_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count,round),values_fill=0)

# no_biallelic genes in permutation
# Figure S2B
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


# Compare the cell expression pattern proportions in RNA, RNA&ATAC and permutation 
# Figure S2C
temp_other_dscore<-final_gene_dscore%>%subset(type=="Others"&gene!="Peg3"&gene!="Kcnq1ot1")
rna<-temp_other_dscore%>%group_by(rna_dscore_type)%>%summarise(count=n())
rna$type<-rep("RNA",3)
rna$dscore_type<-rna$rna_dscore_type
rna$ratio<-rna$count/sum(rna$count)


total<-temp_other_dscore%>%group_by(merge_dscore_type)%>%summarise(count=n())
total$type<-rep("Total",3)
total$dscore_type<-total$merge_dscore_type
total$ratio<-total$count/sum(total$count)

permutation<-sample_100%>%subset(gene%in%temp_other_dscore$gene)%>%group_by(merge_dscore_type)%>%summarise(count=n())
permutation$type<-rep("Permutation",3)
permutation$dscore_type<-permutation$merge_dscore_type
permutation$ratio<-permutation$count/sum(permutation$count)

final_data<-rbind(permutation%>%select(type,dscore_type,ratio),rna%>%select(type,dscore_type,ratio),total%>%select(type,dscore_type,ratio))
final_data$dscore_type<-factor(final_data$dscore_type,levels=c("Paternal_specific","Biallelic","Maternal_specific"))

Total_ratio_p<- ggplot(final_data,aes(x=type,y=ratio,fill=dscore_type)) +geom_bar(position = "fill",stat= "identity")+theme_classic()+scale_fill_manual(values=c("#94b2c2","#dfe2d6","#c68989"))

ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/get_nobiallelic_gene.pdf",count_ratio_cutoff_p+permutation_fre+Total_ratio_p+plot_layout(ncol=3,widths=c(1.5,1,1)), device = "pdf", width = 24, height = 7, units = "cm", dpi = 300,bg="white")


# Get no_biallelic genes
final_count<-final_gene_celltype_count%>%mutate(ratio=cellnumber/count)%>%pivot_wider(names_from="merge_dscore_type",values_from="cellnumber",id_cols=c(type,gene,count),values_fill=0)
final_count$type<-factor(final_count$type,levels=c("OR","chrX","Others"))
temp_final_gene<-merge(final_count,simulation_count,by="count")

no_biallelic<-temp_final_gene%>%subset((Maternal_specific+Paternal_specific)>=cutoff)#349,包括所有的OR gene和chrX gene
no_biallelic_ratio<-final_count_ratio%>%subset(gene%in%no_biallelic$gene)
write.xlsx(no_biallelic_ratio,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/no_biallelic_gene.xlsx")

candicate_gene<-no_biallelic%>%subset(type=="Others"&gene!="Peg3"&gene!="Kcnq1ot1")#281,去除印记基因
candicate_gene_ratio<-final_count_ratio%>%subset(gene%in%candicate_gene$gene)%>%mutate(Parental_specific=Maternal_specific+Paternal_specific)
candicate_gene_ratio<-candicate_gene_ratio%>%mutate(PM_type=case_when(Paternal_specific==0~"Only_Ma",Maternal_specific==0~"Only_Pa",Paternal_specific>0&Maternal_specific>0~"BI_RME"))
table(candicate_gene_ratio$PM_type)
#  BI_RME Only_Ma Only_Pa 
#     257      11      13 


# Figure 2C
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
no_biallelic_gene_p<-ggplot() +
  geom_point(data=candicate_gene_ratio, aes(x = Maternal_specific, y =Paternal_specific,color=PM_type,fill=PM_type),size=1,alpha=0.9) +
  geom_point(data=no_biallelic_ratio%>%subset(type!="Others"),aes(x = Maternal_specific, y =Paternal_specific,color=type,fill=type),size=1,alpha=0.9)+
  scale_fill_manual(values=c("#996277","#364f5e","#597b83","#828e79","#d8bf87"))+
  scale_color_manual(values=c("#828e79","#d8bf87","#996277","#364f5e","#597b83"))+
  scale_x_continuous(limits=c(0,1))+scale_y_continuous(limits=c(0,1))+
  labs(x="% M group",
       y="% P group")+
  theme_classic()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="None",text=element_text(size=10)
)  
candicate_gene_p<-ggMarginal(candicate_gene_p,groupColour=TRUE,groupFill=TRUE,type = "density",alpha = 0.4)
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/nobiallelic_gene.pdf",no_biallelic_count_p+no_biallelic_gene_p+plot_layout(ncol=2,widths=c(1,1.5)), device = "pdf", width = 24, height = 7, units = "cm", dpi = 300,bg="white")


# get candicate RME to filtered by allelic signals around the promoter regions
candicate_RME<-candicate_gene_ratio%>%subset(PM_type=="BI_RME")
# load the fragment file to identify the detail signal of promoter region
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
  signal_max<-which.max(final[[1]])
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
for(i in candicate_RME$gene){
  temp_promoter<-get_max_promoter(i)
  max_promoter_data<-rbind(max_promoter_data,temp_promoter)
}
dist_data<-data.frame()
for(i in candicate_RME$gene[180:257]){
  temp_dis<-get_dis(i)
  dist_data<-rbind(dist_data,temp_dis)
}

save(dist_data,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/dist.RData")
save(max_promoter_data,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/max_promoter_data.RData")

dist_data<-dist_data%>%subset(gene!="Peg3"&gene!="Kcnq1ot1") # because one of the paternal imprinting genes has maternal specific cell, remove it 
dist_data<-dist_data%>%mutate(p_type=case_when(cor_p<0.05~"P<0.05",cor_p>0.05|is.na(cor_p)~"P>0.05|P=NA"))
dist_data<-merge(dist_data,candicate_RME%>%select(gene,count,Biallelic))

dist_data$fc_signacl<-(dist_data$p_signal+10)/(dist_data$m_signal+10)
dist_data$fc_signacl[dist_data$fc_signacl>4]<-4 

# Figure 2D
dist_p<-ggplot() +geom_point(data=dist_data, aes(y = fc_signacl,x=log10(count*(1-Biallelic)),color=cor_r),alpha=0.8) +
  scale_color_gradientn(colors =c("#FCFfF5","#D1Dbbd","#3E606F","#193441"))+
  theme_classic()+theme(plot.title = element_text(hjust = 0.5), text=element_text(size=6), axis.text = element_text(colour = 'black',size = 6))+
  geom_hline(yintercept= c(1.5,0.7),linewidth=0.25,linetype="dashed")+geom_vline(xintercept= c(log10(25)),linewidth=0.25,linetype="dashed")+
  ylab("P group signal / M group signal")+
  xlab(" log10 ( # monoallelic expression cell )")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/dist.pdf",dist_p, device = "pdf", width = 8, height = 6, units = "cm", dpi = 300,bg="white")

# get the gene with similar promoter accessibility signal 
dist_data<-dist_data%>%mutate(dist_type=case_when(count*(1-Biallelic)>=25&fc_signacl<=1.5&fc_signacl>=0.7~"high",count*(1-Biallelic)<25|fc_signacl>1.5|fc_signacl<0.7~"low"))
dist_data%>%group_by(dist_type)%>%summarise(count=n())
#   dist_type count
#   <chr>     <int>
# 1 high        42
# 2 low         215

# get high_confident RME
high<-dist_data%>%subset(dist_type=="high")
RME<-candicate_RME%>%subset(gene%in%high$gene)
RME<-merge(RME,dist_data,by="gene")
write.xlsx(RME,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME.xlsx")
save(RME,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/merge/RME.RData")

# RME rank
# Figure 2E
library(ggprism)
RME_reads_Pa<-final_gene_dscore%>%subset(gene %in%RME$gene&merge_dscore_type=="Paternal_specific")%>%group_by(gene)%>%summarise(Pa_allele_rna=sum(allele_exp),Pa_allele_atac=sum(allele_frag),Pa_allele=sum(allele_exp)+sum(allele_frag))
RME_reads_Ma<-final_gene_dscore%>%subset(gene %in%RME$gene&merge_dscore_type=="Maternal_specific")%>%group_by(gene)%>%summarise(Ma_allele_rna=sum(allele_exp),Ma_allele_atac=sum(allele_frag),Ma_allele=sum(allele_exp)+sum(allele_frag))
RME<-merge(RME,RME_reads_Pa,by="gene")
RME<-merge(RME,RME_reads_Ma,by="gene")
RME<-RME%>%arrange(desc(count),desc(Ma_allele+Pa_allele),desc(Parental_specific),abs((Pa_allele_atac+Ma_allele_atac)/(Ma_allele+Pa_allele)-0.5))
RME$gene<-factor(RME$gene,levels=RME$gene)
RME_rank_p<-ggplot(RME)+
    geom_segment(aes(y=Paternal_specific,yend=Maternal_specific*(-1),
                   x=gene,xend=gene),
                   color="grey",size=0.5)+
    geom_point(aes(x=gene,y=Paternal_specific,size=Paternal_specific*count,color=log10(Pa_allele)))+
    geom_point(aes(x=gene,y=Maternal_specific*(-1),size=Maternal_specific*count,color=log10(Ma_allele)))+
    geom_hline(yintercept = 0, lty=2,color = 'black', lwd=0.2) +
   scale_y_continuous(breaks=seq(0,0.3,0.1), labels = as.character(seq(0,0.3,0.1)),name="% Paternal specific cell ratio",
                    sec.axis=dup_axis(breaks = seq(-0.3,0,0.1),name="% Maternal specific cell",labels = as.character(abs(seq(-0.3,0,0.1)))),
                    limits = c(-0.3, 0.3))+theme_prism(palette = "pearl",
              base_size = 8,
              base_line_size = 0.2,axis_text_angle = 45)+
    scale_color_gradientn(colors =c("#E3CAC3","#b9c8da","#8292a1","857f86"))+
    theme(text=element_text(size=10),
      plot.title = element_text(hjust = 0.5,vjust = 0.5,size=10), 
      axis.text.y=element_text(size=10,color = "black"),
      axis.text.x=element_text(size=10,  color = "black"))    
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME_rank.pdf",RME_rank_p, device = "pdf", width = 20, height = 8, units = "cm", dpi = 300,bg="white")



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

# Figure S2E-F
for(i in high$gene){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/high/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=6,ncol=2,heights = c(1,1.2,1.2,1.2,0.6,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}


