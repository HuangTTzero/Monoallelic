# Obtain the average expression level of the RMEs in different OR subgroups
high_RME_OR_up_exp<-AverageExpression(
  filtered_mature_OSNs,
  assay="RNA",
  layer="data",features=high_RME$gene,
  group.by = "expressed_OR"
)
high_RME_OR_up_exp<-high_RME_OR_up_exp$RNA
high_RME_OR_up_exp<-as.data.frame(high_RME_OR_up_exp)
high_RME_OR_up_exp$gene<-rownames(high_RME_OR_up_exp)
high_RME_OR_up_exp<-high_RME_OR_up_exp%>%pivot_longer(cols = -c("gene"),names_to="expressed_OR",values_to="avg_exp")
high_RME_OR_up_exp<-high_RME_OR_up_exp%>%subset(expressed_OR %in%OR_group_up$expressed_OR)

# CV = (Mean / Standard deviation) * 100%
high_RME_OR_CV<-data.frame(gene=c(),mean=c(),sd=c(),cv=c())
for(i in high_RME$gene){
    temp_data<-high_RME_OR_up_exp%>%subset(gene==i)
    OR_mean<-mean(temp_data$avg_exp)
    OR_sd<-sd(temp_data$avg_exp)
    high_RME_OR_CV<-rbind(high_RME_OR_CV,data.frame(gene=i,mean=OR_mean,sd=OR_sd,cv=OR_sd/OR_mean*100))
}
high_RME_OR_CV<-high_RME_OR_CV%>%arrange(cv)
high_RME_OR_CV$gene<-factor(high_RME_OR_CV$gene,levels=high_RME_OR_CV$gene)

# vlnplot_折线图
# Figure 8a
high_RME_OR_up_exp$gene<-factor(high_RME_OR_up_exp$gene,levels=high_RME_OR_CV$gene)
high_RME_OR_up_exp_p<-ggviolin(high_RME_OR_up_exp, x = "gene", y = "avg_exp", fill="#789197",size=0.2,
              trim=TRUE,
              alpha = 1,width = 1,
              legend = "none",#去掉legend
              xlab="", ylab="",
              add = "boxplot") +
  # 绘制每个基因的变异程度（折线图），使用均值和标准差数据
  geom_line(data = high_RME_OR_CV, aes(x = gene, y = cv/2.5,group=1), 
            color = "grey", size = 0.5) +  # 使用*2缩放标准差，以便与小提琴图匹配
  geom_point(data = high_RME_OR_CV, aes(x = gene, y = cv/2.5,), 
             shape = 17, color = "#192d77", size = 1) +
  # 定义Y轴：主Y轴（小提琴图）和次Y轴（折线图）
  scale_y_continuous(name = "average expression", 
                     sec.axis = sec_axis(~.*2.5, name = "Coefficient of Variation")) +  # 标准差缩放
  labs(title = "",
       x = "")  +
  theme_classic()+ 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text=element_text(size=8),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/high_RME_OR_CV.pdf",high_RME_OR_up_exp_p, device = "pdf", width = 15, height = 6, units = "cm", dpi = 300)



# Obtain the average expression level of the RMEs in sampled different OR subgroups
set.seed(123)
expressed_OR_sample<-sample(filtered_mature_OSNs_up$expressed_OR,size=length(filtered_mature_OSNs_up$expressed_OR))
names(expressed_OR_sample)<-NULL
filtered_mature_OSNs_up$expressed_OR_sample<-expressed_OR_sample

high_RME_OR_up_exp_sample<-AverageExpression(
  filtered_mature_OSNs_up,
  assay="RNA",
  layer="data",features=high_RME$gene,
  group.by = "expressed_OR_sample"
)
high_RME_OR_up_exp_sample<-high_RME_OR_up_exp_sample$RNA
high_RME_OR_up_exp_sample<-as.data.frame(high_RME_OR_up_exp_sample)
high_RME_OR_up_exp_sample$gene<-rownames(high_RME_OR_up_exp_sample)
high_RME_OR_up_exp_sample<-high_RME_OR_up_exp_sample%>%pivot_longer(cols = -c("gene"),names_to="expressed_OR_sample",values_to="avg_exp_sample")

#CV=均值/标准差*100%
high_RME_OR_CV_sample<-data.frame(gene=c(),mean=c(),sd=c(),cv=c())
for(i in high_RME$gene){
    temp_data<-high_RME_OR_up_exp_sample%>%subset(gene==i)
    OR_mean<-mean(temp_data$avg_exp_sample)
    OR_sd<-sd(temp_data$avg_exp_sample)
    high_RME_OR_CV_sample<-rbind(high_RME_OR_CV_sample,data.frame(gene=i,mean_sample=OR_mean,sd_sample=OR_sd,cv_sample=OR_sd/OR_mean*100))
}
high_RME_OR_CV_sample$gene<-factor(high_RME_OR_CV_sample$gene,levels=high_RME_OR_CV$gene)

high_RME_OR_CV<-merge(high_RME_OR_CV,high_RME_OR_CV_sample)

# F_test: compare the RME expression in sampled obj with  truely obj
f_test<-c()
for(i in high_RME_OR_CV$gene){
    temp_data<-high_RME_OR_up_exp%>%subset(gene==i)
    temp_data_sampe<-high_RME_OR_up_exp_sample%>%subset(gene==i)
    test<-var.test(temp_data$avg_exp, temp_data_sampe$avg_exp_sample)
    f_test<-c(f_test,test$p.value)
}
high_RME_OR_CV$f_test<-f_test

high_RME_OR_CV<-high_RME_OR_CV%>%arrange(cv)
high_RME_OR_CV$gene<-factor(high_RME_OR_CV$gene,levels=high_RME_OR_CV$gene)

# RME_OR VS RME_OR_up_exp_sample
# Figure 8b
CV_violin<-function(genename){
  data<-high_RME_OR_up_exp%>%subset(gene==genename)
  data$group<-"Observed"
  data_sample<-high_RME_OR_up_exp_sample%>%subset(gene==genename)
  data_sample$group<-"Random"
  data_sample$avg_exp<-data_sample$avg_exp_sample
  final_data<-rbind(data%>%select(gene,avg_exp,group),data_sample%>%select(gene,avg_exp,group))
  boxplot_p<-ggviolin(final_data, x = "group", y = "avg_exp", fill = "#914350", 
              alpha = 1,width = 0.5,
              trim=TRUE,
              legend = "none",#去掉legend
              xlab="", ylab="avg_exp between OR groups",
              title=genename,
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=6))
  print(boxplot_p)
}

pdf("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/CV_violin.pdf",height=2.5,width=3.5)
DS_gene<-high_RME_OR_CV%>%subset(cv>=53.48619&f_test<0.05) # 5 genes remained
for(gene in DS_gene$gene){
  CV_violin(gene)
}
for(gene in "Adk"){
  CV_violin(gene)
}
dev.off()


# Haploinsufficient and Triplosensitive probability of RME genes
# Figure 8c
gene_disease<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/ds_gene/mmc6.xlsx")[,c(15:17,19)]
colnames(gene_disease)<-c("top_gene","high_confidence","confidence","HP")
gene<-c(gene_disease$high_confidence,gene_disease$confidence)
gene<-gene[!gene=="-"]
gene<-strsplit(gene,split=";")
gene<-sort(unlist(gene)

Human_Mouse<-read.table(file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/ds_gene/Human_Mouse.txt",sep="\t",header=TRUE)[,c(4,10)]
Human_Mouse_subset<-Human_Mouse%>%subset(Gene.name%in%gene)%>%arrange(Gene.name)

Dsgene<-read.xlsx("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/ds_gene/mmc7.xlsx")[,c(1:3)]
Human_Mouse_RME<-Human_Mouse%>%subset(Mouse.gene.name%in%RME$gene)
Dsgene_RME<-merge(Dsgene,Human_Mouse_RME,by.x="Gene",by.y="Gene.name")
Dsgene_RME<-Dsgene_RME%>%mutate(type=case_when(pHaplo>= 0.86&pTriplo<0.94~"haploinsufficient",pTriplo >= 0.94&pHaplo< 0.86~"triplosensitive",pHaplo >= 0.86&pTriplo >= 0.94~"Both"))
write.xlsx(Dsgene_RME,file="/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/ds_gene/Dsgene_RME.xlsx")


Dsgene_RME<-Dsgene_RME%>%arrange(type,desc(trunc(pHaplo * 100) / 100),desc(trunc(pTriplo * 100) / 100))
Dsgene_RME_longer<-Dsgene_RME%>%pivot_longer(cols = -c("Gene","Mouse.gene.name","type"),names_to="score_name",values_to="score")
Dsgene_RME_longer$Mouse.gene.name<-factor(Dsgene_RME_longer$Mouse.gene.name,levels=rev(Dsgene_RME$Mouse.gene.name))

# Dsgene dotplot
# Custom transformation function: Power transformation makes points closer to 1 larger.
custom_trans <- trans_new(
  name = "custom",
  transform = function(x) x^5,   
  inverse = function(x) x^(1/5)
)

Dsgene_dot<-ggplot(data=Dsgene_RME_longer,aes(y=Mouse.gene.name,x=score_name,size=score,fill=score_name,color=score_name))+
  geom_point() +scale_fill_manual(values=c("#7d3d92","#c19cc8"))+scale_color_manual(values=c("#7d3d92","#c19cc8"))+
  labs(title = "",
       x = "",y="")  +scale_size_continuous(range = c(0, 3),trans = custom_trans, limits = c(0, 1),
                        breaks = seq(0, 1, by = 0.25))+
  theme_classic()+ 
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank(),
        text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure5/ds_gene/Dsgene_dot.pdf",Dsgene_dot, device = "pdf", width = 8, height =12, units = "cm", dpi = 300)
