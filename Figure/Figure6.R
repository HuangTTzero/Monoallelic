# Analysis the Epha7 gene
# the promoter track plot
# Figure 6a
for(i in c("Epha7")){
  file=str_c("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/RME/high/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=8,ncol=2,heights = c(1,1.2,1.2,1.2,0.4,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}
write.xlsx(final_gene_dscore%>%subset(gene=="Epha7"),file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Figure6a.xlsx")
# result of Fish
# Figure 6b

Epha7<-data.frame(allele_exp_type=c("P","P","P","M","M","M","B","B","B"),rep=c("rep1","rep2","rep3","rep1","rep2","rep3","rep1","rep2","rep3"),
                ratio=c(29.5,32.6,41.7,39.3,34.9,27.1,31.1,32.6,31.2),
                method=c(rep("RNA_Fish",9)))
Epha7$method<-factor(Epha7$method,levels=c("RNA_Fish"))
Epha7_p<-ggplot(Epha7,aes(x=method,y=ratio,fill=allele_exp_type,color=allele_exp_type))+
    geom_bar(color="black",stat="summary",fun=mean,position="dodge",size=0)+
    geom_point(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9),
             size = 1, shape = 21, stroke = 0.2,color="black") +
    stat_summary(fun.data='mean_sd',geom="errorbar",linewidth=0.15,size=1,width = 0.5,position = position_dodge(0.9))+
    scale_fill_manual(values=c("#dfe2d6","#94b2c2","#c68989"))+
    scale_color_manual(values=c("black","black","black"))+
    scale_y_continuous(limits=c(0,100))+
    labs(y="% cell",
       x="")+
    theme_classic()+
    theme(
        legend.title = element_blank(),
        text=element_text(size=6)
    )
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/Fish/experimental/Epha7_Fishresult.pdf",Epha7_p, device = "pdf", width = 4.5, height = 2, units = "cm", dpi = 60,bg="white")
write.xlsx(Epha7,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Figure6b.xlsx")


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
# Figure 6c
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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Epha7_OR_sample30.pdf",Epha7_OR_sample30, device = "pdf", width = 40, height = 10, units = "cm", dpi = 300)
expression_data <- FetchData(object = filtered_mature_OSNs, 
                             vars = c("Epha7", "expressed_OR"))
write.xlsx(expression_data, 
           file = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Figure6c.xlsx",
           rowNames = TRUE)# barplot of Epha7 expression of OR subgroups
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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Epha7_OR_meanexp.pdf",Epha7_OR_meanexp_p, device = "pdf", width = 10, height = 5, units = "cm", dpi = 300,bg="white")
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
# Figure 6d
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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Epha7_OR_projection.pdf",plot_features(OR_group_up)+plot_layout(ncol=2,width=c(1,1.5)), device = "pdf", width = 25, height = 10, units = "cm", dpi = 300)


# DV and AP scores of genes with two expression patterns at different expression levels
# Figure 6e
OR_group_up_with_OB<-OR_group_up%>%subset(!is.na(x) & !is.na(y))
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
              stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = list(c("high", "low")),label.y = c(0.95))+
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
              stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = list(c("high", "low")),label.y = c(0.95))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Epha7_OR_projection_vln.pdf",Epha7_DV_p+Epha7_AP_p, device = "pdf", width = 18, height = 10, units = "cm", dpi = 300)
# 
stat_results_DV <- OR_group_up_with_OB %>%
  rstatix::wilcox_test(
    scale_y ~ avg_exp_type,
    comparisons = list(c("high", "low")),
    exact = FALSE,
    conf.level = 0.95,
    detailed = TRUE
  )

#
effect_sizes_DV <- OR_group_up_with_OB %>%
  rstatix::wilcox_effsize(
    scale_y ~ avg_exp_type,
    comparisons = list(c("high", "low"))
  )

# 
final_results_DV <- stat_results_DV %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes_DV %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

print("DV score (scale_y) 统计结果：")
print(final_results_DV)

# 2.
stat_results_AP <- OR_group_up_with_OB %>%
  rstatix::wilcox_test(
    scale_x ~ avg_exp_type,
    comparisons = list(c("high", "low")),
    exact = FALSE,
    conf.level = 0.95,
    detailed = TRUE
  )

# 
effect_sizes_AP <- OR_group_up_with_OB %>%
  rstatix::wilcox_effsize(
    scale_x ~ avg_exp_type,
    comparisons = list(c("high", "low"))
  )

# 
final_results_AP <- stat_results_AP %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes_AP %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

print("AP score (scale_x) 统计结果：")
print(final_results_AP)
# [1] "DV score (scale_y) 统计结果："
# # A tibble: 1 × 10
#   group1 group2    n1    n2 statistic     p conf.low conf.high effsize magnitude
#   <chr>  <chr>  <int> <int>     <dbl> <dbl>    <dbl>     <dbl>   <dbl> <ord>    
# 1 low    high     148    35      2480 0.698  -0.0525    0.0357  0.0289 small    
# [1] "AP score (scale_x) 统计结果："
# # A tibble: 1 × 10
#   group1 group2    n1    n2 statistic      p conf.low conf.high effsize magnitude
#   <chr>  <chr>  <int> <int>     <dbl>  <dbl>    <dbl>     <dbl>   <dbl> <ord>    
# 1 low    high     148    35     3262. 0.0173   0.0139     0.137   0.176 small  

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
ggsave(filename =  "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Epha7_expVScell.pdf",Epha7_expVSBicell_p+Epha7_expVSMonocell_p, device = "pdf", width =15, height = 8, units = "cm", dpi = 300,bg="white")

# Figure 7d,7e
# selected OR subgroups with more than three Epha7 allelic informative cells
# compare the foldchange between high and low OR subgroups
ORgroup_Epha7_fc <- ggplot(Epha7_dscore_in_ORgroup %>% subset(total >= 4), 
                           aes(x = avg_exp_type, y = fc, fill = avg_exp_type)) + 
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # 添加抖动点
  stat_compare_means(method = "wilcox.test", label = "p.signif", 
                     comparisons = list(c("low", "high"))) +
  scale_fill_manual(values = c("#cdd2db", "#45486f")) +
  NoLegend() + 
  labs(x = NULL, y = NULL) + 
  theme_classic() +
  theme(text = element_text(size = 6))
ggsave(filename =  "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/ORgroup_Epha7_fc.pdf",ORgroup_Epha7_fc, device = "pdf", width =15, height = 8, units = "cm", dpi = 300,bg="white")


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
                  stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = list(c("high", "low")))+
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
              stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = list(c("high", "low")))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/ORgroup_Epha7_allele.pdf",ORgroup_Epha7_Bicell+ORgroup_Epha7_Monocell+ORgroup_Epha7_fc, device = "pdf", width = 24, height = 7, units = "cm", dpi = 300)

# 1. 双等位基因表达率的统计检验
stat_results_Bi <- Epha7_dscore_in_ORgroup %>%
  rstatix::wilcox_test(
    Bi_rate ~ avg_exp_type,
    comparisons = list(c("high", "low")),
    exact = FALSE,
    conf.level = 0.95,
    detailed = TRUE
  )

# 计算效应量
effect_sizes_Bi <- Epha7_dscore_in_ORgroup %>%
  rstatix::wilcox_effsize(
    Bi_rate ~ avg_exp_type,
    comparisons = list(c("high", "low"))
  )

# 合并结果
final_results_Bi <- stat_results_Bi %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes_Bi %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

print("双等位基因表达率 (Bi_rate) 统计结果：")
print(final_results_Bi)

# 2. 单等位基因表达率的统计检验
stat_results_Mono <- Epha7_dscore_in_ORgroup %>%
  rstatix::wilcox_test(
    Mono_rate ~ avg_exp_type,
    comparisons = list(c("high", "low")),
    exact = FALSE,
    conf.level = 0.95,
    detailed = TRUE
  )

# 计算效应量
effect_sizes_Mono <- Epha7_dscore_in_ORgroup %>%
  rstatix::wilcox_effsize(
    Mono_rate ~ avg_exp_type,
    comparisons = list(c("high", "low"))
  )

# 合并结果
final_results_Mono <- stat_results_Mono %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes_Mono %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

print("单等位基因表达率 (Mono_rate) 统计结果：")
print(final_results_Mono)

# [1] "双等位基因表达率 (Bi_rate) 统计结果："
# # A tibble: 1 × 10
#   group1 group2    n1    n2 statistic          p conf.low conf.high effsize magnitude
#   <chr>  <chr>  <int> <int>     <dbl>      <dbl>    <dbl>     <dbl>   <dbl> <ord>    
# 1 low    high     154    37     1424. 0.00000105  -0.0271   -0.0133   0.353 moderate 
# [1] "单等位基因表达率 (Mono_rate) 统计结果："
# # A tibble: 1 × 10
#   group1 group2    n1    n2 statistic      p   conf.low conf.high effsize magnitude
#   <chr>  <chr>  <int> <int>     <dbl>  <dbl>      <dbl>     <dbl>   <dbl> <ord>    
# 1 low    high     154    37     2424. 0.0505 -0.0000237 0.0000400   0.142 small   

# 3. fc的统计检验
stat_results_fc <- Epha7_dscore_in_ORgroup %>% subset(Epha7_allele_cells >= 4) %>%
  rstatix::wilcox_test(
    fc ~ avg_exp_type,
    comparisons = list(c("high", "low")),
    exact = FALSE,
    conf.level = 0.95,
    detailed = TRUE
  )

# 计算效应量
effect_sizes_fc <- Epha7_dscore_in_ORgroup %>% subset(Epha7_allele_cells >= 4) %>%
  rstatix::wilcox_effsize(
    fc ~ avg_exp_type,
    comparisons = list(c("high", "low"))
  )

# 合并结果
final_results_fc <- stat_results_fc %>%
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  left_join(effect_sizes_fc %>% select(group1, group2, effsize, magnitude), 
            by = c("group1", "group2"))

print("单等位基因表达率 (fc_rate) 统计结果：")
print(final_results_fc)
# [1] "单等位基因表达率 (fc_rate) 统计结果："
# # A tibble: 1 × 10
#   group1 group2    n1    n2 statistic     p   conf.low conf.high effsize magnitude
#   <chr>  <chr>  <int> <int>     <dbl> <dbl>      <dbl>     <dbl>   <dbl> <ord>    
# 1 high   low       11     7      54.5 0.143 -0.0000341      7.18   0.356 moderate 
names(Epha7_dscore_in_ORgroup)[names(Epha7_dscore_in_ORgroup) == "total"] <- "Epha7_allele_cells"
names(Epha7_dscore_in_ORgroup)[names(Epha7_dscore_in_ORgroup) == "count"] <- "cellnumber"
write.xlsx(Epha7_dscore_in_ORgroup%>%select(-scale_x,-scale_y,-fc,-fc1,-Bi_rate,-Mono_rate,-avg_exp_rank,x_rank,y_rank),file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Figure6d-g.xlsx")
Epha7_dscore_in_ORgroup<-read.xlsx("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure6/Figure6d-g.xlsx")