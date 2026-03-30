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

ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/high_RME_OR_CV_new.pdf",high_RME_OR_up_exp_p, device = "pdf", width = 13, height = 7, units = "cm", dpi = 300)



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
f_statistic <- c()
df1 <- c()
df2 <- c()
p_value <- c()
ci_lower <- c()
ci_upper <- c()
variance_ratio <- c()

# 循环计算
for(i in high_RME_OR_CV$gene) {
    temp_data <- high_RME_OR_up_exp %>% subset(gene == i)
    temp_data_sampe <- high_RME_OR_up_exp_sample %>% subset(gene == i)
    
    test <- var.test(temp_data$avg_exp, temp_data_sampe$avg_exp_sample)
    
    # 存储各个统计量
    f_statistic <- c(f_statistic, test$statistic)
    df1 <- c(df1, test$parameter[1])
    df2 <- c(df2, test$parameter[2])
    p_value <- c(p_value, test$p.value)
    ci_lower <- c(ci_lower, test$conf.int[1])
    ci_upper <- c(ci_upper, test$conf.int[2])
    variance_ratio <- c(variance_ratio, test$estimate)
}

# 一次性添加到数据框
high_RME_OR_CV$f_statistic <- f_statistic
high_RME_OR_CV$df1 <- df1
high_RME_OR_CV$df2 <- df2
high_RME_OR_CV$f_test_p <- p_value  # 保留您原来的列名
high_RME_OR_CV$ci_lower <- ci_lower
high_RME_OR_CV$ci_upper <- ci_upper
high_RME_OR_CV$variance_ratio <- variance_ratio

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

pdf("/data/R02/huangtt39/F1_OSN/new_analysis/CV_violin.pdf",height=2.5,width=3.5)
DS_gene<-high_RME_OR_CV%>%subset(cv>=53.48619&f_test<0.05) # 5 genes remained
for(gene in DS_gene$gene){
  CV_violin(gene)
}
dev.off()

write.xlsx(high_RME_OR_CV%>%select(-f_test),file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/Figure8ab.xlsx")

#Figure8b
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

base_path = Path("/data/R02/tanj93/project/RME/human_brain/")

h5ad_file = base_path / "WHB-10Xv3-Neurons-log2.h5ad"
metadata_file = base_path / "cell_metadata_with_celltype.csv"

genes_of_interest = [
"AGO3","BMP6","CAMKMT","CHCHD2","CLEC16A","CLIP4",
"DMXL2","EPHA7","FETUB","HECW1","IMMP2L","KCNH7",
"LARGE1","LINC-PINT","MACROD2","MEGF9","MGST3","PLCL1",
"SESN3","SORBS1","SYN2","TAX1BP1","WDR17","WWOX"
]

adata = sc.read_h5ad(h5ad_file)

print("Cells:", adata.n_obs)
print("Genes:", adata.n_vars)

gene_mask = adata.var["gene_symbol"].isin(genes_of_interest)
adata_sub = adata[:, gene_mask].copy()
gene_names = adata_sub.var["gene_symbol"].values
print("Genes found:", len(gene_names))

expr = adata_sub.X

if hasattr(expr, "toarray"):
    expr = expr.toarray()

expr_df = pd.DataFrame(
    expr,
    columns=gene_names,
    index=adata_sub.obs_names
)

expr_df["barcode"] = expr_df.index.str.split(":").str[1]
expr_df = expr_df.set_index("barcode")

meta = pd.read_csv(metadata_file)
meta = meta.set_index("cell_barcode")
print("Metadata cells:", len(meta))

merged = meta.join(expr_df, how="inner")
print("Merged cells:", len(merged))

genes = genes_of_interest

pct_region = merged.groupby("anatomical_division_label")[genes].apply(
    lambda x: (x > 0).mean() * 100
)
avg_region = merged.groupby("anatomical_division_label")[genes].apply(
    lambda x: x.where(x > 0).mean()
)

pct_celltype = merged.groupby("celltype")[genes].apply(
    lambda x: (x > 0).mean() * 100
)
avg_celltype = merged.groupby("celltype")[genes].apply(
    lambda x: x.where(x > 0).mean()
)

gene_order = avg_region.mean().sort_values(ascending=False).index.tolist()

pct_region.to_csv(base_path / "pct_region.csv")
avg_region.to_csv(base_path / "avg_region.csv")

pct_celltype.to_csv(base_path / "pct_celltype.csv")
avg_celltype.to_csv(base_path / "avg_celltype.csv")

colors = ["#FCFFF5","#D1DBBD","#3E606F","#193441"]
cmap = plt.matplotlib.colors.LinearSegmentedColormap.from_list("custom", colors)

def make_dotplot(pct, avg, x_labels, outfile, figsize=(6,8)):
    pct = pct[gene_order]
    avg = avg[gene_order]
    pct = pct.T
    avg = avg.T
    fig, ax = plt.subplots(figsize=figsize)
    x = np.arange(len(x_labels))
    y = np.arange(len(gene_order))
    X, Y = np.meshgrid(x, y)
    sizes = pct.values.flatten()
    colors_vals = avg.values.flatten()
    sizes = np.nan_to_num(sizes)
    colors_vals = np.nan_to_num(colors_vals)
    sizes_scaled = np.sqrt(sizes/100) * 280
    scatter = ax.scatter(
        X.flatten(),
        Y.flatten(),
        s=sizes_scaled,
        c=colors_vals,
        cmap=cmap,
        vmin=0,
        vmax=12,
        edgecolors="gray",
        linewidths=0.5
    )
    ax.set_xticks(x)
    ax.set_xticklabels(x_labels, rotation=45, ha="right")
    ax.set_yticks(y)
    ax.set_yticklabels(gene_order, fontstyle="italic")
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.4)
    cbar.set_label("Average Expression")
    legend_sizes = [25,50,75,100]
    legend_points = [
        Line2D(
            [0],[0],
            marker='o',
            color='w',
            markerfacecolor='black',
            markersize=np.sqrt(s/100)*15
        )
        for s in legend_sizes
    ]
    ax.legend(
        legend_points,
        [f"{s}%" for s in legend_sizes],
        title="Percent expressed",
        bbox_to_anchor=(1.3,1),
        frameon=False
    )
    plt.tight_layout(rect=[0,0,0.85,1])
    plt.savefig(outfile, dpi=300)
    print("Saved:", outfile)

make_dotplot(
    pct_region,
    avg_region,
    pct_region.index.tolist(),
    base_path / "dotplot_brain_region.pdf",
    figsize=(8,8)
)

make_dotplot(
    pct_celltype,
    avg_celltype,
    pct_celltype.index.tolist(),
    base_path / "dotplot_celltype.pdf",
    figsize=(10,8)
)

# Haploinsufficient and Triplosensitive probability of high_RME genes
# Figure 8d
gene_disease<-read.xlsx("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/ds_gene/mmc6.xlsx")[,c(15:17,19)]
colnames(gene_disease)<-c("top_gene","high_confidence","confidence","HP")
gene<-c(gene_disease$high_confidence,gene_disease$confidence)
gene<-gene[!gene=="-"]
gene<-strsplit(gene,split=";")
gene<-sort(unlist(gene))

Human_Mouse<-read.table(file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/ds_gene/Human_Mouse.txt",sep="\t",header=TRUE)[,c(4,10)]
Human_Mouse_subset<-Human_Mouse%>%subset(Gene.name%in%gene)%>%arrange(Gene.name)

Dsgene<-read.xlsx("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/ds_gene/mmc7.xlsx")[,c(1:3)]
Human_Mouse_RME<-Human_Mouse%>%subset(Mouse.gene.name%in%high_RME$gene)
Dsgene_RME<-merge(Dsgene,Human_Mouse_RME,by.x="Gene",by.y="Gene.name")
Dsgene_RME<-Dsgene_RME%>%mutate(type=case_when(pHaplo>= 0.86&pTriplo<0.94~"haploinsufficient",pTriplo >= 0.94&pHaplo< 0.86~"triplosensitive",pHaplo >= 0.86&pTriplo >= 0.94~"Both"))
write.xlsx(Dsgene_RME,file="/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/Figure8d.xlsx")


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
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/Figure8/Dsgene_dot.pdf",Dsgene_dot, device = "pdf", width = 8, height =8, units = "cm", dpi = 300)
