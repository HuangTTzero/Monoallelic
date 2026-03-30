high_RME<-high_RME%>%arrange(desc(1-Biallelic), desc(Maternal_specific)) 

LiMCA_high_RME_ratio<-LiMCA_10_count_ratio%>%subset(gene%in%high_RME$gene) # 23 genes ,the Bamp6 has no allelic reads
LiMCA_high_RME_ratio$gene<-factor(LiMCA_high_RME_ratio$gene,levels=high_RME$gene)
LiMCA_high_RME_count<-LiMCA_10_count%>%subset(gene%in%high_RME$gene)

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
            file = "/data/R02/huangtt39/F1_OSN/LiMCA/HIC/LiMCA_high_RME_dscore.txt",
            quote = FALSE, 
            row.names = FALSE,sep="\t")
write.table(LiMCA2Hic_dscore, 
            file = "/data/R02/huangtt39/F1_OSN/LiMCA/HIC/LiMCA2Hic_dscore.txt",
            quote = FALSE, 
            row.names = FALSE,sep="\t")
LiMCA2Hic_promoter<-max_promoter_data%>%subset(gene%in%LiMCA2Hic_dscore$gene)
# TAD analysis
TAD<-read.table(file="/data/R02/huangtt39/F1_OSN/LiMCA/HIC/bulk/TAD/TAD_r5kb_contact_domains_list/5000_blocks.bedpe")
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
  gene_dir <- paste0("/data/R02/huangtt39/F1_OSN/LiMCA/HIC/merge_impute/", gene)
  
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
write.xlsx(final_df, "/data/R02/huangtt39/F1_OSN/LiMCA/HIC/merge_impute/TAD_total_contact_counts.xlsx", row.names = FALSE)

final_df$Maternal_specific_dscore<-final_df$Maternal_specific_PP/(final_df$Maternal_specific_MM+final_df$Maternal_specific_PP)-0.5
final_df$Paternal_specific_dscore<-final_df$Paternal_specific_PP/(final_df$Paternal_specific_MM+final_df$Paternal_specific_PP)-0.5
final_df<-final_df%>%mutate(total_counts=(Paternal_specific_PP+Paternal_specific_MM+Maternal_specific_PP+Maternal_specific_MM))
#Extended Data Figure6a
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
  "/data/R02/huangtt39/F1_OSN/analysis/figure/LiMCA/TAD_dscore_scatterplot.pdf",
  p,
  width = 5.5,
  height = 4,
  dpi = 300
)

#Extended Data Figure6b
library(org.Mm.eg.db)
library("plotgardener")
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(plotgardener)
hicFile <- "/data/R02/huangtt39/F1_OSN/LiMCA/HIC/bulk/4DNFIVUGNDD7.hic"
library(plotgardener)
library(strawr)

# 主绘图函数
library(plotgardener)
library(strawr)

# 创建PDF文件
pdf("/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/TAD_plots.pdf", 
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

output_pdf <- "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/all_genes_RNA_expression_log2.pdf"

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
  pdf_file <- paste0("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure/", 
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

# allele-specific Hi-C
plot_allele_specific_hic <- function(gene, chr, extended_start, extended_end, 
                                     sample_type, start_y, plot_x, plot_width) {
  
  gene_dir <- paste0("/data/R02/huangtt39/F1_OSN/LiMCA/HIC/merge_impute/", gene)
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
  
  cat("\nAll PDFs saved in:", "/data/R02/huangtt39/F1_OSN/LiMCA/HIC/merge_impute/\n")
}

# 使用函数
plot_all_genes_separate(final_df)


#Link peak
DefaultAssay(Combined_mOSN) <- "ATAC"
library(BSgenome.Mmusculus.UCSC.mm10)
# first compute the GC content for each peak
Combined_mOSN <- RegionStats(Combined_mOSN, genome = BSgenome.Mmusculus.UCSC.mm10)

# link peaks to genes
Combined_mOSN <- LinkPeaks(
  object = Combined_mOSN,
  peak.assay = "ATAC",
  expression.assay = "RNA",
  genes.use = final_gene$gene
)

# 2. 主处理函数
process_gene_link_peak <- function(gene_name, peak_annotation, final_gene_dscore) {
  # 提取该基因的细胞信息
  data <- subset(final_gene_dscore, gene == gene_name)
  
  # 创建细胞列表：G1和G2版本
  cells_g1 <- paste0("G1_", data$barcode)
  cells_g2 <- paste0("G2_", data$barcode)
  all_cells <- c(cells_g1, cells_g2)
  
  # 提取细胞子集
  temp_Obj <- subset(Obj, cells = all_cells)
  
  # 获取linked peaks
  peaks <- peak_annotation%>%subset(gene==gene_name)
  if (is.null(peaks) || nrow(peaks) == 0) return(NULL)
  
  # 创建细胞类型映射（确保正确对应）
  cell_mapping <- data.frame(
    cell = all_cells,
    original_barcode = rep(data$barcode, 2),
    group = rep(c("G1", "G2"), each = nrow(data))
  )
  
  # 从data中获取merge_dscore_type
  cell_mapping$merge_dscore_type <- data$merge_dscore_type[match(
    cell_mapping$original_barcode, data$barcode
  )]
  
  # 创建最终细胞类型标签
  cell_mapping$cell_type <- paste(cell_mapping$merge_dscore_type, 
                                  cell_mapping$group, sep = "_")

  # 处理每个peak
  results <- lapply(1:nrow(peaks), function(i) {
    peak<-peaks[i,]
    peak_gr <- StringToGRanges(peak$peak, sep = c("-", "-"))
    
    counts <- CountsInRegion(temp_Obj, "ATAC", peak_gr)
    
    # 创建计数数据框
    counts_df <- data.frame(
      cell = names(counts),
      count = as.numeric(counts)
    )
    
    # 添加细胞类型信息
    counts_df$cell_type <- cell_mapping$cell_type[match(counts_df$cell, cell_mapping$cell)]
    
    # 按细胞类型汇总
    summary <- aggregate(count ~ cell_type, counts_df, sum)
    
    # 转换为宽格式
    wide <- data.frame(gene = gene_name, peak = as.character(peak$peak))
    
    for (i in 1:nrow(summary)) {
      wide[[summary$cell_type[i]]] <- summary$count[i]
    }
    wide$peak_type<-peak$peak_type
    return(wide)
  })
  
  # 合并所有peak的结果
  result_df <- do.call(rbind, results)
  
  # 确保所有细胞类型列都存在
  expected_types <- c("Paternal_specific_G1", "Paternal_specific_G2",
                      "Maternal_specific_G1", "Maternal_specific_G2",
                      "Biallelic_G1", "Biallelic_G2")
  
  for (type in expected_types) {
    if (!type %in% colnames(result_df)) {
      result_df[[type]] <- 0
    }
  }
  
  return(result_df)
}

peak_annotation <- function(gene_name, Combined_mOSN, promoter_data) {
  promoter_gr <- GRanges(
    seqnames = promoter_data$chr[promoter_data$gene == gene_name],
    ranges = IRanges(
      start = promoter_data$start[promoter_data$gene == gene_name],
      end = promoter_data$end[promoter_data$gene == gene_name]
    )
  )
  
  all_promoters_gr <- GRanges(
    seqnames = promoter_data$chr,
    ranges = IRanges(start = promoter_data$start, end = promoter_data$end),
    gene = promoter_data$gene
  )

  peaks <- GetLinkedPeaks(Combined_mOSN, gene_name, "ATAC", min.abs.score = 0)
  if (is.null(peaks) || length(peaks) == 0) return(NULL)

  results <- lapply(peaks, function(p) {
    peak_gr <- StringToGRanges(p, sep = c("-", "-"))
    
    # 检查重叠
    overlaps <- findOverlaps(peak_gr, all_promoters_gr)
    
    if (length(overlaps) > 0) {
      overlapping_genes <- unique(all_promoters_gr$gene[subjectHits(overlaps)])
      if (gene_name %in% overlapping_genes) {
        peak_type <- "promoter"
        other_gene <- NA
      } else {
        peak_type <- "other_promoter"
        other_gene <- overlapping_genes[1]
      }
    } else {
      peak_type <- "enhancer"
      other_gene <- NA
    }
     
    data.frame(
      gene = gene_name,
      peak = p,
      peak_type = peak_type,
      other_gene = other_gene,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, results)
}
link_peaks <- Links(Combined_mOSN)
save(link_peaks,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/link_peaks.RData")

link_peaks_df <- as.data.frame(link_peaks)

# 只保留你需要的列
link_peaks_df_subset <- link_peaks_df[, c("gene", "peak", "zscore", "pvalue", "score")]

# 3. 批量处理所有基因
link_peak_anno <- list()
for (gene in final_gene$gene) {
  result <- peak_annotation(gene, Combined_mOSN, promoter_bed)
  if (!is.null(result)) {
    link_peak_anno [[gene]] <- result
  }
}
link_peak_anno_all <- do.call(rbind, link_peak_anno)
gene_list<-final_gene %>%subset(!gene%in%c("Peg3","Kcnq1ot1")&type=="Others")
link_peak_anno_all<-link_peak_anno_all%>%subset(gene%in%gene_list$gene)
link_peak_anno_all$gene_type<-"non-RME"
for (i in 1:nrow(link_peak_anno_all)) { 
  gene <- link_peak_anno_all$gene[i]
  if(gene %in% RME$gene ){
    type="RME"
    link_peak_anno_all[i,]$gene_type=type
  }
} 

save(link_peak_anno_all,file="/data/R02/huangtt39/F1_OSN/analysis/merge/RData/link_peak_anno_all.RData")
link_peak_to_xlsx<-merge(link_peak_anno_all,link_peaks_df_subset,by=c("gene","peak"))
write.xlsx(link_peak_to_xlsx,file="/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/link_peak_anno.xlsx")

link_peak_enhancer_count<-link_peak_anno_all%>%subset(peak_type=="enhancer")%>%group_by(gene)%>%summarise(count=n())
link_peak_enhancer_count<-merge(link_peak_enhancer_count, 
                gene_list %>% select(gene), 
                by = "gene", 
                all.y = TRUE) %>%  # 右连接
  replace(is.na(.), 0) 
link_peak_enhancer_count$gene_type<-"non-RME"
RME<-read.xlsx("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/Figure2h.xlsx")
for (i in 1:nrow(link_peak_enhancer_count)) { 
  gene <- link_peak_enhancer_count$gene[i]
  if(gene %in% RME$gene ){
    type="RME"
    link_peak_enhancer_count[i,]$gene_type=type
  }
} 

link_peak_enhancer_count_p<-ggviolin(link_peak_enhancer_count, x = "gene_type", y = "log2(count+1)", fill = "gene_type", 
              palette=c("#cdd2db","#45486f"),
              ylim=c(0, 7),
              alpha = 1,width = 0.5,
              legend = "none",#去掉legend
              xlab="non_RME vs RME", ylab="enhancer of each gene (log2)",
              font.tickslab = c(15,"plain","black"),
              add = "boxplot", 
              add.params = list(fill = "white", width = 0.1,linetype = 1)) +
              stat_compare_means(method = "wilcox.test", label = "p.format",comparisons = list(c("non-RME", "RME")),label.y = c(6))+
              theme_classic()+ 
              theme(plot.title = element_text(hjust = 0.5), 
                legend.position="right", 
                legend.title = element_blank(),
                text=element_text(size=8))
ggsave(filename = "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/link_peak_enhancer_count.pdf",link_peak_enhancer_count_p, device = "pdf", width = 8, height = 7, units = "cm", dpi = 300)
stat_results_enhancer <- link_peak_enhancer_count %>%
  rstatix::wilcox_test(
    count ~ gene_type,  # y变量 ~ x分组变量（和绘图一致）
    comparisons = list(c("non-RME", "RME")),  # 分组对比（注意和数据里的分组名完全一致）
    exact = FALSE,              # 和参考代码保持一致的参数
    conf.level = 0.95,
    detailed = TRUE             # 输出详细结果（包含置信区间）
  )

# ========== 2. 计算Wilcoxon效应量 ==========
effect_sizes_enhancer <- link_peak_enhancer_count %>%
  rstatix::wilcox_effsize(
    count ~ gene_type,
    comparisons = list(c("non-RME", "RME"))
  )

# ========== 3. 合并检验结果和效应量 ==========
final_results_enhancer <- stat_results_enhancer %>%
  # 选择需要的统计量列（和参考代码一致）
  select(group1, group2, n1, n2, statistic, p, conf.low, conf.high) %>%
  # 合并效应量
  left_join(
    effect_sizes_enhancer %>% select(group1, group2, effsize, magnitude), 
    by = c("group1", "group2")
  )

# ========== 4. 输出结果 ==========
print("增强子计数 (log2(count+1)) 统计结果：")
print(final_results_enhancer)

# 3. 批量处理所有基因
all_results <- list()
for (gene in high_RME$gene) {
  result <- process_gene_link_peak(gene, link_peak_anno_all, final_gene_dscore)
  if (!is.null(result)) {
    all_results[[gene]] <- result
  }
}

# 4. 合并并保存
final_results <- do.call(rbind, all_results)
write.csv(final_results, "/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/linked_peaks_counts.csv", row.names = FALSE)
final_results<-read.table("/data/R02/huangtt39/F1_OSN/analysis/figure/FigureS7/linked_peaks_counts.csv")

# length(unique(final_results$gene)) 16 genes called link peak
# 12 genes was called promoter
peak_link_subset<-final_results %>%subset((Maternal_specific_G1+Maternal_specific_G2)>=2&(Paternal_specific_G1+Paternal_specific_G2)>=2)
peak_link_subset$Maternal_specific_dscore<-peak_link_subset$Maternal_specific_G1/(peak_link_subset$Maternal_specific_G1+peak_link_subset$Maternal_specific_G2)-0.5
peak_link_subset$Paternal_specific_dscore<-peak_link_subset$Paternal_specific_G1/(peak_link_subset$Paternal_specific_G1+peak_link_subset$Paternal_specific_G2)-0.5

peak_link_subset$total_counts <- 
  peak_link_subset$Maternal_specific_G1 + 
  peak_link_subset$Maternal_specific_G2 +
  peak_link_subset$Paternal_specific_G1 + 
  peak_link_subset$Paternal_specific_G2

# 创建散点图
p <- ggplot(peak_link_subset %>% subset(peak_type != "promoter"), 
            aes(x = Paternal_specific_dscore, 
                y = Maternal_specific_dscore,
                size = log10(total_counts))) +
  geom_rect(aes(xmin = 0.15, xmax = 0.5, ymin = -0.15, ymax = -0.5),
            fill = "lightblue", alpha = 0.1, inherit.aes = FALSE) +
  
  # 添加散点
  geom_point() +
  
  # 分panel显示
  facet_wrap(~ peak_type, ncol = 2) +
  scale_x_continuous(limits = c(-0.5, 0.5)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  # 添加坐标轴交叉线（在0,0点）
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  
  # 坐标轴标签
  labs(
    x = "Paternal Group peak D-score",
    y = "Maternal Group peak D-score",
    title = "Peak Accessibility D-scores by Peak Type"
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
  "/data/R02/huangtt39/F1_OSN/new_analysis/peak_dscore_scatterplot.pdf",
  p,
  width = 9,
  height = 4,
  dpi = 300
)

peak_link_diff<-peak_link_subset%>%subset(Maternal_specific_dscore<(-0.1)&Paternal_specific_dscore>0.1&peak_type!="promoter")
#             gene                     peak Biallelic_G1 Biallelic_G2 Maternal_specific_G1 Maternal_specific_G2 Paternal_specific_G1 Paternal_specific_G2
# Chchd2.7  Chchd2 chr5-130256511-130257422           40           27                    2                    4                    7                    1
# Fetub.4    Fetub  chr16-22439296-22440221           59           63                    1                    2                   18                    9
# Tax1bp1  Tax1bp1   chr6-52639938-52640851           25           27                    2                    4                    6                    2
#               peak_type Maternal_specific_dscore Paternal_specific_dscore total_counts
# Chchd2.7 other_promoter               -0.1666667                0.3750000           14
# Fetub.4  other_promoter               -0.1666667                0.1666667           30
# Tax1bp1  other_promoter               -0.1666667                0.2500000           14
peak_link_diff_reserve<-peak_link_subset%>%subset(Maternal_specific_dscore>0.1&Paternal_specific_dscore<(-0.1)&peak_type!="promoter")
#           gene                     peak Biallelic_G1 Biallelic_G2 Maternal_specific_G1 Maternal_specific_G2 Paternal_specific_G1 Paternal_specific_G2
# Ago3.1    Ago3 chr4-126067199-126067866           23           26                    3                    1                    4                   10
# Bmp6.2    Bmp6  chr13-38366540-38367382           22           19                    2                    1                    1                    3
# Fetub.17 Fetub  chr16-23107007-23107843           26           26                    2                    1                    3                    6
# Syn2.4    Syn2 chr6-115127005-115127933           17           13                    2                    1                    3                    6
#          peak_type Maternal_specific_dscore Paternal_specific_dscore total_counts
# Ago3.1    enhancer                0.2500000               -0.2142857           18
# Bmp6.2    enhancer                0.1666667               -0.2500000            7
# Fetub.17  enhancer                0.1666667               -0.1666667           12
# Syn2.4    enhancer                0.1666667               -0.1666667           12
# function of promoter track

# 
calculate_other_gene_expression <- function(gene_name, other_promoter_gene, final_gene_dscore, dscore_data) {
  # 获取当前基因的细胞类型和barcode
  gene_cells <- subset(final_gene_dscore, gene == gene_name)
  if (nrow(gene_cells) == 0) return(NULL)
  
  # 创建结果数据框
  result <- data.frame(
    gene = gene_name,
    other_gene = other_promoter_gene,
    stringsAsFactors = FALSE
  )
  
  # 定义三种细胞类型
  cell_types <- c("Maternal_specific", "Paternal_specific", "Biallelic")
  
  # 对每种细胞类型
  for (cell_type in cell_types) {
    # 获取当前细胞类型的barcode
    type_barcodes <- gene_cells$barcode[gene_cells$merge_dscore_type == cell_type]
    
    if (length(type_barcodes) > 0) {
      # 从dscore_data中提取other_gene在这些barcode中的表达量
      other_gene_data <- subset(dscore_data, 
                                gene == other_promoter_gene & 
                                barcode %in% type_barcodes)
      
      if (nrow(other_gene_data) > 0) {
        # 计算G1总RNA
        g1_total <- sum(other_gene_data$G1rna, na.rm = TRUE)
        # 计算G2总RNA  
        g2_total <- sum(other_gene_data$G2rna, na.rm = TRUE)
        
        result[[paste0(cell_type, "_G1")]] <- g1_total
        result[[paste0(cell_type, "_G2")]] <- g2_total
      } else {
        # 如果没有数据，设为0
        result[[paste0(cell_type, "_G1")]] <- 0
        result[[paste0(cell_type, "_G2")]] <- 0
      }
    } else {
      # 如果没有这种细胞类型，设为0
      result[[paste0(cell_type, "_G1")]] <- 0
      result[[paste0(cell_type, "_G2")]] <- 0
    }
  }
  
  return(result)
}

other_promoter2caculate<-merge(peak_link_diff%>%subset(peak_type=="other_promoter")%>%select(gene,peak),link_peak_anno_all%>%select(peak,gene,other_gene),by=c("peak","gene"))
other_promoter_RNA <- list()
for (i in 1:nrow(other_promoter2caculate)) {
  result <- calculate_other_gene_expression(other_promoter2caculate[i,"gene"],other_promoter2caculate[i,'other_gene'], final_gene_dscore,dscore_data)
  if (!is.null(result)) {
    other_promoter_RNA[[i]] <- result
  }
}

# 4. 合并并保存
other_promoter_RNA <- do.call(rbind, other_promoter_RNA)
other_promoter_RNA$Maternal_specific_dscore<-other_promoter_RNA$Maternal_specific_G1/(other_promoter_RNA$Maternal_specific_G1+other_promoter_RNA$Maternal_specific_G2)-0.5
other_promoter_RNA$Paternal_specific_dscore<-other_promoter_RNA$Paternal_specific_G1/(other_promoter_RNA$Paternal_specific_G1+other_promoter_RNA$Paternal_specific_G2)-0.5

other_promoter_RNA$total_counts <- 
  other_promoter_RNA$Maternal_specific_G1 + 
  other_promoter_RNA$Maternal_specific_G2 +
  other_promoter_RNA$Paternal_specific_G1 + 
  other_promoter_RNA$Paternal_specific_G2
# r$> other_promoter_RNA
#      gene other_gene Maternal_specific_G1 Maternal_specific_G2 Paternal_specific_G1 Paternal_specific_G2 Biallelic_G1 Biallelic_G2 Maternal_specific_dscore
# 1   Fetub       Etv5                    2                    1                   17                   11           50           54                0.1666667
# 2  Chchd2       Sbds                    1                    0                    0                    0            5            7                0.5000000
# 3 Tax1bp1     Hibadh                    0                    5                    1                    2           16           13               -0.5000000
#   Paternal_specific_dscore total_counts
# 1                0.1071429           31
# 2                      NaN            1
# 3               -0.1666667            8
DefaultAssay(Obj)<-"ATAC"
Link_track <- function(gene_name) {
  # 获取基因信息和细胞数据
  gene_id <- id2gene %>% subset(gene == gene_name)
  gene_id <- gene_id[1, 1]
  promoter_bed <- max_promoter_data %>% 
    subset(gene == gene_name) %>% 
    select(chr, start, end)
  gene_bed <- bt.merge(subset(trans_bed, gene == gene_id))
  gene_bed <- data.frame(chr = c(gene_bed[1, 1]), 
                         start = min(gene_bed[, 2]), 
                         end = max(gene_bed[, 3]))
  
  # 确定原始基因的total_bed
  if (promoter_bed$start < gene_bed$start) {
    gene_total_bed <- str_c(gene_bed[1, 1], "-", promoter_bed$start, "-", promoter_bed$end + 10000)
  } else {
    gene_total_bed <- str_c(gene_bed[1, 1], "-", promoter_bed$start - 10000, "-", promoter_bed$end)
  }
  
  # 获取peak信息
  other_promoter_info <- other_promoter2caculate %>% 
    subset(gene == gene_name)
  
  if (nrow(other_promoter_info) == 0) {
    warning(paste("No peak information found for gene:", gene_name))
    # 如果没有peak信息，只画基因的track
    return(heat_track(gene_name))
  }
  
  peak_info <- other_promoter_info[1, ]
  other_gene_name <- peak_info$other_gene
  
  # 解析peak区域
  peak_parts <- strsplit(peak_info$peak, "-")[[1]]
  chr <- peak_parts[1]
  start <- as.numeric(peak_parts[2])
  end <- as.numeric(peak_parts[3])
  
  # 创建peak_ranges（使用原始peak范围，不包含5kb扩展）
  peak_ranges <- GRanges(
    seqnames = chr,
    ranges = IRanges(start = start, end = end)
  )
  
  # 扩增peak区域用于显示（扩增5kb）
  extended_start <- max(0, start - 5000)
  extended_end <- end + 5000
  peak_total_bed <- str_c(chr, "-", extended_start, "-", extended_end)
  
  # 获取细胞数据（只获取一次）
  heat_data <- final_gene_dscore %>% 
    subset(gene == gene_name, select = c(gene, barcode, mergescore, merge_dscore_type))
  rownames(heat_data) <- heat_data$barcode
  
  if (nrow(heat_data) == 0) {
    stop(paste("No cell data found for gene:", gene_name))
  }
  
  # 准备所有细胞的barcode
  all_barcodes <- c(str_c("G1_", heat_data$barcode), str_c("G2_", heat_data$barcode))
  Obj_subset <- subset(x = Obj, cells = all_barcodes)
  max_exp<-max(Obj_subset[["RNA"]]@data[gene_name,])
  max_exp<-ceiling(max_exp)
  max_exp_other<-max(Obj_subset[["RNA"]]@data[other_gene_name ,])
  max_exp_other<-ceiling(max_exp_other)
  # 创建存储plot对象的列表
  plot_objects <- list()
  
  # 处理三种细胞类型
  cell_types <- c("Paternal_specific", "Maternal_specific", "Biallelic")
  
  for (cell_type in cell_types) {
    barcode <- heat_data %>% subset(merge_dscore_type == cell_type)
    
    if (nrow(barcode) > 0) {
      # 准备细胞
      barcode_list <- barcode$barcode
      barcode_list <- c(str_c("G1_", barcode_list), str_c("G2_", barcode_list))
      temp_obj <- subset(x = Obj_subset, cells = barcode_list)
      
      # 添加类型标签
      label <- do.call(rbind, strsplit(WhichCells(temp_obj), "_"))[, 1]
      temp_obj@meta.data$type <- str_c(cell_type, "_", label)
      # 1. 基因区域的Coverage plot
      gene_cov <- CoveragePlot(
        object = temp_obj,
        region = gene_total_bed,
        group.by = "type",
        annotation = FALSE,
        tile = FALSE,
        peaks = FALSE
      ) & 
      scale_fill_manual(values = c("#505f84", "#b45a55")) &
      theme(
        axis.title = element_blank(),  
        axis.text = element_blank(),  
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank()
      )
      
      # 2. 基因区域的Expression plot
      if (gene_name %in% rownames(temp_obj[["RNA"]]@data)) {
        gene_exp <- ExpressionPlot(temp_obj, features = gene_name, group.by = "type", assay = "RNA") + 
          scale_fill_manual(values = c("#505f84", "#b45a55")) +
          theme(text = element_text(size = 8))+scale_x_continuous(limits = c(0, max_exp))
      } else {
        gene_exp <- ggplot() + theme_void()
      }
      
      # 3. Peak区域的Coverage plot
      peak_cov <- CoveragePlot(
        object = temp_obj,
        region = peak_total_bed,
        group.by = "type",
        annotation = FALSE,
        tile = FALSE,
        peaks = FALSE,
        region.highlight = peak_ranges
      ) & 
      scale_fill_manual(values = c("#505f84", "#b45a55")) &
      theme(
        axis.title = element_blank(),  
        axis.text = element_blank(),  
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        plot.title = element_blank()
      )
      
      # 4. Other gene的Expression plot
      if (other_gene_name %in% rownames(temp_obj[["RNA"]]@data)) {
        other_exp <- ExpressionPlot(temp_obj, features = other_gene_name, group.by = "type", assay = "RNA") + 
          scale_fill_manual(values = c("#505f84", "#b45a55")) +
          theme(text = element_text(size = 8))+scale_x_continuous(limits = c(0, max_exp_other))
      } else {
        other_exp <- ggplot() + theme_void()
      }
      
      # 添加到plot列表
      plot_objects[[paste0(cell_type, "_gene_cov")]] <- gene_cov
      plot_objects[[paste0(cell_type, "_gene_exp")]] <- gene_exp
      plot_objects[[paste0(cell_type, "_peak_cov")]] <- peak_cov
      plot_objects[[paste0(cell_type, "_other_exp")]] <- other_exp
    }
  }
  
  # 5. Total coverage plot - 基因区域
  gene_total_cov <- CoveragePlot(
    object = Combined_mOSN,
    region = gene_total_bed,
    annotation = FALSE,
    tile = FALSE,
    peaks = FALSE
  ) & 
  scale_fill_manual(values = c("#c5cd94")) &
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),  
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
  
  # 6. Total coverage plot - peak区域
  peak_total_cov <- CoveragePlot(
    object = Combined_mOSN,
    region = peak_total_bed,
    annotation = FALSE,
    tile = FALSE,
    peaks = FALSE,
    region.highlight = peak_ranges
  ) & 
  scale_fill_manual(values = c("#c5cd94")) &
  theme(
    axis.title = element_blank(),  
    axis.text = element_blank(),  
    axis.line = element_blank(),
    axis.ticks = element_blank()
  )
  
  # 7. Annotation plots
  gene_ap <- AnnotationPlot(object = Obj, region = gene_total_bed) + 
    theme(
      text = element_text(size = 8),
      axis.title.y = element_blank(),  
      axis.text.y = element_blank(),   
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  peak_ap <- AnnotationPlot(object = Obj, region = peak_total_bed) + 
    theme(
      text = element_text(size = 8),
      axis.title.y = element_blank(),  
      axis.text.y = element_blank(),   
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank()
    )
  
  # 组合所有图形 - 使用类似于原始heat_track的布局方式
  # 左列：基因区域
  left_column <- gene_total_cov +plot_spacer()+
      plot_objects[["Paternal_specific_gene_cov"]] +  plot_objects[["Paternal_specific_gene_exp"]] +
    plot_objects[["Maternal_specific_gene_cov"]]  +plot_objects[["Maternal_specific_gene_exp"]] + 
    plot_objects[["Biallelic_gene_cov"]]+plot_objects[["Biallelic_gene_exp"]] +
    gene_ap+plot_spacer()
  
  # 右列：peak区域
  right_column <- peak_total_cov +plot_spacer()+ 
    plot_objects[["Paternal_specific_peak_cov"]] + plot_objects[["Paternal_specific_other_exp"]] +
     plot_objects[["Maternal_specific_peak_cov"]]  +plot_objects[["Maternal_specific_other_exp"]] +
    plot_objects[["Biallelic_peak_cov"]]+plot_objects[["Biallelic_other_exp"]] +
    peak_ap +plot_spacer()
  # 横向组合两列
  return(list(
    left_plot = left_column,
    right_plot = right_column
  ))
}
for(i in other_promoter2caculate$gene){
  file <- str_c("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure2/RME/link/same/", i, ".pdf")
  
  # 获取两个图
  plots <- Link_track(i)
  
  # 创建PDF文件
  pdf(file, width = 7, height = 8, bg = "white")
  
  # 第一页：基因区域
  print(plots$left_plot  +plot_layout(nrow=5,ncol=2,heights = c(1.2,1.2,1.2,1.2,0.6),widths=c(4,1)))
  
  # 第二页：peak区域
  print(plots$right_plot +plot_layout(nrow=5,ncol=2,heights = c(1.2,1.2,1.2,1.2,0.6),widths=c(4,1)))
  
  dev.off()
  
  cat("Saved", i, "to", file, "\n")
}


heat_track_for_chip<-function(gene_name){
  gene_id<-id2gene%>%subset(gene==gene_name)
  gene_id<-gene_id[1,1]
  promoter_bed<-max_promoter_data%>%subset(gene==gene_name)%>%select(chr,start,end)
  gene_bed<-bt.merge(subset(trans_bed,gene==gene_id))
  gene_bed<-data.frame(chr=c(gene_bed[1,1]),start=min(gene_bed[,2]),end=max(gene_bed[,3]))
  if(promoter_bed$start<gene_bed$start){
    total_bed<-str_c(gene_bed[1,1],"-",promoter_bed$start,"-",gene_bed[1,3])
  }else{
    total_bed<-str_c(gene_bed[1,1],"-",gene_bed[1,2],"-",promoter_bed$end)
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

# Add chip-track

# track plot of histone modification
seqlength_file<-read.table(file="/data/R02/huangtt39/F1_OSN/mapping/reference/mm10-mask/fasta/genome.fa.fai")
seqlength<-seqlength_file$V2
names(seqlength)<-seqlength_file$V1
chip_dir<-"/data/R02/huangtt39/data/OSN/chip/"
chip_type<-c("H3K9me3","H3K4me3")

# function of bigwig track
cell_BigwigTrack <- function(
  region,
  bigwig,
  smooth = 100,
  y_label = "bigWig",
  bigwig.scale = "common",
  ymax = NULL,
  max.downsample = 3000,
  downsample.rate = 0.1,
  type="coverage"
) {
  if (!inherits(x = bigwig, what = "list")) {
    bigwig <- list("bigWig" = bigwig)
  }
  if (.Platform$OS.type == "windows") {
    message("BigwigTrack not supported on Windows")
    return(NULL)
  }
  if (!requireNamespace("rtracklayer", quietly = TRUE)) {
    message("Please install rtracklayer. http://www.bioconductor.org/packages/rtracklayer/")
    return(NULL)
  }
  if (!inherits(x = region, what = "GRanges")) {
    stop("region should be a GRanges object")
  }
  all.data <- data.frame()
  for (i in seq_along(bigwig)) {
    region_data <- rtracklayer::import(
      con = bigwig[[i]],
      which = region,
      as = "NumericList"
    )[[1]]
    if (!is.null(x = smooth)) {
      region_data <- roll_mean(x = region_data, n = smooth, fill = 0L)
    }
    region_data <- data.frame(
      position = start(x = region):end(x = region),
      score = region_data,
      stringsAsFactors = FALSE,
      bw = names(x = bigwig)[[i]]
    )
    if (bigwig.scale == "separate") {
      # scale to fraction of max for each separately
      file.max <- max(region_data$score, na.rm = TRUE)
      region_data$score <- region_data$score / file.max
    }
    all.data <- rbind(all.data, region_data)
  }
  all.data$bw <- factor(x = all.data$bw, levels = names(x = bigwig))
  window.size = width(x = region)
  sampling <- ceiling(x = max(max.downsample, window.size * downsample.rate))
  coverages <- slice_sample(.data = all.data, n = sampling)
  
  covmax <- signif(x = max(coverages$score, na.rm = TRUE), digits = 0)
  if (is.null(x = ymax)) {
    ymax <- covmax
  } else if (is.character(x = ymax)) {
    if (!startsWith(x = ymax, prefix = "q")) {
      stop("Unknown ymax requested. Must be NULL, a numeric value, or 
           a quantile denoted by 'qXX' with XX the desired quantile value,
           e.g. q95 for 95th percentile")
    }
    percentile.use <- as.numeric(
      x = sub(pattern = "q", replacement = "", x = as.character(x = ymax))
    ) / 100
    ymax <- covmax * percentile.use
  }
  ymin <- 0
  
  # perform clipping
  coverages$score[coverages$score > ymax] <- ymax 
  if (type == "coverage") {
    p <- ggplot(
      data = coverages,
      mapping = aes_string(x = "position", y = "score", fill = "bw")
    ) + geom_area() +
      facet_wrap(facets = ~bw, strip.position = "left", ncol = 1)
  }
  chromosome <- as.character(x = seqnames(x = region))
  # p <- p + Signac::theme_browser(axis.text.y = TRUE) +
  #   xlab(label = paste0(chromosome, " position (bp)")) +
  #   ylab(label = y_label)
  p <- p + Signac::theme_browser(axis.text.y = TRUE)+scale_fill_manual(values=c())
  return(p)
}

setbigwig<-function(gene_name){
  gene_id<-id2gene%>%subset(gene==gene_name)
  gene_id<-gene_id[1,1]
  promoter_bed<-max_promoter_data%>%subset(gene==gene_name)%>%select(chr,start,end)
  gene_bed<-bt.merge(subset(trans_bed,gene==gene_id))
  gene_bed<-data.frame(chr=c(gene_bed[1,1]),start=min(gene_bed[,2]),end=max(gene_bed[,3]))

  if(promoter_bed$start<gene_bed$start){
    total_bed<-c(gene_bed[1,1],as.numeric(promoter_bed$start),gene_bed[1,3])
  }else{
    total_bed<-c(gene_bed[1,1],gene_bed[1,2],as.numeric(promoter_bed$end))
  }
  bwlist<-list("/data/R02/huangtt39/data/OSN/chip/H3K9me3/GSM2454912_mOSN_H3K9me3_Mnase-ChIP_rep1.bigWig","/data/R02/huangtt39/data/OSN/chip/H3K9me3/GSM2454913_mOSN_H3K9me3_Mnase-ChIP_rep2.bigWig")
  bwlist1<-list("/data/R02/huangtt39/data/OSN/chip/H3K4me3/GSM7947099_OSN_WT_H3K4me3_Rep1.bw","/data/R02/huangtt39/data/OSN/chip/H3K4me3/GSM7947100_OSN_WT_H3K4me3_Rep2.bw")
  names(bwlist)<-c("mOSN_H3K9me3_rep1","mOSN_H3K9me3_rep2")
  names(bwlist1)<-c("mOSN_H3K4me3_rep1","mOSN_H3K4me3_rep2")
  bwp<-cell_BigwigTrack(GRanges(seqnames= total_bed[1],IRanges(start= as.numeric(total_bed[2])+1,end= as.numeric(total_bed[3]))),bwlist,smooth=NULL)+scale_fill_manual(values=c("#996598","#996598"))+NoLegend()
  bwp1<-cell_BigwigTrack(GRanges(seqnames= total_bed[1],IRanges(start= as.numeric(total_bed[2])+1,end= as.numeric(total_bed[3]))),bwlist1,smooth=NULL)+scale_fill_manual(values=c("#9b80b3","#9b80b3"))+NoLegend()
  return(bwp+bwp1)
}

# the promoter track plot
for(i in c("Ago3","Wdr17")){
  file=str_c("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure/",i,".pdf")
  ggsave(filename = file,heat_track_for_chip(i)+reads_bar_gene(i)+plot_layout(nrow=8,ncol=2,heights = c(1,1.2,1.2,1.2,0.4,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}
# the histone modification plot
for(i in c("Ago3","Wdr17")){
  file=str_c("/data/R02/huangtt39/F1_OSN/analysis/figure/Figure/",i,"_chip.pdf")
  ggsave(filename = file,setbigwig(i)+plot_layout(ncol=1), device = "pdf", width = 20, height =10, units = "cm", dpi = 300,bg="white")
}