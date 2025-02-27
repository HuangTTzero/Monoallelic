# the histone modification of RME
# load result of chip-seq calling
chip_dir<-"/data/R02/huangtt39/data/OSN/chip/"
chip_type<-c("H3K9me3","H3K27ac")
chip_peak<-lapply(chip_type,function(type){
    temp<-lapply(1:2,function(rep){
      file<-str_c(chip_dir,type,"/callpeak/","mOSN_",type,"_Mnase-ChIP_rep",rep,".bed")
      ori_bed<-read.table(file,sep = '\t',header = FALSE,quote = "")%>%select(V1,V2,V3,V4)
      colnames(ori_bed)<-c("chr","start","end","id")
      ori_bed
    })
    temp
})
# merge the peak of duplicates
chip_peak_merge<-lapply(1:2,function(i){
    temp<-bt.merge(bt.sort(rbind(chip_peak[[i]][[1]],chip_peak[[i]][[2]])))
    temp$id<-paste(chip_type[i],1:nrow(temp),sep="_")
    temp
})
# get the peak conserved in duplicates
chip_peak_final<-lapply(1:2,function(i){
    a<-bt.intersect(a=chip_peak_merge[[i]],b=chip_peak[[i]][[1]],wa=TRUE)
    b<-bt.intersect(a=chip_peak_merge[[i]],b=chip_peak[[i]][[2]],wa=TRUE)
    both<-unique(c(a$V4,b$V4)[duplicated(c(a$V4,b$V4))])
    chip_peak_merge[[i]]%>%subset(id%in%both)
})
# extend the promoter region to 2000 bp, allow incomplete overlap of peaks 
promoter_bed_2000<-tss_bed%>%mutate(start=case_when(strand=="+" ~ end-2000-1,
                         strand=="-" ~end-2000),end=case_when(strand=="+" ~ end+2000-1,
                         strand=="-" ~end+2000))
promoter_bed_2000$start[promoter_bed_2000$start<0]<-0
# overlap with extended promoter 
H3K9me3_promoter<-bt.intersect(a=promoter_bed_2000,b=chip_peak_final[[1]],wa=TRUE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",geneid="V4",score="V5",strand="V6"))
H3K27ac_promoter<-bt.intersect(a=promoter_bed_2000,b=chip_peak_final[[2]],wa=TRUE)%>%dplyr::rename(c(chr="V1",start="V2",end="V3",geneid="V4",score="V5",strand="V6"))
# the promoter both with H3K27ac and H3K9me3
chip_both<-merge(H3K27ac_promoter,H3K9me3_promoter,by=c("chr","start","end","geneid","score","strand"))
chip_both<-merge(chip_both,id2gene,by="geneid")

merge(chip_both%>%dplyr::select(gene,geneid)%>%unique(),candicate_RME,by="gene")
#    gene             geneid   type count Biallelic Maternal_specific Paternal_specific Parental_specific PM_type
# 1 Clip4 ENSMUSG00000024059 Others   521 0.7600768        0.22072937        0.01919386         0.2399232  BI_RME
# 2 Fetub ENSMUSG00000022871 Others   275 0.7781818        0.02909091        0.19272727         0.2218182  BI_RME
# 3 Hsbp1 ENSMUSG00000031839 Others    15 0.6666667        0.20000000        0.13333333         0.3333333  BI_RME
chip_gene<-c("Clip4","Fetub","Hsbp1")

# track plot of histone modification
seqlength_file<-read.table(file="/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-mask/fasta/genome.fa.fai")
seqlength<-seqlength_file$V2
names(seqlength)<-seqlength_file$V1
chip_dir<-"/data/R02/huangtt39/data/OSN/chip/"
chip_type<-c("H3K9me3","H3K27ac")
# get the bigwig of Clip4, Fetub, Hsbp1
chip_bedgraph<-lapply(chip_type,function(type){
  repeat_bed<-lapply(c(1,2),function(i){
    filename<-str_c(chip_dir,type,"/mOSN_",type,"-ChIP_rep",i,".bedGraph")
    ori_bed<-read.table(filename,col.names = c("chr","start","end","score"),sep = '\t',header = FALSE,quote = "")
    ori_bed<-ori_bed%>%subset(chr%in%str_c("chr",c(1:19,"X","Y")))
    ori_bed<-GRanges(seqnames= ori_bed$chr,ranges=IRanges(start= as.numeric(ori_bed$start)+1,end= as.numeric(ori_bed$end)),score=ori_bed$score)
    seqlengths(ori_bed)<- seqlength[names(seqlengths(ori_bed))]
    export.bw(ori_bed, str_c(chip_dir,type,"/mOSN_",type,"-ChIP_rep",i,"_seqlength.bw"))
  })
  return()
})
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
    total_bed<-c(gene_bed[1,1],as.numeric(promoter_bed$start),as.numeric(promoter_bed$end)+10000)
  }else{
    total_bed<-c(gene_bed[1,1],as.numeric(promoter_bed$start)-10000,as.numeric(promoter_bed$end))
  }
  bwlist<-list("/data/R02/huangtt39/data/OSN/chip/H3K9me3/mOSN_H3K9me3-ChIP_rep1_seqlength.bw","/data/R02/huangtt39/data/OSN/chip/H3K9me3/mOSN_H3K9me3-ChIP_rep2_seqlength.bw")
  bwlist1<-list("/data/R02/huangtt39/data/OSN/chip/H3K27ac/mOSN_H3K27ac-ChIP_rep1_seqlength.bw","/data/R02/huangtt39/data/OSN/chip/H3K27ac/mOSN_H3K27ac-ChIP_rep2_seqlength.bw")
  names(bwlist)<-c("mOSN_H3K9me3_rep1","mOSN_H3K9me3_rep2")
  names(bwlist1)<-c("mOSN_H3K27ac_rep1","mOSN_H3K27ac_rep2")
  bwp<-cell_BigwigTrack(GRanges(seqnames= total_bed[1],IRanges(start= as.numeric(total_bed[2])+1,end= as.numeric(total_bed[3]))),bwlist,smooth=NULL)+scale_fill_manual(values=c("#996598","#996598"))+NoLegend()
  bwp1<-cell_BigwigTrack(GRanges(seqnames= total_bed[1],IRanges(start= as.numeric(total_bed[2])+1,end= as.numeric(total_bed[3]))),bwlist1,smooth=NULL)+scale_fill_manual(values=c("#9b80b3","#9b80b3"))+NoLegend()
  return(bwp+bwp1)
}

# Figure 3A-D
# the promoter track plot
for(i in c("Clip4","Fetub","Hsbp1")){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/chip/",i,".pdf")
  ggsave(filename = file,heat_track(i)+reads_bar_gene(i)+plot_layout(nrow=8,ncol=2,heights = c(1,1.2,1.2,1.2,0.4,3),widths=c(4,1)), device = "pdf", width = 20, height =20, units = "cm", dpi = 300,bg="white")
}
# the histone modification plot
for(i in chip_gene){
  file=str_c("/data/R02/huangtt39/ATAC-RNAseq/analysis/figure/Figure2/RME/chip/",i,"_chip.pdf")
  ggsave(filename = file,setbigwig(i)+plot_layout(ncol=1), device = "pdf", width = 20, height =10, units = "cm", dpi = 300,bg="white")
}
