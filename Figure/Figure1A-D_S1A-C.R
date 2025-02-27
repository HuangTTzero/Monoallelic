## Load required packages
library(Signac)
library(Seurat)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(MASS)
library(viridis)
library(dplyr)
library(sctransform)
library(SummarizedExperiment)
#library(glmGamPoi)
library(stringr)
library(RColorBrewer)
library(cowplot)
library(tidydr)
set.seed(100)


# get gene annotations for mm10, build from the cellranger gtf reference, a newer version than EnsDb.Mmusculus.v79 
annotation<-import('/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/gencode.vM23.primary_assembly.annotation.gtf')
annotation<-annotation[,c("type","transcript_id","gene_name","gene_id","gene_type")][annotation$type%in%c("exon","CDS","UTR"),]
annotation$tx_id<-annotation$transcript_id
annotation$gene_biotype<-annotation$gene_type
annotation$type<-as.character(annotation$type)
annotation$type<-replace(annotation$type,annotation$type=="CDS","cds")
annotation$type<-replace(annotation$type,annotation$type=="UTR","utr")
annotation$type<-factor(annotation$type,levels = c("utr","cds","exon"))
annotation<-annotation[,c("type","tx_id","gene_name","gene_id","gene_biotype")]


## Load the data
dirs <- c("/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_mask/joint202201_mask/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_mask/jointF1OSN_mask/outs/","/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_mask/jointF18w_mask/outs/")
samples <- c("dataset1","dataset2","dataset3")
objList <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_feature_bc_matrix.h5")) 
  fragpath <- str_c(dirs[i],"atac_fragments.tsv.gz")
  metadata <- read.csv(str_c(dirs[i],"per_barcode_metrics.csv"),row.names=1) %>% dplyr::filter(is_cell==1) 
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
  obj
})

## Quality control
for ( i in 1:length(objList)){
  DefaultAssay(objList[[i]]) <- "RNA"
  objList[[i]] <- PercentageFeatureSet(objList[[i]], pattern = "^mt-", col.name = "percent.mt")
  DefaultAssay(objList[[i]]) <- "ATAC"
  # Computes the ratio of fragments between 147 bp and 294 bp (mononucleosome) to fragments < 147 bp (nucleosome-free)
  objList[[i]] <- NucleosomeSignal(objList[[i]]) 
  # Compute the transcription start site (TSS) enrichment score for each cell
  objList[[i]] <- TSSEnrichment(objList[[i]],fast = FALSE)
  # add fraction of reads in peaks
  objList[[i]]$pct_reads_in_peaks <- objList[[i]]$atac_peak_region_fragments / objList[[i]]$atac_fragments * 100
  objList[[i]]$nucleosome_group <- ifelse(objList[[i]]$nucleosome_signal > 1, 'NS > 1', 'NS <= 1')
}

# qc plot, Figure S1A-B

merge_obj<-merge(objList[[1]],objList[[2]])
merge_obj<-merge(merge_obj,objList[[3]])

qc_data1=VlnPlot(
  object =merge_obj,
  features = c('nFeature_RNA',
               'percent.mt',
               'TSS.enrichment', 'nucleosome_signal','pct_reads_in_peaks'),
  cols=c("#517e8e", "#517e8e", "#517e8e"),
  split.by =  "orig.ident",
  pt.size =0,combine = FALSE)
hline<-list()
hline[[1]]<-c(1500,3000,4500,100,200,500)
hline[[2]]<-c(15,10,5)
hline[[3]]<-c(5)
hline[[4]]<-c(1)
hline[[5]]<-c(40)
for ( i in 1:5){
  qc_data1[[i]]<-qc_data1[[i]]+geom_boxplot()+NoLegend()+geom_hline(yintercept= hline[[i]],linetype="dashed")+labs(x=NULL,y=NULL)
}
vln_all<-plot_grid(plotlist = qc_data1, align = "h", nrow = 3)
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/qc_all.pdf", vln_all, device = "pdf", width = 30, height = 30, units = "cm", dpi = 300,bg="white")

qc_data2=VlnPlot(
  object =merge_obj,
  features = c('nCount_RNA','atac_peak_region_fragments'),
  cols=c("#517e8e", "#517e8e", "#517e8e"),
  split.by =  "orig.ident",
  pt.size =0,combine = FALSE,log=TRUE)
hline1<-list()
hline1[[1]]<-c(2000,8000,20000)
hline1[[2]]<-c(30000,45000,60000,2000,3000)
for ( i in 1:2){
  qc_data2[[i]]<-qc_data2[[i]]+geom_boxplot()+NoLegend()+geom_hline(yintercept=hline1[[i]],linetype="dashed")+labs(x=NULL,y=NULL)
}
vln_all1<-plot_grid(plotlist = qc_data2, align = "h", nrow = 1)
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/qc_all1.pdf", vln_all1, device = "pdf", width = 30, height = 10, units = "cm", dpi = 300,bg="white")


##qc
objList[[1]] <- subset(
  x = objList[[1]],
  subset = nCount_RNA <= 2000 &
  nFeature_RNA >= 100 &
  nFeature_RNA <= 1500 & 
  percent.mt < 15 &
  nCount_ATAC >= 2000 &
  nCount_ATAC <= 30000 &
  TSS.enrichment >= 5 &
  nucleosome_signal < 1 &
  pct_reads_in_peaks >= 40
  )

objList[[2]] <- subset(
  x = objList[[2]],
  subset = nCount_RNA <= 8000 &
  nFeature_RNA >= 200 &
  nFeature_RNA <= 3000 & 
  percent.mt < 10 &
  nCount_ATAC >= 3000 &
  nCount_ATAC <= 45000 &
  TSS.enrichment >= 5 &
  nucleosome_signal < 1 &
  pct_reads_in_peaks >= 40
  )
objList[[3]] <- subset(
  x = objList[[3]],
  subset = nCount_RNA <= 20000 &
  nFeature_RNA >= 500 &
  nFeature_RNA <= 4500 & 
  percent.mt < 5 &
  nCount_ATAC >= 3000 &
  nCount_ATAC <= 60000 &
  TSS.enrichment >= 5 &
  nucleosome_signal < 1 &
  pct_reads_in_peaks >= 40
  )

## Find doublet
## doublet detection in single-cell ATAC-seq
# reference : https://www.bioconductor.org/packages/devel/bioc/vignettes/scDblFinder/inst/doc/scATAC.html
library(scDblFinder)
dataset1.sce <- as.SingleCellExperiment(objList[[1]],assay="ATAC")
dataset2.sce <- as.SingleCellExperiment(objList[[2]],assay="ATAC")
dataset3.sce <- as.SingleCellExperiment(objList[[3]],assay="ATAC")

# Applying the scDblFinder method
dataset1.sce <- scDblFinder(dataset1.sce, aggregateFeatures=TRUE, processing="normFeatures")
#1126 (8.9%) doublets called
dataset2.sce <- scDblFinder(dataset2.sce, aggregateFeatures=TRUE, processing="normFeatures")
#1027 (7.5%) doublets called
dataset3.sce <- scDblFinder(dataset3.sce, aggregateFeatures=TRUE, processing="normFeatures")
# 867(7.3%) doublets called

# Using the Amulet method
# assumption : cells with loci covered by more than two fragments are indicative of the droplet being a doublet
dataset1_fragfile <- "/data/R02/huangtt39/ATAC-RNAseq/mapping/joint202201/mapping_mask/joint202201_mask/outs/atac_fragments.tsv.gz"
dataset2_fragfile <- "/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF1OSN/mapping_mask/jointF1OSN_mask/outs/atac_fragments.tsv.gz"
dataset3_fragfile <- "/data/R02/huangtt39/ATAC-RNAseq/mapping/jointF18w/mapping_mask/jointF18w_mask/outs/atac_fragments.tsv.gz"

# excluding the mitochondrial and sex chromosomes, as well as repetitive regions
mm10_blacklist_bed <- read.table("/data/R02/huangtt39/data/mm10/annotation/mm10-blacklist.v2.bed",sep="\t")
mm10_blacklist.gr <- GRanges(mm10_blacklist_bed$V1,
  IRanges(start = mm10_blacklist_bed$V2+1, end = mm10_blacklist_bed$V3))

mm10.chrom.sizes <- read.table("/data/R02/huangtt39/data/mm10/annotation/mm10.chrom.sizes")
excludedChroms.sizes_01 <- mm10.chrom.sizes[grep("random",mm10.chrom.sizes$V1),]
excludedChroms.sizes_01$V1 <- gsub("chr\\w\\_","",excludedChroms.sizes_01$V1)
excludedChroms.sizes_01$V1 <- gsub("_random","",excludedChroms.sizes_01$V1)
excludedChroms.sizes_01$V1 <- str_c(excludedChroms.sizes_01$V1,".1")

excludedChroms.sizes_02 <- mm10.chrom.sizes[grep("chrUn",mm10.chrom.sizes$V1),]
excludedChroms.sizes_02$V1 <- gsub("chrUn_","",excludedChroms.sizes_02$V1)
excludedChroms.sizes_02$V1 <- str_c(excludedChroms.sizes_02$V1,".1")

excludedChroms.sizes_03 <- mm10.chrom.sizes[which(mm10.chrom.sizes$V1 %in% str_c("chr",c("X","Y","M"))),]

excludedChroms.gr <- GRanges(c(excludedChroms.sizes_01$V1,excludedChroms.sizes_02$V1,excludedChroms.sizes_03$V1),
  IRanges(start = 1, end = c(excludedChroms.sizes_01$V2,excludedChroms.sizes_02$V2,excludedChroms.sizes_03$V2))) 

toExclude.gr <- c(mm10_blacklist.gr, excludedChroms.gr)

# launch the Amulet method
library(BiocParallel)
multicore <- MulticoreParam(workers = 8)
dataset1_res <- amulet(dataset1_fragfile, regionsToExclude=toExclude.gr,barcodes=colnames(objList[[1]]),BPPARAM=multicore)
dataset2_res <- amulet(dataset2_fragfile, regionsToExclude=toExclude.gr,barcodes=colnames(objList[[2]]),BPPARAM=multicore)
dataset3_res <- amulet(dataset3_fragfile, regionsToExclude=toExclude.gr,barcodes=colnames(objList[[3]]),BPPARAM=multicore)

## Combining mehtods
# The Amulet method tends to perform best with datasets that have homotypic doublets and where cells have a high library size (i.e. median library size per cell of 10-15k reads),while the scDblFinder-based approach works better for heterotypic doublets

# The Clamulet method (Classification-powered Amulet-like method)
# It generates artificial doublets by operating on the fragment coverages. This has the advantage that the number of loci covered by more than two reads can be computed for artificial doublets, enabling the use of this feature (along with the kNN-based ones) in a classification scheme. 
# multicore <- MulticoreParam(workers = 8)
# dataset1_clamulet_res <- clamulet(dataset1_fragfile,regionsToExclude=toExclude.gr,barcodes=colnames(objList[[1]]),BPPARAM=multicore)

# Simple p-value combination
dataset1_res$scDblFinder.p <- 1-colData(dataset1.sce)[row.names(dataset1_res), "scDblFinder.score"]
dataset1_res$combined <- apply(dataset1_res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
  x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
  suppressWarnings(aggregation::fisher(x))
})
dataset1_res$doublet <- ifelse(dataset1_res$combined<0.05,"Yes","No")
write.table(dataset1_res,"/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/dataset1_scDblFinder_result.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
objList[[1]] <- AddMetaData(objList[[1]],dataset1_res[,"doublet",drop=FALSE])
objList[[1]]$doublet <- factor(objList[[1]]$doublet,levels=c("Yes","No"),labels=c("Doublet","Singlet"))
#1386/12675

dataset2_res$scDblFinder.p <- 1-colData(dataset2.sce)[row.names(dataset2_res), "scDblFinder.score"]
dataset2_res$combined <- apply(dataset2_res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
  x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
  suppressWarnings(aggregation::fisher(x))
})
dataset2_res$doublet <- ifelse(dataset2_res$combined<0.05,"Yes","No")
write.table(dataset2_res,"/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/dataset2_scDblFinder_result.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
objList[[2]] <- AddMetaData(objList[[2]],dataset2_res[,"doublet",drop=FALSE])
objList[[2]]$doublet <- factor(objList[[2]]$doublet,levels=c("Yes","No"),labels=c("Doublet","Singlet"))
#1712/13730

dataset3_res$scDblFinder.p <- 1-colData(dataset3.sce)[row.names(dataset3_res), "scDblFinder.score"]
dataset3_res$combined <- apply(dataset3_res[,c("scDblFinder.p", "p.value")], 1, FUN=function(x){
  x[x<0.001] <- 0.001 # prevent too much skew from very small or 0 p-values
  suppressWarnings(aggregation::fisher(x))
})
dataset3_res$doublet <- ifelse(dataset3_res$combined<0.05,"Yes","No")
write.table(dataset3_res,"/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/dataset3_scDblFinder_result.txt",row.names=TRUE,col.names=TRUE,sep="\t",quote=FALSE)
objList[[3]] <- AddMetaData(objList[[3]],dataset3_res[,"doublet",drop=FALSE])
objList[[3]]$doublet <- factor(objList[[3]]$doublet,levels=c("Yes","No"),labels=c("Doublet","Singlet"))
#1769/11831

## Doublet umap
for ( i in 1:length(objList)){
  DefaultAssay(objList[[i]]) <- "RNA"
  objList[[i]] <- NormalizeData(objList[[i]])
  objList[[i]] <- FindVariableFeatures(objList[[i]], selection.method = "vst", nfeatures = 2000)
  objList[[i]] <- ScaleData(objList[[i]])
  objList[[i]] <- RunPCA(objList[[i]])
#plot_pca <- ElbowPlot(osn_filter, ndims=50, reduction="pca") 
  objList[[i]] <- RunUMAP(objList[[i]],dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
  DefaultAssay(objList[[i]]) <- "ATAC"
  objList[[i]] <- FindTopFeatures(objList[[i]], min.cutoff = "q0")
  objList[[i]] <- RunTFIDF(objList[[i]])
  objList[[i]] <- RunSVD(objList[[i]])
  objList[[i]] <- RunUMAP(objList[[i]], reduction = 'lsi', dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
  DefaultAssay(objList[[i]]) <- "RNA"
  #聚类和非线性降维
  objList[[i]] <- FindMultiModalNeighbors(objList[[i]], reduction.list = list("pca", "lsi"), dims.list = list(1:30, 2:20))
  objList[[i]] <- FindClusters(objList[[i]], graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution = 0.8)
  objList[[i]] <- RunUMAP(objList[[i]], nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
}

data1_rna_umap <- DimPlot(objList[[1]], reduction = "umap.rna",group.by="doublet", repel = TRUE, pt.size = 0.5) + ggtitle("RNA")+NoLegend()
data1_atac_umap <- DimPlot(objList[[1]], reduction = "umap.atac",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("ATAC")+NoLegend()
data1_wnn_umap <- DimPlot(objList[[1]], reduction = "wnn.umap",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("WNN")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/data1_doublet.pdf",data1_rna_umap +data1_atac_umap + data1_wnn_umap & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 24, height = 8, units = "cm", dpi = 300)


data2_rna_umap <- DimPlot(objList[[2]], reduction = "umap.rna",group.by="doublet", repel = TRUE, pt.size = 0.5) + ggtitle("RNA")+NoLegend()
data2_atac_umap <- DimPlot(objList[[2]], reduction = "umap.atac",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("ATAC")+NoLegend()
data2_wnn_umap <- DimPlot(objList[[2]], reduction = "wnn.umap",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("WNN")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/data2_doublet.pdf",data2_rna_umap +data2_atac_umap + data2_wnn_umap & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 24, height = 8, units = "cm", dpi = 300)


data3_rna_umap <- DimPlot(objList[[3]], reduction = "umap.rna",group.by="doublet", repel = TRUE, pt.size = 0.5) + ggtitle("RNA")+NoLegend()
data3_atac_umap <- DimPlot(objList[[3]], reduction = "umap.atac",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("ATAC")+NoLegend()
data3_wnn_umap <- DimPlot(objList[[3]], reduction = "wnn.umap",group.by="doublet",repel = TRUE, pt.size = 0.5) + ggtitle("WNN")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/data3_doublet.pdf",data3_rna_umap +data3_atac_umap + data3_wnn_umap & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 24, height = 8, units = "cm", dpi = 300)

##remove doublet
data1_cell<-WhichCells(subset(x=objList[[1]],subset=doublet=="Singlet"))
data2_cell<-WhichCells(subset(x=objList[[2]],subset=doublet=="Singlet"))
data3_cell<-WhichCells(subset(x=objList[[3]],subset=doublet=="Singlet"))

obj.list <- lapply(1:length(dirs),function(i){
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
  obj
})

obj.list[[1]]<-subset(x = obj.list[[1]], cells = data1_cell)
obj.list[[2]]<-subset(x = obj.list[[2]], cells = data2_cell)
obj.list[[3]]<-subset(x = obj.list[[3]], cells = data3_cell)

# callpeaks
DefaultAssay(obj.list[[1]])<-"ATAC"
DefaultAssay(obj.list[[2]])<-"ATAC"
DefaultAssay(obj.list[[3]])<-"ATAC"


combined.peaks <- UnifyPeaks(object.list = list(obj.list[[1]][["ATAC"]], obj.list[[2]][["ATAC"]],obj.list[[3]][["ATAC"]]), mode = "reduce")
data1.new_counts <- FeatureMatrix(
  fragments = Fragments(obj.list[[1]][["ATAC"]]),
  features = combined.peaks,
  sep = c(":", "-"),
  cells = colnames(obj.list[[1]][["ATAC"]])
)
data2.new_counts <- FeatureMatrix(
  fragments = Fragments(obj.list[[2]][["ATAC"]]),
  features = combined.peaks,
   sep = c(":", "-"),
  cells = colnames(obj.list[[2]][["ATAC"]])
)
data3.new_counts <- FeatureMatrix(
  fragments = Fragments(obj.list[[3]][["ATAC"]]),
  features = combined.peaks,
   sep = c(":", "-"),
  cells = colnames(obj.list[[3]][["ATAC"]])
)
obj.list[[1]][['ATAC']] <- CreateChromatinAssay(counts =data1.new_counts,sep = c(":", "-"),fragments = str_c(dirs[1],"atac_fragments.tsv.gz"),annotation = annotation)
obj.list[[2]][['ATAC']] <- CreateChromatinAssay(counts =data2.new_counts,sep = c(":", "-"),fragments = str_c(dirs[2],"atac_fragments.tsv.gz"),annotation = annotation)
obj.list[[3]][['ATAC']] <- CreateChromatinAssay(counts =data3.new_counts,sep = c(":", "-"),fragments = str_c(dirs[3],"atac_fragments.tsv.gz"),annotation = annotation)

##integrate RNA
for ( i in 1:length(obj.list)){
  DefaultAssay(obj.list[[i]]) <- "RNA"
  obj.list[[i]] <- NormalizeData(obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(obj.list[[i]], selection.method = "vst", nfeatures = 2000)
  DefaultAssay(obj.list[[i]]) <- "ATAC"
  obj.list[[i]] <- RunTFIDF(obj.list[[i]])
  obj.list[[i]] <- FindTopFeatures(obj.list[[i]], min.cutoff = "q0")
  obj.list[[i]] <- RunSVD(obj.list[[i]])
}
#find RNA integration features
DefaultAssay(obj.list[[1]]) <- "RNA"
DefaultAssay(obj.list[[2]]) <- "RNA"
DefaultAssay(obj.list[[3]]) <- "RNA"
rna.features <- SelectIntegrationFeatures(object.list = obj.list)
rna.anchors <- FindIntegrationAnchors(object.list = obj.list, anchor.features = rna.features)#CCA
#integrate
rna.combined <- IntegrateData(anchorset = rna.anchors,new.assay.name = "inteRNA")

##integrate ATAC
#构建merge object，获得未校正的lsi
obj.merge <- merge(
  x =obj.list[[1]],
  y =list(obj.list[[2]],obj.list[[3]]),
  merge.data = TRUE,
  project = "OSN"
)
DefaultAssay(obj.merge) <- "ATAC"
obj.merge <- FindTopFeatures(obj.merge, min.cutoff = 'q0')
obj.merge <- RunTFIDF(obj.merge)
obj.merge <- RunSVD(obj.merge)



#构建ATAC anchors
DefaultAssay(obj.list[[1]]) <- "ATAC"
DefaultAssay(obj.list[[2]]) <- "ATAC"
DefaultAssay(obj.list[[3]]) <- "ATAC"
atac.anchors <- FindIntegrationAnchors(
  object.list = obj.list,
  anchor.features = combined.peaks,
  reduction = "rlsi",
  dims = 2:30
)
#integrate
atac.combined <- IntegrateEmbeddings(
  anchorset = atac.anchors,
  reductions = obj.merge[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

#将rna和atac合并
combined<-rna.combined
combined[["integrated_lsi"]]<-atac.combined$integrated_lsi
#Gene expression data processing
DefaultAssay(combined) <- "inteRNA"
combined <- ScaleData(combined)
combined <- RunPCA(combined)
combined <- RunUMAP(combined,dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

#DNA accessibility data processing
DefaultAssay(combined) <- "ATAC"
combined <- RunUMAP(combined, reduction = "integrated_lsi", dims = 2:20, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

DefaultAssay(combined) <- "RNA"
#聚类和非线性降维
combined <- FindMultiModalNeighbors(combined, reduction.list = list("pca", "integrated_lsi"), dims.list = list(1:30, 2:20))
combined <- FindClusters(combined, graph.name = "wsnn", algorithm = 3, verbose = FALSE,resolution = 1)
combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
#注释前的cluster umap
rna_umap <- DimPlot(combined, reduction = "umap.rna",label = TRUE, label.size = 2, repel = TRUE, pt.size = 0.3) + ggtitle("RNA")
atac_umap <- DimPlot(combined, reduction = "umap.atac",  label = TRUE, label.size = 2, repel = TRUE, pt.size = 0.3) + ggtitle("ATAC")
wnn_umap <- DimPlot(combined, reduction = "wnn.umap",label = TRUE, label.size = 2, repel = TRUE, pt.size = 0.3) + ggtitle("WNN")
omp<-FeaturePlot(combined, reduction = "wnn.umap", feature='Omp')+ ggtitle("Omp")
Arc<-FeaturePlot(combined, reduction = "wnn.umap", feature='Arc')+ ggtitle("Arc")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/umap_after.pdf",rna_umap + atac_umap + wnn_umap & NoLegend() & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 30, height = 10, units = "cm", dpi = 300)
batch_umap <- DimPlot(combined, reduction = "wnn.umap",group.by = 'orig.ident',  pt.size = 0.3) + ggtitle("WNN")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/batch.pdf",batch_umap & theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 10, height = 10, units = "cm", dpi = 300)

#find markers
DefaultAssay(combined) <- "RNA"
combined.markers <- FindAllMarkers(combined,only.pos=TRUE)
filter.markers = combined.markers %>% select(gene, everything()) %>% subset(p_val<0.05)
filter.markers=filter.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC) -> top10

new.cluster.ids <- c("Mature OSNs", "Mature OSNs", "Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Immature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Mature OSNs","Macrophage","Neutrophils","Periglomerular cells","Sustentacular cells","INPs",
    "Gucy1b2+ cells", "Ensheathing glia","Mature OSNs")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
combined$celltype <- Idents(combined)
save(combined,file="~/ATAC-RNAseq/F1joint_mask/combined.RData")

# Figure 1C
# Umap after celltype annotation
celltype_color<-c("#96b6c7","#80a39f","#88898b","#b5aaab","#deb7b7","#bcb5d9","#8787b7","#886c93","#825176")
names(celltype_color)<-levels(combined)
p1 <- DimPlot(combined, reduction = "umap.rna",group.by="celltype",label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.3,cols=celltype_color) + ggtitle("RNA")
p2 <- DimPlot(combined, reduction = "umap.atac",group.by="celltype",  label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.3,cols=celltype_color) + ggtitle("ATAC")
p3 <- DimPlot(combined, reduction = "wnn.umap",group.by="celltype",  label = TRUE, label.size = 4, repel = TRUE, pt.size = 0.3,cols=celltype_color) + ggtitle("WNN")
p<-p1 + p2 + p3 & 
  theme_void()&
  labs(x="UMAP_1",y="UMAP_2") &
  tidydr::theme_dr(xlength = 0.2, 
           ylength = 0.2,
           arrow = arrow(length = unit(0.1, "inches"),type = "closed")) &
  theme(legend.position="right",plot.title = element_text(hjust = 0.5,size=10),text=element_text(size=8),axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.x = element_blank(),panel.grid = element_blank())
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/umap_after_celltype.pdf",p, device = "pdf", width =40, height = 10, units = "cm", dpi = 300)


# Figure S1C
cd_genes <- c("Omp","Gap43","Arhgap15","S100a9","Cntnap2","Cyp2g1","Sox11","Gucy1b2","Ptn")
dotplot<-DotPlot(object = combined, features = cd_genes,group.by="celltype")+scale_color_gradientn(colors =c("#FCFfF5","#D1Dbbd","#3E606F","#193441"))
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/dotplot.pdf",dotplot+theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8),axis.text.y=element_text(size=8,color = "black"), axis.text.x=element_text(size=8,  color = "black")), device = "pdf", width = 18, height = 8, units = "cm", dpi = 300)

# Figure 1D
# cell of each celltype
cellnumber<-table(Idents(combined))
cellnumber<-as.data.frame(cellnumber)
cellnumber$Group<-rep("total",9)
cell_number_p<-ggplot(cellnumber,aes(x=Var1,y=log10(Freq),fill=Var1))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = celltype_color) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Cellnumber")
ggsave(filename = "/data/R02/huangtt39/ATAC-RNAseq/F1joint_mask/celltype_pct.pdf",cell_number_p, device = "pdf", width =15, height = 10, units = "cm", dpi = 300)


#  Mature OSNs  recallpeak
mOSN_barcode<-WhichCells(combined, idents = "Mature OSNs")
data1_mOSN<-gsub('_1','',mOSN_barcode[grep('1$',mOSN_barcode)])
data2_mOSN<-gsub('_2','',mOSN_barcode[grep('2$',mOSN_barcode)])
data3_mOSN<-gsub('_3','',mOSN_barcode[grep('3$',mOSN_barcode)])

obj.list <- lapply(1:length(dirs),function(i){
  counts <- Read10X_h5(str_c(dirs[i],"filtered_feature_bc_matrix.h5")) 
  fragpath <- str_c(dirs[i],"atac_fragments.tsv.gz")
  metadata <- read.csv(str_c(dirs[i],"per_barcode_metrics.csv"),row.names=1) %>%dplyr::filter(is_cell==1) 
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
  DefaultAssay(obj)<-"ATAC"
  obj
})

obj.list[[1]]<-subset(x = obj.list[[1]], cells = data1_mOSN)
obj.list[[2]]<-subset(x = obj.list[[2]], cells = data2_mOSN)
obj.list[[3]]<-subset(x = obj.list[[3]], cells = data3_mOSN)

obj.list <- lapply(1:length(dirs),function(i){
  call_peak<-CallPeaks(obj.list[[i]])
  peak_counts <- FeatureMatrix(
    fragments = Fragments(obj.list[[i]][["ATAC"]]),
    features = call_peak,
    sep = c(":", "-"),
    cells = colnames(obj.list[[i]][["ATAC"]])
  )
  obj.list[[i]][['peak']] <- CreateChromatinAssay(counts =peak_counts,sep = c(":", "-"),fragments = str_c(dirs[i],"atac_fragments.tsv.gz"),annotation = annotation)
  obj.list[[i]]
})

# Callpeak merge
combined.peaks <- UnifyPeaks(object.list = list(obj.list[[1]][["peak"]], obj.list[[2]][["peak"]],obj.list[[3]][["peak"]]), mode = "reduce")
objList <- lapply(1:length(dirs),function(i){
  combined_peak_counts <- FeatureMatrix(
    fragments = Fragments(obj.list[[i]][["peak"]]),
    features = combined.peaks,
    sep = c(":", "-"),
    cells = colnames(obj.list[[i]][["peak"]])
  )
  obj.list[[i]][['peak']] <- CreateChromatinAssay(counts =combined_peak_counts,sep = c(":", "-"),fragments = str_c(dirs[i],"atac_fragments.tsv.gz"),annotation = annotation)
  obj.list[[i]]
})
Idents(combined)<-combined@meta.data$celltype
Combined_mOSN<-subset(combined,idents = "Mature OSNs")
merge_mOSN<-merge(x =obj.list[[1]],y =list(obj.list[[2]],obj.list[[3]]),merge.data = TRUE)
Combined_mOSN[["peak"]]<-merge_mOSN[["peak"]]
save(Combined_mOSN,file="~/ATAC-RNAseq/F1joint_mask/Combined_mOSN.RData")
