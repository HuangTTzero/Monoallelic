#!/usr/bin/env Rscript
# packages ----------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))


Sys.time()
print("allele expression split")

# number crunching function -----------------------------------------------

load_reads <- function(filename){
  
  cols_to_read <- c(1,3,5,6)
  colname_vec <- c("readID","pos","cigar","seq")

  reads <- fread(file = filename,
                 sep = "\t",
                 header = F, fill = T,
                 select = cols_to_read, #read only necessary cls
                 col.names = colname_vec)
  return(reads)
}

variant_parsing <- function(reads, variant_positions){
  #parse all cigars to reference seq
  ops <- c("M", "=", "X")
  ranges_on_ref <- cigarRangesAlongReferenceSpace(reads$cigar, pos=reads$pos, ops=ops)
  ranges_on_query <- cigarRangesAlongQuerySpace(reads$cigar, ops=ops)
  gc(verbose = F)
  range_group <- togroup(PartitioningByWidth(ranges_on_ref))
  ranges_on_ref <- unlist(ranges_on_ref, use.names=FALSE)
  ranges_on_query <- unlist(ranges_on_query, use.names=FALSE)
  query2ref_shift <- start(ranges_on_ref) - start(ranges_on_query)

  var_pos <- variant_positions
  hits <- findOverlaps(var_pos, ranges_on_ref)
  hits_at_in_x <- var_pos[queryHits(hits)] - query2ref_shift[subjectHits(hits)]
  hits_group <- range_group[subjectHits(hits)]
  fetched_bases <- subseq(reads[hits_group,]$seq, start=hits_at_in_x, width=1L)

  #now add everything together in the output data.table
  out_vars <- data.table(
    obs_base = fetched_bases,
    pos = var_pos[queryHits(hits)]
  )
  out_vars[, c("readID") := reads[hits_group, c("readID"), with = F] ]
  out_vars <- out_vars[obs_base %in% c("A","C","G","T") ]
  setnames(out_vars,"pos","POS")

  return(out_vars)
}

calc_coverage_new <- function(vcf_chunk, out){
  chr <- unique(vcf_chunk$CHROM)
  print(paste("Starting to read data for chr ", chr))
  Sys.time()
  reads <- load_reads(filename = paste0(out,"var_overlap/",chr,".var_overlap.readsout"))

  print("Reading complete, processing reads & cigar values...")
  Sys.time()

  out_vars <- variant_parsing(reads, variant_positions = as.integer(vcf_chunk$POS))

  #crunch the numbers :-)
  out_vars <- merge(out_vars,vcf_chunk,by = "POS" )


  out_vars[          , basecall := "other"][
          obs_base == REF, basecall := "G1"][
          str_detect(string = c(ALT),pattern =obs_base), basecall := "G2"]

  out_reads <- out_vars[, .(readcall = read_decision(basecall)), by = c("readID")]
  rm(out_vars)
  print("Done!")
  return(out_reads)
}

read_decision <- function(basecalls){
    ux <- unique(basecalls)
    basecall_summary<-c()
    basecall_summary[c("G1","G2","other")]<-0
    basecall_summary[ux] <- tabulate(match(basecalls, ux))
    if(basecall_summary["G1"] > 0 & basecall_summary["G2"] > 0){
        return("CF")
      }
      else if(basecall_summary["G1"] > basecall_summary["G2"]){
        return("G1")
      }
      else if(basecall_summary["G1"] < basecall_summary["G2"]){
        return("G2")
      }
      else if(basecall_summary["G1"] ==0 & basecall_summary["G2"] == 0){
        return("UA")
      }
      else{
        print(paste0("The decision of read is error with G1:",basecall_summary["G1"],"G2:",basecall_summary["G2"]))
      }
}
# read_decision <- function(basecalls){
#   if(length(basecalls) == 1){
#     return(basecalls)
#   }else{
#     ux <- unique(basecalls)
#     basecall_summary <- tabulate(match(basecalls, ux))
#     names(basecall_summary) <- ux
#     majority_basecall <- ux[which.max(basecall_summary)]
#     if(basecall_summary[majority_basecall]/sum(basecall_summary) >= 0.66){
#       return(majority_basecall)
#     }else{
#       return("other")
#     }
#   }
# }


# startup variables -------------------------------------------------------
option_list <- list(
  make_option(c("-i", "--input"), type="character",
              help="Input sorted bam file. Mandatory"),
  make_option(c("-v", "--vcf"), type="character",
              help="SNP position list (VCF file) with variant annotation. Mandatory"),
  make_option(c("-o", "--out_dir"), type="character",
              help="Output directory. Mandatory"),
  make_option(c("-t", "--num_threads"), type="character",
              help="Number of threads to use. Mandatory")
)
opt <- parse_args(OptionParser(option_list=option_list))

if (any(is.null(opt$input),is.null(opt$out_dir),is.null(opt$num_threads),is.null(opt$vcf))) {
  stop("All mandatory parameters must be provided. See script usage (--help)")
}


path_snps <- opt$vcf
outpath <- opt$out_dir

if(!dir.exists(outpath)){
  try(system(paste("mkdir",outpath)))
}

ncores <- opt$num_threads

setwd(opt$out_dir)
setDTthreads(ncores)

# read stuff --------------------------------------------------------------
print("Reading Variants...")
if(grepl(path_snps, pattern = ".gz$")){
  vcf <- fread(cmd = paste("zcat",path_snps," | grep -v '^#'","| cut -f2,3,5,6"), col.names = c("CHROM","POS","REF","ALT"))
}else{
  vcf <- fread(cmd = paste("grep -v '^#'",path_snps,"|cut -f2,3,5,6"), col.names = c("CHROM","POS","REF","ALT"))
}

print("Done!")
Sys.time()


path_bam <- opt$input


# extract unique maps per chromosome  -------------------------------------
print("Extracting reads...")
if(!dir.exists(paste0(outpath,"var_overlap/"))){
  try(system(paste("mkdir",outpath,"var_overlap/")))
}
samtools_cmd1 <- paste0("samtools view -@ ",ncores)
samtools_cmd2 <- paste0(" | cut -f1,3,4,5,6,10 |awk '{if($4 == \"255\"){print > \"",outpath,"var_overlap/\"$2\".var_overlap.readsout\"}}'")
samtools_cmd <- paste(samtools_cmd1,path_bam,samtools_cmd2)
system(samtools_cmd)
print("Done")
Sys.time()


# crunch data ---------------------------------------------------------------

vcf_list <- split(vcf, by = "CHROM")
out_list <- lapply(vcf_list, function(x) calc_coverage_new(vcf_chunk = x, out = outpath))
out_reads <- rbindlist(out_list)


print("Finalizing converting & output ...")
Sys.time()
print(paste0("Total ",nrow(out_reads)," reads with strain snp"))
G1_reads <- out_reads%>%subset(readcall=="G1")%>%select(readID)
print(paste0(nrow(G1_reads)," reads is G1 specific"))
G2_reads <- out_reads%>%subset(readcall=="G2")%>%select(readID)
print(paste0(nrow(G2_reads)," reads is G2 specific"))
CF_reads <- out_reads%>%subset(readcall=="CF")%>%select(readID)
print(paste0(nrow(CF_reads)," reads is confliting"))
UA_reads <- out_reads%>%subset(readcall=="UA")%>%select(readID)
print(paste0(nrow(UA_reads)," reads is Unassignable"))
rm(out_reads)

print("Processing complete, writing output...")
Sys.time()
fwrite(G1_reads, file = paste0(outpath,"G1_reads.txt" ), sep= "\t", quote = F,col.names = FALSE)
fwrite(G2_reads, file = paste0(outpath,"G2_reads.txt" ), sep= "\t", quote = F,col.names = FALSE)
fwrite(CF_reads, file = paste0(outpath,"CF_reads.txt" ), sep= "\t",na = "NA", quote = F,col.names = FALSE)
fwrite(UA_reads, file = paste0(outpath,"UN_reads.txt" ), sep= "\t",na = "NA", quote = F,col.names = FALSE)
paste("DONE")
Sys.time()