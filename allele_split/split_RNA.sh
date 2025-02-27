# For RNA-seq, the bam file from cellranger output can be used directly, and then allele-specific assignment should be perform
SNP="snp/align_snp_RNA.txt"
file='mapping_mask/dateset_mask/outs/gex_possorted_bam.bam'
outdir='mapping/jointF1OSN/split/RNA/'
Rscript software/SNPsplit_RNA/SNPsplit_RNA.R -i $file -v $SNP -t 24 -o $outdir

