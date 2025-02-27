# For ATAC-seq, the sequences need to be mapped to the SNP-masked reference using Bowtie2, and then allele-specific assignment should be perform
SNP="snp/align_snp.txt"
file='mapping/dataset/mapping_ATAC/ATAC_total.bam'
outfile='mapping/dataset/split/ATAC/SNPsplit/ATAC_SNPsplit'
outdir='mapping/dataset/split/ATAC/SNPsplit'
mkdir $outdir
software/SNPsplit/SNPsplit_ATAC --no_sort --snp_file $SNP --paired  --conflicting $file  -o $outdir 1>$outfile.log 2>$outfile.error

