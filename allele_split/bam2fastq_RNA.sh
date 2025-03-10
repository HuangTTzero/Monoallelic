outdir='/data/R02/huangtt39/ATAC-RNAseq/mapping/dataset/split/RNA/'
ori_head='/data/R02/huangtt39/ATAC-RNAseq/mapping/dataset/fastq/RNA/dataset_RNA_total'
out_head1='dataset-scRNA-G1_S1'
out_head2='dataset-scRNA-G2_S1'
seq_outdir='/data/R02/huangtt39/ATAC-RNAseq/mapping/dataset/split/RNA/'
mkdir $seq_outdir'g1_fq/'
mkdir $seq_outdir'g2_fq/'
#根据ID 提取序列
seqkit grep --pattern-file $outdir'G1_reads.txt' ${ori_head}_I1.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g1_fq/'$out_head1'_L002_I1_001.fastq.gz'
seqkit grep --pattern-file $outdir'G1_reads.txt' ${ori_head}_I2.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g1_fq/'$out_head1'_L002_I2_001.fastq.gz'
seqkit grep --pattern-file $outdir'G1_reads.txt' ${ori_head}_R1.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g1_fq/'$out_head1'_L002_R1_001.fastq.gz'
seqkit grep --pattern-file $outdir'G1_reads.txt' ${ori_head}_R2.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g1_fq/'$out_head1'_L002_R2_001.fastq.gz'

seqkit grep --pattern-file $outdir'G2_reads.txt' ${ori_head}_I1.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g2_fq/'$out_head2'_L002_I1_001.fastq.gz'
seqkit grep --pattern-file $outdir'G2_reads.txt' ${ori_head}_I2.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g2_fq/'$out_head2'_L002_I2_001.fastq.gz'
seqkit grep --pattern-file $outdir'G2_reads.txt' ${ori_head}_R1.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g2_fq/'$out_head2'_L002_R1_001.fastq.gz'
seqkit grep --pattern-file $outdir'G2_reads.txt' ${ori_head}_R2.fastq.gz|seqkit sort -n|gzip >$seq_outdir'g2_fq/'$out_head2'_L002_R2_001.fastq.gz'
