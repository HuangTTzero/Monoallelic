import sys
from Bio import SeqIO
fasta=[]
with open(sys.argv[1], "r") as handle:
    i=0
    for seq_record in SeqIO.parse(handle, "fasta"):
        i+=1
        seq_record.id="chr"+i
        fasta.append(seq_record) 
#outfile='/data/R02/huangtt39/ATAC-RNAseq/snp/rf/mm10_modified.fa'
outfile='/data/R02/huangtt39/ATAC-RNAseq/snp/rf/chr_fa/mm10_chr.fa'
SeqIO.write(fasta, outfile, "fasta")
