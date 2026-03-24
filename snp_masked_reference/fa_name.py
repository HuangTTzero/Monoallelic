import sys
from Bio import SeqIO
fasta=[]
with open(sys.argv[1], "r") as handle:
    for seq_record in SeqIO.parse(handle, "fasta"):
        description=seq_record.description
        chrom=description.split(' ')[1]
        chrom=chrom.split(':')[0]
        seq_record.id=chrom
        seq_record.description = seq_record.id
        fasta.append(seq_record)
outfile='/data/R02/huangtt39/F1_OSN/mapping/reference/rf/mm10-masked-modified.fa'
SeqIO.write(fasta, outfile, "fasta")
