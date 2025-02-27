~/software/gatk-4.3.0.0/gatk IndexFeatureFile -I /data/R02/huangtt39/ATAC-RNAseq/snp/merge_snp/merge_snp2mask.vcf
samtools faidx /data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/Mus_musculus.GRCm38.dna.primary_assembly.modified.fa
picard CreateSequenceDictionary R=/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/Mus_musculus.GRCm38.dna.primary_assembly.modified.fa O=/data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/Mus_musculus.GRCm38.dna.primary_assembly.modified.dict
~/software/gatk-4.3.0.0/gatk FastaAlternateReferenceMaker -R /data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10-2020-A-build/Mus_musculus.GRCm38.dna.primary_assembly.modified.fa -O /data/R02/huangtt39/ATAC-RNAseq/mapping/reference/rf/mm10_masked.fa -V /data/R02/huangtt39/ATAC-RNAseq/snp/merge_snp/merge_snp2mask.vcf
# FastaAlternateReferenceMaker will change the name of chromsome,the fa_name.py is to change the name of each chromsome from 1 to chr1
python /data/R02/huangtt39/ATAC-RNAseq/snp/myscript/rf/fa_name.py /data/R02/huangtt39/ATAC-RNAseq/mapping/reference/mm10_masked.fa
