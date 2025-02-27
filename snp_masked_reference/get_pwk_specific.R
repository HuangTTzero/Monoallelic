pwk_snp = read.table("/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/all_PWK_SNPs.txt",sep="\t")#读取文本数据
names(pwk_snp)=c("ID","chr","pos","strand","Allele")

p2_snp = read.table("/data/R02/huangtt39/ATAC-RNAseq/snp/129/all_129_SNPs.txt",sep="\t")#读取文本数据
names(p2_snp)=c("ID","chr","pos","strand","Allele")

pwk_vcf = read.table("/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/origin_PWK_SNPs.txt",sep="\t")#读取文本数据
names(pwk_vcf)=c("chr","pos","id","ref","alt","qual","filter","info","format","pwk_phj")

library(sqldf)
pwk_specific=sqldf("select pwk_snp.* from pwk_snp left join p2_snp on pwk_snp.chr=p2_snp.chr and pwk_snp.pos=p2_snp.pos and pwk_snp.Allele=p2_snp.Allele where p2_snp.ID is NULL")
pwk_specific_vcf=sqldf("select pwk_vcf.* from pwk_vcf join pwk_specific on pwk_vcf.chr=pwk_specific.chr and pwk_vcf.pos=pwk_specific.pos")
pwk_specific_p2=sqldf("select p2_snp.* from pwk_specific join p2_snp on pwk_specific.chr=p2_snp.chr and pwk_specific.pos=p2_snp.pos")
write.table(x=pwk_specific, file = "/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/pwk_specific_snp.txt", quote = FALSE, sep = "\t")
write.table(x=pwk_specific_vcf, file = "/data/R02/huangtt39/ATAC-RNAseq/snp/pwk/pwk_specific.vcf", quote = FALSE, sep = "\t")
write.table(x=pwk_specific_p2, file = "/data/R02/huangtt39/ATAC-RNAseq/snp/129/pwk_specific_129_snp.txt", quote = FALSE, sep = "\t")