args <- commandArgs(TRUE)
genefile <- args[1]
snpfile <- args[2]
eqtlfile <- args[3]
outfile <- args[4]

library(tidyverse)
library(data.table)

##########

d_pc_genes=fread(genefile)
d_pc_genes=d_pc_genes %>% select(GeneSymbol,gene,hgnc_id)
colnames(d_pc_genes)=c("GeneSymbol","GeneName","hgnc_id")

d_snp=fread(snpfile)
d_snp <- (d_snp) %>% separate(SNP,c("chr","pos","A1","A2"),sep=":",remove=F)

d_eqtl=fread(eqtlfile)
d_eqtl_pc=d_eqtl[d_eqtl$Gene %in% d_pc_genes$GeneSymbol,]
d_eqtl_pc_qc1=d_eqtl_pc[d_eqtl_pc$SNP %in% d_snp$rsid,]

d_eqtl_pc_qc1_snps=d_eqtl_pc_qc1[!duplicated(d_eqtl_pc_qc1$SNP),] %>% select(SNP,SNPChr,SNPPos,AssessedAllele,OtherAllele)
colnames(d_eqtl_pc_qc1_snps)=c("rsid","Chr","Pos","B1","B2")

dx=left_join(d_eqtl_pc_qc1_snps,d_snp,by="rsid")
dx$Chr=as.numeric(as.character(dx$Chr))
dx$chr=as.numeric(as.character(dx$chr))
dx$Pos=as.numeric(as.character(dx$Pos))
dx$pos=as.numeric(as.character(dx$pos))

dx1=dx[dx$Chr==dx$chr & dx$Pos==dx$pos & dx$A1==dx$B1 & dx$A2==dx$B2,]
dx2=dx[dx$Chr==dx$chr & dx$Pos==dx$pos & dx$A1==dx$B2 & dx$A2==dx$B1,]
dx_qc=rbind(dx1,dx2)
dx_qc=dx_qc %>% select(SNP,rsid)

d_eqtl_pc_qc2=d_eqtl_pc_qc1[d_eqtl_pc_qc1$SNP %in% dx_qc$rsid,]
d_eqtl_pc_qc2$rsid=d_eqtl_pc_qc2$SNP
d_eqtl_pc_qc2$GeneSymbol=d_eqtl_pc_qc2$Gene
d_eqtl_pc_qc2=d_eqtl_pc_qc2 %>% select(-SNP,-Gene)
d_eqtl_pc_qc2=left_join(d_eqtl_pc_qc2,dx_qc,by="rsid")
d_eqtl_pc_qc2=d_eqtl_pc_qc2 %>% select(-SNPChr,-SNPPos,-AssessedAllele,-OtherAllele)

write.table(d_eqtl_pc_qc2,file=outfile,quote = FALSE,sep="\t",row.names = F)
