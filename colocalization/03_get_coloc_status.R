set.seed(0)

library(data.table)
library(tidyverse)
library(dplyr)

indir="colocalization"
snpfile="snp_annotations/maf01.kgp.snps"
pcfile="gene_annotations/genes.protein_coding.v39.gtf"
eqtlgenfile="extended_eqtl_analyses/eqtlgen/eQTLs.pc_genes_ukb_snps.txt"
gtexfile="colocalization/Whole_Blood.v8.signif_variant_gene_pairs.hg19.txt" #GTEx whole blood eQTLs, no clumping
outfile=paste0(indir,"/coloc_status.txt")

d_qc_snps=fread(snpfile)
d_pc_genes=fread(pcfile)
d_gtex=fread(gtexfile,header=F)
d_gen=fread(eqtlgenfile)

d_snp=fread(paste0(indir,"/blood_trait.snps.txt")) #blood-trait GWAS SNPs
d_snp <- (d_snp) %>% separate(SNP,c("chr","pos","A1","A2"),sep=":",remove=F) 

#qc GTEx eQTLs
colnames(d_gtex)=c("chr","Gene","SNP","Pval")
d_gtex <- (d_gtex) %>% separate(Gene,c("GeneSymbol","nn"),sep="\\.",remove=F) %>% select(-nn)
d_gtex <- d_gtex %>% select(GeneSymbol,SNP,Pval)
d_gtex=d_gtex[d_gtex$GeneSymbol %in% d_pc_genes$GeneSymbol,]
d_gtex=d_gtex[d_gtex$SNP %in% d_qc_snps$SNP,]
d_gtex=d_gtex %>% arrange(Pval)

#qc eQTLGen eQTLs
d_gen <- d_gen %>% select(GeneSymbol,SNP,Pvalue)
colnames(d_gen)=c("GeneSymbol","SNP","Pval")
d_gen=d_gen[d_gen$GeneSymbol %in% d_pc_genes$GeneSymbol,]
d_gen=d_gen[d_gen$SNP %in% d_qc_snps$SNP,]
d_gen=d_gen %>% arrange(Pval)

pval_qs=quantile(d_gen$Pval,seq(0,1,0.1))[2:11] #p-val deciles

###

outfile_temp=paste0(indir,"/gtex_eqtls.txt")
write.table(d_gtex,file=outfile_temp,quote=F,sep="\t",row.names=F)

outfile_temp=paste0(indir,"/eqtlgen_eqtls.txt")
write.table(d_gen,file=outfile_temp,quote=F,sep="\t",row.names=F)

###
#function to evalue coloc status with GWAS SNPs
get_coloc=function(d_ld,d_eqtl){

	d_temp=d_eqtl[d_eqtl$SNP %in% d_ld$SNP_B,] %>% select(SNP,Pval)
	colnames(d_temp)=c("SNP_B","Pval")
	dx=left_join(d_ld,d_temp,by="SNP_B")

	d1=dx[is.na(dx$Pval),]
	d2=dx[!(is.na(dx$Pval)),] %>% arrange(Pval)
	d2=d2[!duplicated(d2$SNPs),]

	d1=d1[!duplicated(d1$SNP_A),]
	d2=d2[!duplicated(d2$SNP_A),]
	d1=d1[!(d1$SNP_A %in% d2$SNP_A),]
	 
	d_out=rbind(d1,d2) %>% select(SNP_A,Pval)
	colnames(d_out)=c("GWAS_SNP","eQTL_P")

	return(d_out)
}


#function to evalue coloc status with GWAS SNPs by decile of p-values in eQTLGen
get_coloc_by_pcat=function(d_ld,d_eqtl,pval_qs){

	p_temp=1
	d_temp=d_eqtl[d_eqtl$Pval<=pval_qs[p_temp],]
	d_out=get_coloc(d_ld,d_temp)
	colnames(d_out)=c("GWAS_SNP",paste0("eQTL_P",as.character(p_temp)))

	for (p_temp in 2:10){

		d_temp=d_eqtl[d_eqtl$Pval<=pval_qs[p_temp],]
		d_out_temp=get_coloc(d_ld,d_temp)
		colnames(d_out_temp)=c("GWAS_SNP",paste0("eQTL_P",as.character(p_temp)))
		d_out=left_join(d_out,d_out_temp,by="GWAS_SNP")

	}

	return(d_out)
}

###
#loop over chromosomes

chr=1
d_ld=fread(paste0(indir,"/chr",as.character(chr),".ld"))
d_ld$SNPs=paste0(d_ld$SNP_A,"_",d_ld$SNP_B)

d_out1=get_coloc(d_ld,d_gtex)
colnames(d_out1)=c("GWAS_SNP","GTEx_P")
d_out2=get_coloc_by_pcat(d_ld,d_gen,pval_qs)
d_out=left_join(d_out1,d_out2,by="GWAS_SNP")

#
for (chr in 2:22) {

	print(as.character(chr))

	d_ld=fread(paste0(indir,"/chr",as.character(chr),".ld"))
	d_ld$SNPs=paste0(d_ld$SNP_A,"_",d_ld$SNP_B)

	d_out1=get_coloc(d_ld,d_gtex)
	colnames(d_out1)=c("GWAS_SNP","GTEx_P")
	d_out2=get_coloc_by_pcat(d_ld,d_gen,pval_qs)
	d_out_temp=left_join(d_out1,d_out2,by="GWAS_SNP")

	#####
	d_out=rbind(d_out,d_out_temp)

}

write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)


