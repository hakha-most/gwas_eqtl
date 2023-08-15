library(data.table)
library(tidyverse)
library(dplyr)
library(topGO)
library(biomaRt)


#####from GENCODE

#set of ~18K protein-coding genes to be used as background
d_pc_genes=fread("genes.protein_coding.v39.gtf")
bg_genes=unique(as.character(d_pc_genes$hgnc_id)) #using HGNC IDs as gene identifier


#####use biomaRt

#get most sepecific GO annotations per gene
db=useMart('ENSEMBL_MART_ENSEMBL',dataset='hsapiens_gene_ensembl', host="www.ensembl.org")
go_ids=getBM(attributes=c('go_id', 'hgnc_id','external_gene_name', 'namespace_1003'), filters='hgnc_id', values=bg_genes, mart=db)

write.table(go_ids,file="12_GO_annotations/all_GO_annots.txt",quote=F,sep="\t",row.names=F)

#####run topGO

gene_2_GO=unstack(go_ids[,c(1,2)])
candidate_list=sample(bg_genes,2000) #random set of genes, selected to build GO DAGs
geneList=factor(as.integer(bg_genes %in% candidate_list))
names(geneList)= bg_genes

#build DAG, for "Biological Processes" (ontology) conditioning on GO categories with >400 genes (nodeSize)
GOdata=new('topGOdata', ontology="BP", allGenes = geneList, nodeSize=400, annot = annFUN.gene2GO, gene2GO = gene_2_GO)

#test enrichment run
classic_fisher_result=runTest(GOdata, algorithm='classic', statistic='fisher') 
allGO=usedGO(GOdata)
d_BP=GenTable(GOdata, Fisher=classic_fisher_result, orderBy='Fisher', topNodes=length(allGO))
d_BP=d_BP[c("GO.ID","Term","Annotated")]

write.table(d_BP,file="12_GO_annotations/GO_names_BP_min400.txt",quote=F,sep="\t",row.names=F)

#####get genes per GO category
allGO_genes = genesInTerm(GOdata)

d_out=data.frame(gene=bg_genes)
GO_list=names(allGO_genes)

i=0
for (GO_temp in GO_list){
  i=i+1
  print(i)
  
  i=which(GO_list==GO_temp)
  d_list_temp=unname(allGO_genes[i])[[1]]
  d_temp=data.frame(gene=bg_genes,go=0)
  d_temp[d_temp$gene %in% d_list_temp,]$go=1
  colnames(d_temp)=c("gene",GO_temp)
  d_out=left_join(d_out,d_temp,by="gene")
}

outfile="12_GO_annotations/genes_BP_min400.txt"
write.table(d_out,file=outfile,quote=F,sep="\t",row.names=F)
