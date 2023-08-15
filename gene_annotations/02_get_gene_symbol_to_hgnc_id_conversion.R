library(data.table)
library(tidyverse)
library(stringr)

df=fread("Homo_sapiens.gene_info.gz")

df_filter=df %>% filter(str_detect(dbXrefs, 'HGNC')) #filter genes linked with HGNC ID
df_filter$hgnc_id=paste0("HGNC:",sub(".*HGNC:", "", df_filter$dbXrefs)) #extract HGNC ID

#redo extraction of HGNC ID for the subset of genes for which the previous step did not work due to different structure of the dbXrefs column
df_bad=df_filter %>% filter(str_detect(hgnc_id, '\\|')) 
df_bad$hgnc_id=sub("\\|.*", "", df_bad$hgnc_id)

df_good=df_filter[!(df_filter$hgnc_id %in% df_bad$hgnc_id),]
df_test=df_good %>% filter(str_detect(hgnc_id, '\\|'))

df_all=rbind(df_good %>% select(Symbol,Synonyms,hgnc_id) , df_bad %>% select(Symbol,Synonyms,hgnc_id))

df1=df_all[df_all$Synonyms=="-",] #genes with no Synonyms
df2=df_all[df_all$Synonyms!="-",] #genes with >1 Synonyms


###expand HGNC IDs to Synonyms
df2_rest=data.frame(Symbol=as.character(),hgnc_id=as.character())
for (i in 1:nrow(df2)){
    
  gene_rest=unlist(str_split(df2$Synonyms[i], pattern = "\\|"))
  
  d_out_temp=data.frame(Symbol=gene_rest,hgnc_id=df2$hgnc_id[i])
  df2_rest=rbind(df2_rest,d_out_temp)
  
}

# merge all
d_out_all=rbind((df1 %>% select(Symbol,hgnc_id)), (df2 %>% select(Symbol,hgnc_id)),df2_rest)

outfile="ncbi_genes.hgnc_id.txt"
write.table(d_out_all,file=outfile,quote=F,sep="\t",row.names=F)