args <- commandArgs(TRUE)
infile <- args[1]
study <- args[2]
expr_type <- args[3]
outfile <- args[4]

library(tidyverse)
library(data.table)

df=fread(infile)

df=df %>% arrange(-pip)
df=df[!duplicated(df$cs_id),]
df=df %>% select(molecular_trait_id,SNP,pip,z)

df$study=study
df$type=expr_type

write.table(df,file=outfile,quote = FALSE,sep="\t",row.names = F)
