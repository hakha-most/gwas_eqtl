args <- commandArgs(TRUE)
infile <- args[1]
outfile <- args[2]

df=read.csv(infile, sep='')
df=df[c("SNP","P")]

write.table(df,file=outfile,quote = FALSE,sep=" ",row.names = F, col.names = F)
