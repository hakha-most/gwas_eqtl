
infile="gencode.v39lift37.basic.annotation.gtf" #see URLs for download link

#############
#subset genes

awk 'BEGIN{FS="\t"}($3=="gene" && ($2=="HAVANA" || $2=="ENSEMBL"))' $infile > "gencode.genes"

cut -f1,4,5,7 "gencode.genes" > "temp1"
awk 'BEGIN{FS="gene_id"}{print $2}' "gencode.genes" | cut -d';' -f1 | cut -d'"' -f2 | cut -d"." -f1 > "temp2"
awk 'BEGIN{FS="gene_type"}{print $2}' "gencode.genes" | cut -d';' -f1 | cut -d'"' -f2 > "temp3"
awk 'BEGIN{FS="gene_name"}{print $2}' "gencode.genes" | cut -d';' -f1 | cut -d'"' -f2  > "temp4"
awk 'BEGIN{FS="hgnc_id"}{print $2}' "gencode.genes" | cut -d';' -f1 | cut -d'"' -f2  > "temp5"

paste "temp1" "temp2" "temp3" "temp4" "temp5" > "gencode.v39lift37.basic.annotation.processed.gtf"
rm "temp1" "temp2" "temp3" "temp4" "temp5" "gencode.genes"

#############
#subset protein-coding genes

ml R
script="filter_protein_coding_genes.R"
Rscript $script "gencode.v39lift37.basic.annotation.processed.gtf" "genes.protein_coding.v39.gtf" "genes.protein_coding.v39.1Mb_interval.bed"