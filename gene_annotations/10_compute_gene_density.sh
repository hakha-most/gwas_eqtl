
intervalfile="genes.protein_coding.v39.1Mb_interval.bed" #.bed file for +/- 1Mb of TSS of genes
annotfile="TSS.protein_coding.v39.bed" #.bed file for TSS of genes
outfile="genic_density_around_TSS.txt"

$bedtools intersect -a $intervalfile -b $annotfile -c | awk 'BEGIN{FS="\t"}{print $4,$5}' > $outfile