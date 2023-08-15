#!/bin/sh

tissue=$1  	#performed for all tissue/cell type pairs
cell_type=$2

codedir="other_qtls/ieqtls"
outdir="other_qtls/ieqtls/"$tissue"."$cell_type
data_dir="GTEx/v8/ieQTL" #link to GWAS data, see URLs

#Get ieGenes==========

infile=$data_dir/$tissue"."$cell_type".ieQTL.eigenMT.annotated.txt.gz"
outfile=$outdir/$tissue"."$cell_type".v8.ieGenes.txt"

zcat $infile | cut -f2 | tail -n +2 | sort | uniq > $outfile

#lift hg38 stats file to hg19==========

infile=$data_dir/$tissue"."$cell_type".ieQTL.eigenMT.annotated.txt.gz"
outfile=$outdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg38.txt"
zcat $infile | cut -f1,2,16 | tail -n +2 | awk 'BEGIN{FS="\t"}{print $1,$2,$3,"snp"NR}' > $outfile

infile=$outdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg38.txt"
outfile=$outdir/$tissue"."$cell_type".pos_hg38.bed"
cat $infile | awk 'BEGIN{FS="_"}{print $1,$2,$5}' | awk '{print $1,$2-1,$2,$6}' > $outfile

liftover="liftOver" #link to software
changefile="hg38ToHg19.over.chain.gz" #link to chain file used in liftOver

infile=$outdir/$tissue"."$cell_type".pos_hg38.bed"
outfile1=$outdir/$tissue"."$cell_type".pos_hg19.bed"
outfile2=$outdir/$tissue"."$cell_type".pos_hg38.failed"
$liftover $infile $changefile $outfile1 $outfile2

maffile="snp_annotations/filter_snps.txt"
infile1=$outdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg38.txt"
infile2=$outdir/$tissue"."$cell_type".pos_hg38.bed"
infile3=$outdir/$tissue"."$cell_type".pos_hg19.bed"
outfile=$outdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg19.txt"

ml R
script=$codedir/"lift_snps.R"
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile

#clump==========

clumpdir=$outdir/"clump"
mkdir $clumpdir

cd $clumpdir
rm *

infile=$outdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg19.txt"
outfile=$clumpdir/$tissue"."$cell_type".v8.signif_variant_gene_pairs.hg19.txt"
cp $infile $outfile
awk -F' ' '{print > $2".txt"}' $outfile

rm $outfile
cd ../

plink="plink" #link to plink software
genelist=$outdir/$tissue"."$cell_type".v8.ieGenes.txt"
lddir="ldrefdir" #link to LD reference
script=$codedir/"get_clump_assoc.R"

concatfile=$outdir/$tissue"."$cell_type".clumped.snps"
rm $concatfile

#loop over genes
while IFS= read -r line 
do

gene="`echo $line`"

infile=$clumpdir/$gene".txt"
chr="`head -n 1 $infile | awk '{print $1}'`"
tempfile=$clumpdir/$gene".temp.txt"
echo "CHR GENE SNP P" > $tempfile
cat $infile >> $tempfile
mv $tempfile $infile

ldfile=$lddir/"chr"$chr
gwasfile=$clumpdir/$gene".txt"
outfile=$clumpdir/$gene
$plink --bfile $ldfile --clump $gwasfile --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out $outfile

infile=$clumpdir/$gene".clumped"
temp=$clumpdir/$gene".clumped.temp"
Rscript $script $infile $temp
awk -v tissue=$tissue -v cell_type=$cell_type -v gene=$gene '{print tissue,cell_type,gene,$0,NR}' $temp >> $concatfile
rm $temp

done < $genelist


