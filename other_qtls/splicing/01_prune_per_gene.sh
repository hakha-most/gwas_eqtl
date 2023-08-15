#!/bin/sh

tissue=$1 #loop over GTEx tissues

codedir="other_qtls/splicing"
outdir="other_qtls/splicing/"$tissue
data_dir="GTEx/v8/EUR/sqtls" #link to GWAS data, see URLs

########

#lift hg38 stats file to hg19==========

infile=$data_dir/$tissue".v8.EUR.signif_pairs.txt.gz"
outfile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
zcat $infile | awk 'BEGIN{FS="\t"}{print $2,$1,$7}' | tail -n +2 | awk '{print $1,$2,$3,"snp"NR}' > $outfile

infile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
outfile=$outdir/$tissue".pos_hg38.bed"
cat $infile | awk 'BEGIN{FS="_"}{print $1,$2,$6}' | awk '{print $1,$2-1,$2,$5}' > $outfile

liftover="liftOver" #link to software
changefile="hg38ToHg19.over.chain.gz" #link to chain file used in liftOver

infile=$outdir/$tissue".pos_hg38.bed"
outfile1=$outdir/$tissue".pos_hg19.bed"
outfile2=$outdir/$tissue".pos_hg38.failed"
$liftover $infile $changefile $outfile1 $outfile2

maffile="snp_annotations/maf01.kgp.snps"
infile1=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
infile2=$outdir/$tissue".pos_hg38.bed"
infile3=$outdir/$tissue".pos_hg19.bed"
outfile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg19.txt"

ml R
script=$codedir/"lift_snps.R"
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile

############

script=$codedir/"filter_sGenes.R"

snpfile="snp_annotations/filter_snps.txt"
pcfile="gene_annotations/genes.protein_coding.v39.gtf"
infile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg19.txt"
outfile1=$outdir/$tissue".v8.sqtl_signifpairs.filter.txt"
outfile2=$outdir/$tissue".v8.sGenes.txt"

ml R
Rscript $script $snpfile $pcfile $infile $outfile1 $outfile2

#clump==========

clumpdir=$outdir/"clump"
mkdir $clumpdir

cd $clumpdir
rm *

infile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg19.txt"
outfile=$clumpdir/$tissue".v8.signif_variant_gene_pairs.hg19.txt"
cp $infile $outfile
awk -F' ' '{print > $2".txt"}' $outfile

rm $outfile
cd ../

plink="plink" #link to plink software
genelist=$outdir/$tissue".v8.sGenes.txt"
lddir="ldrefdir" #link to LD reference
script=$codedir/"get_clump_assoc.R"

concatfile=$outdir/$tissue".clumped.snps"
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
awk -v tissue=$tissue -v gene=$gene '{print tissue,gene,$0,NR}' $temp >> $concatfile
rm $temp

done < $genelist




