### Script to LD clump eQTL hits for each eGene in each tissue

tissue=$1 #loop over tissues

supp_data_dir="path to the supp data" # see Supplementary Data
lddir="path to LD ref data" #UKB data used as LD ref panel
data_dir="path to the GTEx data"
codedir="path to the directory of scripts"
outdirs="path to the directory of output files"
outdir=$outdirs/$tissue
mkdir -p $outdir

#Get eGenes==========

infile=$data_dir/$tissue".v8.EUR.egenes.txt.gz"
outfile=$outdir/$tissue".v8.eGenes.txt"

zcat $infile | awk 'BEGIN{FS="\t"}($18<=0.05){print $1}' | sort | uniq > $outfile

#lift hg38 stats file to hg19==========

infile=$data_dir/$tissue".v8.EUR.signif_pairs.txt.gz"
outfile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
zcat $infile | awk 'BEGIN{FS="\t"}{print $2,$1,$7}' | tail -n +2 | awk '{print $1,$2,$3,"snp"NR}' > $outfile

infile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
outfile=$outdir/$tissue".pos_hg38.bed"
cat $infile | awk 'BEGIN{FS="_"}{print $1,$2,$5}' | awk '{print $1,$2-1,$2,$6}' > $outfile

changefile=$supp_data_dir/"auxiliary_files/hg38ToHg19.over.chain.gz" 
maffile=$supp_data_dir/"snp_annots/maf01.kgp.snps"

infile=$outdir/$tissue".pos_hg38.bed"
outfile1=$outdir/$tissue".pos_hg19.bed"
outfile2=$outdir/$tissue".pos_hg38.failed"
$liftover $infile $changefile $outfile1 $outfile2 #using LiftOver software

infile1=$outdir/$tissue".v8.signif_variant_gene_pairs.hg38.txt"
infile2=$outdir/$tissue".pos_hg38.bed"
infile3=$outdir/$tissue".pos_hg19.bed"
outfile=$outdir/$tissue".v8.signif_variant_gene_pairs.hg19.txt"

ml R
script=$codedir/"lift_snps.R" #restrict to eQTLs that match KGP SNPs with MAF>0.01
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile

#LD clumping==========

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

genelist=$outdir/$tissue".v8.eGenes.txt"
script=$codedir/"get_clump_assoc.R"

concatfile=$outdir/$tissue".clumped.snps"
rm $concatfile

#loop over eGenes
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
$plink --bfile $ldfile --clump $gwasfile --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out $outfile #using plink software

infile=$clumpdir/$gene".clumped"
temp=$clumpdir/$gene".clumped.temp"
Rscript $script $infile $temp
awk -v tissue=$tissue -v gene=$gene '{print tissue,gene,$0,NR}' $temp >> $concatfile
rm $temp

done < $genelist


