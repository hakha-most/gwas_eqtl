### Script to LD clump GWAS hits for categorical traits (i.e., not irnt traits)

iter=$1 #loop over trait index in a list of traits

supp_data_dir="path to the supp data" # see Supplementary Data
lddir="path to LD ref data" #UKB data used as LD ref panel
codedir="path to the directory of scripts"
outdirs="path to the directory of output files" 
trait_list=$supp_data_dir/"ukb_neale_lab/cat.phenos.txt" # categorical traits info in Neale lab data
outdir=$outdirs/"cat_"$iter

mkdir -p $outdir
cd $outdir

line=`head -n $iter $trait_list | tail -n 1`
pheno=`echo $line | awk '{print $4}' | awk 'BEGIN{FS=".gwas"}{print $1}'`

#restrict to eQTLs that match KGP SNPs with MAF>0.01
maffile=$supp_data_dir/"snp_annots/maf01.kgp.snps"
infile=$outdir/$pheno"_cat.hits.txt"
outfile=$outdir/$pheno"_cat.hits.txt"
script=$codedir/"prep_for_clump.R"
Rscript $script $infile $maffile $outfile


#LD clumping
gwasfile=$outdir/$pheno"_cat.hits.txt"
for chr in {1..22}
do

ldfile=$lddir/"chr"$chr
outfile=$outdir/"chr"$chr
$plink --bfile $ldfile --clump $gwasfile --clump-p1 0.00000005 --clump-r2 0.1 --clump-kb 1000 --out $outfile

done

#concatenate chrs
clumpdir=$outdirs/"clumped_hits"
outfile=$clumpdir/$pheno"_cat.hits.txt"
rm $outfile

script=$codedir/"get_clump_assoc.R"

for chr in {1..22}
do

infile=$outdir/"chr"$chr".clumped"
temp=$outdir/"chr"$chr".snps.temp"

Rscript $script $infile $temp

cat $temp >> $outfile
rm $temp

done #chr

