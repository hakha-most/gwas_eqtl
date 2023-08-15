### Script to extract GWAS hits from Neale lab data for categorical traits (i.e., not irnt traits)

iter=$1 #loop over trait index in a list of traits

codedir="path to the directory of scripts"
outdirs="path to the directory of output files" 
trait_list=$supp_data_dir/"ukb_neale_lab/cat.phenos.txt" # categorical traits info in Neale lab data (see Supp Data)
outdir=$outdirs/"cat_"$iter # trait index

mkdir -p $outdir
cd $outdir

line=`head -n $iter $trait_list | tail -n 1`
eval $line # download Neale lab GWAS sum stats

pheno=`echo $line | awk '{print $4}' | awk 'BEGIN{FS=".gwas"}{print $1}'`
outname=`echo $line | awk '{print $4}' | awk 'BEGIN{FS=".bgz"}{print $1}'`
mv $outdir/$outname".bgz" $outdir/$outname".gz" 

n_cols=`zcat $outdir/$outname".gz" | awk '{print NF}' | head -n 1`
outfile=$outdir/$pheno"_cat.hits.txt"

#extract hits and p-values
if [ "$n_cols" -eq 11 ]
then
zcat $outdir/$outname".gz" | awk -v pheno=$pheno 'BEGIN{FS="\t"}($11<0.00000005){print pheno,$1,$11}' > $outfile
fi

if [ "$n_cols" -eq 12 ]
then
zcat $outdir/$outname".gz" | awk -v pheno=$pheno 'BEGIN{FS="\t"}($12<0.00000005){print pheno,$1,$12}' > $outfile
fi

rm $outdir/$outname".gz"