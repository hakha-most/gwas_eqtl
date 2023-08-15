### Script to extract extract GWAS hits from Neale lab data for quantitative traits (rank-normalized, both sexes)

iter=$1 #loop over trait index in a list of traits

codedir="path to the directory of scripts"
outdirs="path to the directory of output files" 
trait_list=$supp_data_dir/"ukb_neale_lab/irnt.phenos.txt" # quantitative traits info in Neale lab data (see Supp Data)
outdir=$outdirs/"irnt_"$iter # trait index

mkdir -p $outdir
cd $outdir

line=`head -n $iter $trait_list | tail -n 1`
eval $line # download Neale lab GWAS sum stats

pheno=`echo $line | awk '{print $4}' | awk 'BEGIN{FS="_irnt"}{print $1}'`
outname=`echo $line | awk '{print $4}' | awk 'BEGIN{FS=".bgz"}{print $1}'`
mv $outdir/$outname".bgz" $outdir/$outname".gz" 

#extract hits and p-values
outfile=$outdir/$pheno"_irnt.hits.txt"
zcat $outdir/$outname".gz" | awk -v pheno=$pheno 'BEGIN{FS="\t"}($11<0.00000005){print pheno,$1,$11}' > $outfile

rm $outdir/$outname".gz"