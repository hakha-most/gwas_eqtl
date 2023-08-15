#!/bin/sh

codedir="other_qtls/ieqtls"
outdirs="other_qtls/ieqtls"

tissue_list=$outdirs/"tissues_list.txt"

#===========concat all tissue/cell type pairs

outfile=$outdirs/"all.clumped.assoc"
rm $outfile

while IFS= read -r line 
do

params="`echo $line`"
tissue="`echo $params | cut -d ' ' -f 1`"
cell_type="`echo $params | cut -d ' ' -f 2`"

echo $tissue"."$cell_type

infile=$outdirs/$tissue"."$cell_type/$tissue"."$cell_type".clumped.snps"

cat $infile >> $outfile

done < $tissue_list


#======

infile=$outdirs/"all.clumped.assoc"
outfile=$outdirs/"ieqtls.clumped.txt"

ml R
script=$codedir/"filter_merge.R"
Rscript $script $infile $outfile
