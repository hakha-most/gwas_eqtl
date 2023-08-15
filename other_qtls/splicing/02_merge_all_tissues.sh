#!/bin/sh

outdirs="other_qtls/splicing"
codedir="other_qtls/splicing"

tissue_list=$outdirs/"tissues_list.txt"

###

outfile=$outdirs/"all_tissues.clumped.assoc"
rm $outfile

while IFS= read -r line 
do

tissue="`echo $line`"
echo $tissue

infile=$outdirs/$tissue/$tissue".clumped.snps"
cat $infile >> $outfile

done < $tissue_list


####

infile=$outdirs/"all_tissues.clumped.assoc"
outfile=$outdirs/"sqtls.clumped.txt"

ml R
script=$codedir/"filter_merge.R"
Rscript $script $infile $outfile
