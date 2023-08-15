#!/bin/sh

iter=$SLURM_ARRAY_TASK_ID #performed for all batches of eGenes

codedir="extended_eqtl_analyses/eqtlgen"
outdirs="extended_eqtl_analyses/eqtlgen"

eqtlfile="supp_note_data/eqtlgen/eQTLs.pc_genes_ukb_snps.txt" #data on Zenodo

genelist=$outdirs/"batches/batch_"$iter".txt"
clumpdir=$outdirs/"clump/batch_"$iter
mkdir $clumpdir

plink="plink" #link to plink software
lddir="ldrefdir" #link to LD reference

ml R

#####

cd $clumpdir
rm *

script=$codedir/"get_clump_assoc.R"

concatfile=$outdirs/"clump/batch_"$iter".clumped.snps"
rm $concatfile

#loop over genes
while IFS= read -r line 
do

gene="`echo $line`"

infile=$clumpdir/$gene".txt"
grep $gene $eqtlfile > $infile

chr="`head -n 1 $infile | cut -f4`"
tempfile=$clumpdir/$gene".temp.txt"
echo "CHR GENE SNP P" > $tempfile
awk 'BEGIN{FS="\t"}{print $4,$3,$11,$1}' $infile >> $tempfile
mv $tempfile $infile

ldfile=$lddir/"chr"$chr
gwasfile=$clumpdir/$gene".txt"
outfile=$clumpdir/$gene
$plink --bfile $ldfile --clump $gwasfile --clump-p1 1 --clump-p2 1 --clump-r2 0.1 --clump-kb 1000 --out $outfile

infile=$clumpdir/$gene".clumped"
temp=$clumpdir/$gene".clumped.temp"
Rscript $script $infile $temp
awk -v gene=$gene '{print gene,$0,NR}' $temp >> $concatfile
rm $temp

done < $genelist


