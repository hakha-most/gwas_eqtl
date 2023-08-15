#!/bin/sh


plink="plink" #link to the plink software

plinkdir="ld_dir" #link to LD ref files
outdir="colocalization"

snplist="colocalization/blood_trait.snps.txt"

for chr in {1..22}
do


ldfile=$plinkdir/"chr"$chr
outfile=$outdir/"chr"$chr

$plink --bfile $ldfile --r2 --ld-window 99999 --ld-window-kb 1000 --ld-window-r2 0.8 --ld-snp-list $snplist --out $outfile

done