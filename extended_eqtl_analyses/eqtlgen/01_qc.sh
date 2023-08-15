#!/bin/sh

codedir="extended_eqtl_analyses/eqtlgen"
outdir="extended_eqtl_analyses/eqtlgen"
eqtl_file="supp_note_data/eqtlgen/2019-12-11-cis-eQTLsFDR0.05-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz" #data on Zenodo

#filter SNPs and genes

genefile="gene_annotations/genes.protein_coding.v39.gtf"
snpfile="snp_annotations/maf01.kgp.snps" #data on Zenodo
outfile="supp_note_data/eqtlgen/eQTLs.pc_genes_ukb_snps.txt" #deposited on Zenodo

ml R
script=$codedir/"filter_snps.R"
Rscript $script $genefile $snpfile $eqtl_file $outfile

#Get eGenes=

infile=$outdir/"eQTLs.pc_genes_ukb_snps.txt"
outfile=$outdir/"eGenes.txt"

cat $infile | awk 'BEGIN{FS="\t"}{print $3}' | tail -n +2 | sort | uniq > $outfile

#split egenes into batches to perfrom LD clumping pre batch

save_dir=$outdir/"batches"
mkdir $save_dir

infile=$outdir/"eGenes.txt"

script=$codedir/"split_eGenes.R"
Rscript $script $save_dir $infile

