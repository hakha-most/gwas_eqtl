#!/bin/sh

###########
## stratify annotated regions by the constraint metric of the linked gene 

codedir="06_get_roadmap_gene_enhancer_links_per_gene"
outdirs="Liu_et_al_roadmap_gene_enhancer_links/processed"

outdir=$outdirs"/by_gene_constraint"
mkdir $outdir

script=$codedir/"stratify_by_gene_annot.R"
Rscript $script 

###########
## merge regions by the constraint metric across all genes

for annot in {"LOEUF1","LOEUF2","LOEUF3","LOEUF4","LOEUF5"}
do

echo $annot

infile=$outdir/"all_regions."$annot"_genes.txt"
outfile=$outdir/"all_regions.merged_"$annot"_genes.bed"
$bedtools merge -i $infile > $outfile

done

