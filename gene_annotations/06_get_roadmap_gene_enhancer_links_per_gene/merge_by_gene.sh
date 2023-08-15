#!/bin/sh

job_index=$SLURM_ARRAY_TASK_ID

outdirs="Liu_et_al_roadmap_gene_enhancer_links"
datadir=$outdirs/"raw_data"
outdir=$outdirs/"processed"

annot_file=$outdirs/"biosamples.txt"
annot="`tail -n+$job_index $annot_file | head -1`"

#temp dir per biosample to store intermediate files
annot_dir=$outdir/$annot
mkdir $annot_dir

cat $datadir/"links_"$annot"_6_2.5.txt" $datadir/"links_"$annot"_7_2.5.txt" $datadir/"links_"$annot"_12_2.5.txt" > $annot_dir/$annot".temp1"

#######
#filter protein-coding genes
pc_genes="genes.protein_coding.v39.genes" #ensembl IDs in "genes.protein_coding.v39.gtf"

cut -f4 $annot_dir/$annot".temp1" | sort | uniq > $annot_dir/$annot".genes"
comm -12 $annot_dir/$annot".genes" $pc_genes > $annot_dir/$annot".pc_genes"

#######
#merge annotated regions per gene
gene_list=$annot_dir/$annot".genes"

outfile=$outdir/$annot".txt"
rm $outfile

#loop over genes
while IFS= read -r line 
do

gene_temp="`echo $line`"

#get all annotated regions for the focal gene
grep -w $gene_temp $annot_dir/$annot".temp1" | sort -k1,1 -k2,2n > $annot_dir/$gene_temp".bed"

#merge regions
$bedtools merge -i $annot_dir/$gene_temp".bed" > $annot_dir/$gene_temp".bed_merge"

#write regions for all genes into one file
awk -v gene=$gene_temp -v annot=$annot 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3,gene,annot}' $annot_dir/$gene_temp".bed_merge" >> $outfile

#delete intermediate files
rm $annot_dir/$gene_temp".bed" $annot_dir/$gene_temp".bed_merge"

done < $gene_list

