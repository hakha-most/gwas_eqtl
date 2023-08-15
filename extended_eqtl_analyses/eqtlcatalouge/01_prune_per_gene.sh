#!/bin/sh

iter=$SLURM_ARRAY_TASK_ID #performed for each study

codedir="extended_eqtl_analyses/eqtlcatalouge"
outdirs="extended_eqtl_analyses/eqtlcatalouge"
datadir="supp_note_data/eqtl_catalogue/data_download" #data on Zenodo, downloaded from https://www.ebi.ac.uk/eqtl/

liftover="liftOver" #link to software
changefile="hg38ToHg19.over.chain.gz" #link to chain file used in liftOver
ml R

datafile=$outdirs/"studies.txt"
line=`head -n $iter $datafile | tail -n 1`
study=`echo $line | awk '{print $1}'`
type=`echo $line | awk '{print $2}'`
study_type=$study"_"$type

outdir=$outdirs/$study_type
mkdir $outdir
cd $outdir

#loop over eQTL types

for expr_type in {"ge","exon"}
do

echo $expr_type

#lift hg38 stats file to hg19==========

infile=$datadir/$study"_"$expr_type"_"$type".purity_filtered.txt.gz"
outfile=$outdir/$expr_type".variant_gene_pairs.hg38.txt"
zcat $infile | awk 'BEGIN{FS="\t"}{print $1,$3,$4}' | tail -n +2 | awk '{print $1,$2,$3,"snp"NR}' > $outfile

infile=$outdir/$expr_type".variant_gene_pairs.hg38.txt"
outfile=$outdir/$expr_type".pos_hg38.bed"
cat $infile | awk '{print $2,$3-1,$3,$4}' > $outfile

temp=$outdir/$expr_type"_tmp"
awk 'BEGIN{OFS=""}{print "chr",$0}' $outfile > $temp
mv $temp $outfile

infile=$outdir/$expr_type".pos_hg38.bed"
outfile1=$outdir/$expr_type".pos_hg19.bed"
outfile2=$outdir/$expr_type".pos_hg38.failed"
$liftover $infile $changefile $outfile1 $outfile2

maffile="snp_annotations/maf01.kgp.snps" #data on Zenodo
infile1=$datadir/$study"_"$expr_type"_"$type".purity_filtered.txt.gz"
infile2=$outdir/$expr_type".pos_hg38.bed"
infile3=$outdir/$expr_type".pos_hg19.bed"
outfile=$outdir/$expr_type".variant_gene_pairs.hg19.txt"

script=$codedir/"lift_snps.R"
Rscript $script $maffile $infile1 $infile2 $infile3 $outfile

#get top variant for each credible set==========

infile=$outdir/$expr_type".variant_gene_pairs.hg19.txt"
outfile=$outdir/$expr_type".clumped.assoc.txt"

script=$codedir/"get_clump_assoc.R"
Rscript $script $infile $study_type $expr_type $outfile


done

