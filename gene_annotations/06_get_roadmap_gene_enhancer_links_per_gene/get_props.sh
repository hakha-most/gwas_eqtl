#!/bin/sh

iter=$SLURM_ARRAY_TASK_ID

indir="Liu_et_al_roadmap_gene_enhancer_links/processed"
outdir=$indir"/all_props"

enh_file=$indir/"all_regions.txt"
outfile=$outdir/"iter_"$iter".props"

ml R
codedir="/oak/stanford/groups/pritch/users/hakha/data_sets/roadmap_gene_enhancer_links/codes"
script=$codedir/"get_props.R"
Rscript $script $enh_file $iter $outfile