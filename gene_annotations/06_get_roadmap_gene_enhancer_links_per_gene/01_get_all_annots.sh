###########
## for each gene merge enhancer-related annotated regions per biosamples 
## all input and ouput data files are on Zenodo in dir "gene_annots"

codedir="06_get_roadmap_gene_enhancer_links_per_gene"
outdirs="Liu_et_al_roadmap_gene_enhancer_links"
outdir=$outdirs/"processed"

logdir=$outdir/"logs"
mkdir -p $logdir
cd $logdir

script=$codedir/"merge_by_gene.sh"
job_name="merge_by_gene"
log_name=$job_name"_"%a".log"
sbatch -a 1-111 -o $log_name --cpus-per-task=1 --time=12:00:00 --mem-per-cpu=20gb --job-name=$job_name $script


###########
## concatenate all regions across genes and biosamples into one file 

outfile=$outdir/"all_regions.txt"
rm $outfile

while IFS= read -r line 
do

annot="`echo $line`"
echo $annot

cat $outdir/$annot".txt" >> $outfile

rm $outdir/$annot".txt"
rm -r $outdir/$annot

done < $outdirs/"biosamples.txt"


###########
## for each gene merge regions across biosample 

bed_file=$outdir/"all_regions.txt"

annot_dir=$outdir/"all_biosamples"
mkdir $annot_dir

#get list of all genes
gene_list=$annot_dir/"all.genes"
cut -f4 $bed_file | sort | uniq > $gene_list

#loop over genes

outfile=$outdir/"all_regions.merged_biosamples.txt"
rm $outfile

while IFS= read -r line 
do

gene_temp="`echo $line`"

grep -w $gene_temp $bed_file | sort -k1,1 -k2,2n > $annot_dir/$gene_temp".bed"
$bedtools merge -i $annot_dir/$gene_temp".bed" > $annot_dir/$gene_temp".bed_merge"
awk -v gene=$gene_temp 'BEGIN{FS="\t"; OFS="\t"}{print $1,$2,$3,gene}' $annot_dir/$gene_temp".bed_merge" >> $outfile

rm $annot_dir/$gene_temp".bed" $annot_dir/$gene_temp".bed_merge"

done < $gene_list

sort -k1,1 -k2,2n $outdir/"all_regions.merged_biosamples.txt" > $outdir/"temp"
mv $outdir/"temp" $outdir/"all_regions.merged_biosamples.txt" 

rm -r $annot_dir




