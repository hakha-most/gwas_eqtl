###########
## compute enhancer features per gene 

codedir="06_get_roadmap_gene_enhancer_links_per_gene"
script=$codedir/"get_props.sh"

outdirs="Liu_et_al_roadmap_gene_enhancer_links/processed"
outdir=$outdirs"/all_props"

logdir=$outdir/"logs"
mkdir -p $logdir
cd $logdir

N=97 #97 jobs, 200 genes per job, for a total of 19240 genes 
job_name="get_props"
log_name=$job_name"_"%a".log"
sbatch -a 1-$N -o $log_name --cpus-per-task=1 --time=03:00:00 --mem-per-cpu=30gb --job-name=$job_name $script 


###########
## concatenate all 97 batches

outfile=$outdirs/"all_regions.props"
rm $outfile

for iter in {1..97}
do

infile=$outdir/"iter_"$iter".props"
cat $infile >> $outfile

done

rm -r $outdirs"/all_props"