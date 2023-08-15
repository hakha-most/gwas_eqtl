
outdirs="extended_eqtl_analyses/eqtlgen"
indir=$outdirs/"clump"

#

outfile="supp_note_data/eqtlgen/eqtlgen_eqtls.clumped.txt" #data on Zenodo
rm $outfile

for iter in {1..N_batch} #N_batch=number of batches
do

infile=$indir/"batch_"$iter".clumped.snps"
N=`wc -l $infile`

cat $infile >> $outfile
echo $iter $N

done


#remove intermediate files
cd $outdirs
rm -r batches clump