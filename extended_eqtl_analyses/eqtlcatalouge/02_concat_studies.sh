#!/bin/sh

codedir="extended_eqtl_analyses/eqtlcatalouge"
outdirs="extended_eqtl_analyses/eqtlcatalouge"
datafile=$outdirs/"studies.txt"

#loop over eQTL types
#outputs removed eventually as intermediate files

for expr_type in {"ge","exon"}
do

outfile=$outdirs/$expr_type".clumped.assoc.txt"


#loop over studies
iter=1

echo $expr_type": "$iter

line=`head -n $iter $datafile | tail -n 1`
study=`echo $line | awk '{print $1}'`
type=`echo $line | awk '{print $2}'`
study_type=$study"_"$type
infile=$outdirs/$study_type/$expr_type".clumped.assoc.txt"
cat $infile > $outfile

for iter in {2..105}
do

echo $expr_type": "$iter

line=`head -n $iter $datafile | tail -n 1`
study=`echo $line | awk '{print $1}'`
type=`echo $line | awk '{print $2}'`
study_type=$study"_"$type
infile=$outdirs/$study_type/$expr_type".clumped.assoc.txt"
cat $infile | tail -n +2 >> $outfile

done

done