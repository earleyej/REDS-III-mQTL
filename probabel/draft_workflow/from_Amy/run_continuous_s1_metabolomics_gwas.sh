#!/bin/bash 

cd /share/storage/REDS-III/almoore/metabolomics/input/phenotypes/continuous/s1/
for file in *.txt
do 
dir="${file%.txt}"
mkdir -- /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/s1/"$dir"
cd /share/storage/REDS-III/almoore/metabolomics/scripts/
./_pipeline.association.out_stats_files.all.s1.continuous_v09272021.run.sh $dir 
#./statstotable_continuous.sh $dir 
#./plots.qsub_continuous.sh $dir 
done 
