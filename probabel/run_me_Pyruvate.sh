#!/usr/bin/bash

# Eric J. Earley
# 2023
# RTI International

#this code will:
#1. create a phenotype file from pre-existing metadata
#2. run a ProbABEL GWAS using NextFlow (MAF > 1% filtered .gen files)
#3. combine, filter (Rsq >0.8), and sort into a single .stats file
#4. manhattan plot & qq plot
#
# it will also automatically delete all intermediate files, including the absurdly large 'work' directory which is like 500Gb


export METABOLITE="Pyruvate"
export scriptDir=/share/storage/REDS-III/RBCOmics/metabolomics/workflow/codebase
export phenoDir=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes
export phenoOut=$phenoDir/phenotype_${METABOLITE}.txt
export pcFile=$phenoDir/phenotype_prep/PCA/combined.allChr.pca.eigenvec
export metaboliteFile=$phenoDir/phenotype_prep/redsiii_norm_impute_int_20230505.txt
export demographicFile=$phenoDir/phenotype_prep/rbcomics_metadata_for_metabolomics.csv
export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/${METABOLITE}




#### Create phenotype file ####
# This creates a phenotype file in the $phenoDir which is formatted for ProbABEL
echo "Creating phenotype file: $phenoOut"


Rscript $scriptDir/01_create_phenotype.R \
  --metabolite $METABOLITE \
  --metaboliteFile $metaboliteFile \
  --pc $pcFile \
  --demographic $demographicFile \
  --out $phenoOut


#### Run NextFlow ####
# This runs the NextFlow pipeline and creates a 'work' directory wherever it is run. The 'cd $processingDir' is meant to isolate this 'work' dir to the correct parent dir
[ ! -d $processingDir ] &&  mkdir $processingDir
cd $processingDir
sh $scriptDir/02_nextflow_probabel.sh

echo "\n"
echo "##############################################\n"
echo "	GWAS Complete!	\n"
echo "##############################################\n"
echo "\n"

#### Combine and filter results ####
# Since the previous step needs to be run on an interactive compute node, no need to pause here, but this how you'd do it

# check for status of qsub jobs
# when grep finds something, status is 0
#qstat | grep prob &> /dev/null
#sleep 1
#while [ $? -eq 0 ]; 
#do
#    sleep 1
#    qstat | grep prob &> /dev/null
#    #echo "working..."
#done;
#echo "DONE"

sh $scriptDir/03_process_probabel_results.sh



#### Plotting ####
sh $scriptDir/04_plotting.sh



#### Locus Zoom ####
# to be developed


