#!/usr/bin/bash

# Eric J. Earley
# 2023
# RTI International

#this code will:
#1. create a phenotype file from pre-existing metadata
#2. run a saigeGDS GWAS on pre-filtered, MAF > 1%, imputed .gds files
#3. combine, filter (Rsq >0.8), and sort into a single .stats file
#4. manhattan plot & qq plot
#


export METABOLITE="Acetoacetate"
export scriptDir=/share/storage/REDS-III/RBCOmics/metabolomics/workflow/codebase
export phenoDir=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes
export phenoOut=$phenoDir/phenotype_${METABOLITE}.txt
export pcFile=$phenoDir/phenotype_prep/PCA/combined.allChr.pca.eigenvec
export metaboliteFile=$phenoDir/phenotype_prep/redsiii_norm_impute_int_20230505.txt
export demographicFile=$phenoDir/phenotype_prep/rbcomics_metadata_for_metabolomics.csv
export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/${METABOLITE}




#### Create phenotype file ####
# This creates a phenotype file in the $phenoDir which is formatted for ProbABEL
if [ -f "$phenoOut" ]; then
  echo "Phenotype file exists - $phenoOut -- skipping this step"
else
  echo "Creating phenotype file: $phenoOut"
  Rscript $scriptDir/01_create_phenotype.R \
    --metabolite $METABOLITE \
    --metaboliteFile $metaboliteFile \
    --pc $pcFile \
    --demographic $demographicFile \
    --out $phenoOut
fi


#### Run saigeGDS ####
# This runs saigeGDS and creates a processingDir if it doesn't existi already
if [ -f "$processingDir/_saige_complete" ]; then
  echo "Found the file $processingDir/_saige_complete -- skipping this step"
else 
  [ ! -d $processingDir ] &&  mkdir $processingDir
  echo "Performing GWAS..."
  sh $scriptDir/02_run_saige.sh


  # pause here to wait for tasks to complete
  status=`qstat | grep rbc`
  sleep 1
  while [ ! -z "$status" ];
  do
      sleep 1
      status=`qstat | grep rbc`
      echo "working..."
  done
  echo "GWAS DONE"




  touch $processingDir/_saige_complete
fi


#### Combine and filter ####
if [ -f "$processingDir/${METABOLITE}_rbc.ALL.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz" ]; then
  echo "Found the combined stats file $processingDir/${METABOLITE}_rbc.ALL.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz -- skipping this step"
else
  echo "combining results into one stats file"
  sh $scriptDir/03_process_saige_results.sh
fi


#### Plotting ####
if [ -f "$processingDir/${METABOLITE}_rbc.ALL.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.plot.snps+indels.manhattan.png" ]; then
  echo "Found plot file $processingDir/${METABOLITE}_rbc.ALL.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz.plot.snps+indels.manhattan.png -- skipping this step"
else
  echo "Plotting..."
  sh $scriptDir/04_plotting.sh
fi

#### Locus Zoom ####
# to be developed



#### Cleanup ####
status=`ls $processingDir | grep -c "assoc.txt"`
#if [ ]; then
echo "Cleaning up the $processingDir"
sh $scriptDir/05_cleanup.sh






#### COJO Conditional GWAS ####
if [ -f "${processingDir}/${METABOLITE}_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.cojo.ma" ]; then
  echo "Found file ${processingDir}/${METABOLITE}_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.cojo.ma -- skipping this step"
else
  echo "Running cojo step-wise conditional regressions"
  sh $scriptDir/06_cojo_conditional_gwas.sh
fi




