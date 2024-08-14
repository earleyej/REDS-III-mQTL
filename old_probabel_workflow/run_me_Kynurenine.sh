export METABOLITE="Kynurenine"
export scriptDir=/share/storage/REDS-III/RBCOmics/metabolomics/workflow/codebase
export phenoDir=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes
export phenoOut=$phenoDir/phenotype_${METABOLITE}.txt
export pcFile=$phenoDir/phenotype_prep/PCA/combined.allChr.pca.eigenvec
export metaboliteFile=$phenoDir/phenotype_prep/redsiii_norm_impute_int_20230501.txt
export demographicFile=$phenoDir/phenotype_prep/rbcomics_metadata_for_metabolomics.csv
export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/${METABOLITE}




#### Create phenotype file ####
# This creates a phenotype file in the $phenoDir which is formatted for ProbABEL


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


#### Combine and filter results ####




#### Plotting ####




#### Locus Zoom ####
