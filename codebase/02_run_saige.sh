# These variables are set in the environment within the config file
# METABOLITE 

[ -z "$imputeDir" ] && export imputeDir=/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01
#export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/$METABOLITE
#export pheno=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/phenotype_${METABOLITE}.txt


#if ancestry is undefined, then make it "ALL"
[ -z "$ancestry" ] && export ancestry=ALL
[ -z "$pheno" ] && export pheno=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/phenotype_${METABOLITE}.txt

echo "Imputed genotypes: $imputeDir"
echo "Phenotype: $pheno"
echo "Ancestry: $ancestry"


for (( chr=1; chr<24; chr++ )); do
  for chunk in $(grep -P "^$chr\s+" /share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_by_race/final_chunks | perl -lane 'print $F[1];'); do
    export baseName=rbc.${ancestry}.1000G_p3.chr$chr.$chunk #rbc.ALL.1000G_p3.chr22.0.gds
    /home/eearley/apps/qsub_job.sh \
      --job_name ${baseName}.${chr}_${chunk} \
      --script_prefix ${processingDir}/${baseName}.${chr}.${chunk} \
      --mem 8 \
      --priority 0 \
      --cpu 1 \
      --program /home/eearley/apps/saigeGDS.R \
      --geno ${imputeDir}/chr${chr}/${baseName}.gds \
      --pheno $pheno \
      --grm ${imputeDir}/chr${chr}/${baseName}.grm_geno.gds \
      --out $processingDir/${baseName}.assoc.txt
  done
done

