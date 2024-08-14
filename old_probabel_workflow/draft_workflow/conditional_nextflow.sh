
#export phenoDir=$processingDir
export imputation_root=/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01

#export phenotype_file=$working_dir/phenotypes/$PHENO_FILE
export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/DHexosephosphate/conditional_rs144986629
export METABOLITE=DHexosephosphate
export method=palinear
export phenoFile=$processingDir/phenotype_DHexosephosphate_rs144986629.txt

#/share/storage/REDS-III/RBCOmics/metabolomics/processing/DHexosephosphate/conditional_rs144986629


#for (( chr=1; chr<24; chr++ )); do
chr=17
        out_file=$processingDir/${METABOLITE}_rbc.all.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
        geno_prefix=rbc.ALL.1000G_p3.chr

                /share/storage/REDS-III/common/software/nextflow/nextflow-0.25.1-all \
                /share/storage/REDS-III/common/software/pipelines/_pipeline.association.out_stats_files.v0.1.nf \
                --final_chunks $imputation_root/chunks/final_chunks.chr$chr \
                --input_pheno $phenoFile \
                --imputation_dir $imputation_root/chr$chr \
                --example_mldose $imputation_root/chr$chr/rbc.ALL.1000G_p3.chr${chr}.0.mach.mldose.gz \
                --geno_prefix $geno_prefix \
                --working_dirs $processingDir \
                --out $out_file \
                --method $method
#done


