


export METABOLITE="Lactate"
export scriptDir=/share/storage/REDS-III/RBCOmics/metabolomics/workflow/probabel/codebase
export phenoDir=/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes
export phenoOut=$phenoDir/phenotype_${METABOLITE}.txt
export pcFile=$phenoDir/phenotype_prep/PCA/combined.allChr.pca.eigenvec
export metaboliteFile=$phenoDir/phenotype_prep/redsiii_norm_impute_int_20230505.txt
export demographicFile=$phenoDir/phenotype_prep/rbcomics_metadata_for_metabolomics.csv
export processingDir=/share/storage/REDS-III/RBCOmics/metabolomics/processing/${METABOLITE}/ancestry_stratified
export imputeDir=/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01


#### Create phenotype ####
# just subset the existing pheno file
cut -f1 $imputeDir/chr22/rbc.AFRAMRCN.1000G_p3.chr22.0.fam > $processingDir/samples_AFAM.txt
cut -f1 $imputeDir/chr22/rbc.ASIAN.1000G_p3.chr22.0.fam > $processingDir/samples_ASIAN.txt
cut -f1 $imputeDir/chr22/rbc.CAUCASIAN.1000G_p3.chr22.0.fam > $processingDir/samples_CAUCASIAN.txt
cut -f1 $imputeDir/chr22/rbc.HISPANIC.1000G_p3.chr22.0.fam > $processingDir/samples_HISPANIC.txt

head -n 1 $phenoOut > $processingDir/head
grep -f $processingDir/samples_AFAM.txt $phenoOut > $processingDir/bottom; cat $processingDir/head $processingDir/bottom > $processingDir/pheno_AFAM.txt
grep -f $processingDir/samples_ASIAN.txt $phenoOut > $processingDir/bottom; cat $processingDir/head $processingDir/bottom > $processingDir/pheno_ASIAN.txt
grep -f $processingDir/samples_CAUCASIAN.txt $phenoOut > $processingDir/bottom; cat $processingDir/head $processingDir/bottom > $processingDir/pheno_CAUCASIAN.txt
grep -f $processingDir/samples_HISPANIC.txt $phenoOut > $processingDir/bottom; cat $processingDir/head $processingDir/bottom > $processingDir/pheno_HISPANIC.txt


rm head bottom



#### run Probabel GWAS ####

#export imputation_root=/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_by_race
export method=palinear


##############
#### AFAM ####
##############
export ancestry=AFAM
export phenoFile=$processingDir/pheno_${ancestry}.txt
for (( chr=1; chr<24; chr++ )); do
        out_file=$processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
        geno_prefix=rbc.AFRAMRCN.1000G_p3.chr
        /share/storage/REDS-III/common/software/nextflow/nextflow-0.25.1-all \
          /share/storage/REDS-III/common/software/pipelines/_pipeline.association.out_stats_files.v0.1.nf \
          --final_chunks $imputeDir/chunks/final_chunks.chr$chr \
          --input_pheno $phenoFile \
          --imputation_dir $imputeDir/chr$chr \
          --example_mldose $imputeDir/chr$chr/rbc.AFRAMRCN.1000G_p3.chr${chr}.0.mach.mldose.gz \
          --geno_prefix $geno_prefix \
          --working_dirs $processingDir \
          --out $out_file \
          --method $method
done

# WAIT FOR THIS PRIOR STEP TO COMPLETE
out_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  stats_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
  # these are un-sorted. Sort them by position here
  sort -k2 -n $stats_file > ${stats_file}.sorted
  # add chromosome column and SNP type
  # only allow Rsq > 0.8
  #awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")'
  awk -v var="$chr" 'BEGIN {OFS=" "} {if ($8 > 0.8) {print $0,var}}' ${stats_file}.sorted | awk '{if ( ($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' >> $out_file
done
# add header to combined fileprint $0,"INDEL"nprint $0,"INDEL"
#name chrom position A1 A2 Freq1 MAF Quality Rsq n Mean_predictor_allele beta_SNP_add se_beta_SNP_add chi2_SNP chi p or_95_percent_ci
echo -e "VARIANT_ID POSITION A1 A2 Freq1 MAF Quality Rsq N Mean_predictor_allele beta_SNP_add se_beta_SNP_add chr2_SNP chi P OR_95_percent_ci CHR TYPE" > ${processingDir}/header
cat ${processingDir}/header ${out_file} > ${processingDir}/tmp
mv ${processingDir}/tmp $out_file
echo -e "compressing\n"
gzip $out_file
echo -e "subsetting\n"
portable_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats.forFUMA.gz
zcat ${out_file}.gz | awk '{FS=" "} {print $1,$17,$2,$3,$4,$6,$8,$11,$12,$15}' | tr " " "\t" | gzip > $portable_file



export INPUT_FILE=${out_file}.gz
zcat $INPUT_FILE | awk '$15 < 0.00000005' > ${INPUT_FILE}.sig
/home/eearley/apps/qsub_job.sh \
  --job_name plot_${METABOLITE}_${ancestry} \
  --script_prefix $processingDir/plot_${METABOLITE}_${ancestry} \
  --mem 8 \
  --cpu 1 \
  --program /share/apps/R-4.0.3/bin/Rscript /share/storage/REDS-III/common/software/R/generate_gwas_plots.v10.R \
    --in $INPUT_FILE \
    --in_chromosomes autosomal_nonPAR \
    --in_header \
    --out ${INPUT_FILE}.plot \
    --col_id VARIANT_ID \
    --col_chromosome CHR \
    --col_position POSITION \
    --col_p P \
    --generate_snp_indel_manhattan_plot \
    --manhattan_odd_chr_color red \
    --manhattan_even_chr_color blue \
    --manhattan_points_cex 1.5 \
    --manhattan_cex_axis 2 \
    --manhattan_cex_lab 2 \
    --highlight_list ${INPUT_FILE}.sig \
    --col_variant_type TYPE \
    --generate_snp_indel_qq_plot \
    --qqlines \
    --qq_points_bg black \
    --qq_lambda





#############################################
####             ASIAN                   ####
#############################################
export ancestry=ASIAN
export phenoFile=$processingDir/pheno_${ancestry}.txt
for (( chr=1; chr<24; chr++ )); do
        out_file=$processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
        geno_prefix=rbc.ASIAN.1000G_p3.chr
        /share/storage/REDS-III/common/software/nextflow/nextflow-0.25.1-all \
          /share/storage/REDS-III/common/software/pipelines/_pipeline.association.out_stats_files.v0.1.nf \
          --final_chunks $imputeDir/chunks/final_chunks.chr$chr \
          --input_pheno $phenoFile \
          --imputation_dir $imputeDir/chr$chr \
          --example_mldose $imputeDir/chr$chr/rbc.ASIAN.1000G_p3.chr${chr}.0.mach.mldose.gz \
          --geno_prefix $geno_prefix \
          --working_dirs $processingDir \
          --out $out_file \
          --method $method
done

export out_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  stats_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
  # these are un-sorted. Sort them by position here
  sort -k2 -n $stats_file > ${stats_file}.sorted
  # add chromosome column and SNP type
  # only allow Rsq > 0.8
  #awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")'
  awk -v var="$chr" 'BEGIN {OFS=" "} {if ($8 > 0.8) {print $0,var}}' ${stats_file}.sorted | awk '{if ( ($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' >> $out_file
done
# add header to combined fileprint $0,"INDEL"nprint $0,"INDEL"
#name chrom position A1 A2 Freq1 MAF Quality Rsq n Mean_predictor_allele beta_SNP_add se_beta_SNP_add chi2_SNP chi p or_95_percent_ci
echo -e "VARIANT_ID POSITION A1 A2 Freq1 MAF Quality Rsq N Mean_predictor_allele beta_SNP_add se_beta_SNP_add chr2_SNP chi P OR_95_percent_ci CHR TYPE" > ${processingDir}/header
cat ${processingDir}/header ${out_file} > ${processingDir}/tmp
mv ${processingDir}/tmp $out_file
echo -e "compressing\n"
gzip $out_file
echo -e "subsetting\n"
portable_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats.forFUMA.gz
zcat ${out_file}.gz | awk '{FS=" "} {print $1,$17,$2,$3,$4,$6,$8,$11,$12,$15}' | tr " " "\t" | gzip > $portable_file

export INPUT_FILE=${out_file}.gz
zcat $INPUT_FILE | awk '$15 < 0.00000005' > ${INPUT_FILE}.sig
/home/eearley/apps/qsub_job.sh \
  --job_name plot_${METABOLITE}_${ancestry} \
  --script_prefix $processingDir/plot_${METABOLITE}_${ancestry} \
  --mem 8 \
  --cpu 1 \
  --program /share/apps/R-4.0.3/bin/Rscript /share/storage/REDS-III/common/software/R/generate_gwas_plots.v10.R \
    --in $INPUT_FILE \
    --in_chromosomes autosomal_nonPAR \
    --in_header \
    --out ${INPUT_FILE}.plot \
    --col_id VARIANT_ID \
    --col_chromosome CHR \
    --col_position POSITION \
    --col_p P \
    --generate_snp_indel_manhattan_plot \
    --manhattan_odd_chr_color red \
    --manhattan_even_chr_color blue \
    --manhattan_points_cex 1.5 \
    --manhattan_cex_axis 2 \
    --manhattan_cex_lab 2 \
    --highlight_list ${INPUT_FILE}.sig \
    --col_variant_type TYPE \
    --generate_snp_indel_qq_plot \
    --qqlines \
    --qq_points_bg black \
    --qq_lambda












###################
#### CAUCASIAN ####
###################
export ancestry=CAUCASIAN
export phenoFile=$processingDir/pheno_${ancestry}.txt
for (( chr=1; chr<24; chr++ )); do
        out_file=$processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
        geno_prefix=rbc.CAUCASIAN.1000G_p3.chr
        /share/storage/REDS-III/common/software/nextflow/nextflow-0.25.1-all \
          /share/storage/REDS-III/common/software/pipelines/_pipeline.association.out_stats_files.v0.1.nf \
          --final_chunks $imputeDir/chunks/final_chunks.chr$chr \
          --input_pheno $phenoFile \
          --imputation_dir $imputeDir/chr$chr \
          --example_mldose $imputeDir/chr$chr/rbc.CAUCASIAN.1000G_p3.chr${chr}.0.mach.mldose.gz \
          --geno_prefix $geno_prefix \
          --working_dirs $processingDir \
          --out $out_file \
          --method $method
  done
out_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  stats_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats
  # these are un-sorted. Sort them by position here
  sort -k2 -n $stats_file > ${stats_file}.sorted
  # add chromosome column and SNP type
  # only allow Rsq > 0.8
  #awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")'
  awk -v var="$chr" 'BEGIN {OFS=" "} {if ($8 > 0.8) {print $0,var}}' ${stats_file}.sorted | awk '{if ( ($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' >> $out_file
done
# add header to combined fileprint $0,"INDEL"nprint $0,"INDEL"
#name chrom position A1 A2 Freq1 MAF Quality Rsq n Mean_predictor_allele beta_SNP_add se_beta_SNP_add chi2_SNP chi p or_95_percent_ci
echo -e "VARIANT_ID POSITION A1 A2 Freq1 MAF Quality Rsq N Mean_predictor_allele beta_SNP_add se_beta_SNP_add chr2_SNP chi P OR_95_percent_ci CHR TYPE" > ${processingDir}/header
cat ${processingDir}/header ${out_file} > ${processingDir}/tmp
mv ${processingDir}/tmp $out_file
echo -e "compressing\n"
gzip $out_file
echo -e "subsetting\n"
portable_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_evs_SNP.stats.forFUMA.gz
zcat ${out_file}.gz | awk '{FS=" "} {print $1,$17,$2,$3,$4,$6,$8,$11,$12,$15}' | tr " " "\t" | gzip > $portable_file








