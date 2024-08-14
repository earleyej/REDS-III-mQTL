#!/bin/bash

INPUT_FILE=${processingDir}"/"$METABOLITE"_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz"


# create list of significant SNPs
echo "Marking significant SNPs..."
zcat $INPUT_FILE | awk '{FS=" "} {if ($15 < 0.00000005) print $0}' > ${INPUT_FILE}.sig


# plot
echo "plotting..."
/home/eearley/apps/qsub_job.sh \
  --job_name plot_${METABOLITE} \
  --script_prefix $processingDir/plot_${METABOLITE} \
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


