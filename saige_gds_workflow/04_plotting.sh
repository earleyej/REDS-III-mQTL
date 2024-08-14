#!/bin/bash

[ -z "$ancestry" ] && export ancestry=ALL
export INPUT_FILE=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz #23Bisphosphoglycerate_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz


# create list of significant SNPs
echo "Marking significant SNPs..."
zcat $INPUT_FILE | awk '{FS=" "} {if ($12 < 0.00000005) print $0}' > ${INPUT_FILE}.sig


# plot
echo "plotting..."
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
    --col_id rs.id \
    --col_chromosome chr \
    --col_position pos \
    --col_p pval \
    --generate_snp_indel_manhattan_plot \
    --manhattan_odd_chr_color red \
    --manhattan_even_chr_color blue \
    --manhattan_points_cex 1.5 \
    --manhattan_cex_axis 2 \
    --manhattan_cex_lab 2 \
    --col_variant_type type \
    --generate_snp_indel_qq_plot \
    --qqlines \
    --qq_points_bg black \
    --qq_lambda 

#    --highlight_list ${INPUT_FILE}.sig \

