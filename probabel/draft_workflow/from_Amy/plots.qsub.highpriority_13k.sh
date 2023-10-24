#!/bin/bash 
cd /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/
for dir2 in *; do 
#dir2="kynurenine_13k_day42_narm"

echo "hello world"
#pulls out significant rows
cat /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir2"/rbc.all.1000G_p3.allchr.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table | awk '$6<0.00000005' > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir2"/rbc.all.1000G_p3.allchr.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table.GWASsig

#sed -i "s/$/\t$dir2/" /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/x1/"$dir2"/v09272021.rbc.all.1000G_p3.allchr.${dir2}~age+gender+hub+donate+10evs+SNP.stats_maf_gt_0.05_rsq_gt_0.90_noindels.table.GWASsig

#cat /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/x1/"$dir2"/v09272021.rbc.all.1000G_p3.allchr.${dir2}~age+gender+hub+donate+10evs+SNP.stats_maf_gt_0.05_rsq_gt_0.90_noindels.table.GWASsig >> /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/x1/v09272021.rbc.all.1000G_p3.allchr.x1~age+gender+hub+donate+10evs+SNP.stats_maf_gt_0.05_rsq_gt_0.90_noindels.table.GWASsig

/share/apps/R-4.0.3/bin/Rscript /share/storage/REDS-III/common/software/R/generate_gwas_plots.v9.R --in /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir2"/rbc.all.1000G_p3.allchr.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table --in_chromosomes autosomal_nonPAR --in_header --out /share/storage/REDS-III/almoore/metabolomics/final/highpriority/rbc.all.1000G_p3.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.maf_gt_0.01_rsq_gt_0.80_noindels --col_id VARIANT_ID --col_chromosome CHR --col_position POSITION --col_p P --col_variant_type TYPE --generate_snp_indel_manhattan_plot --manhattan_odd_chr_color red --manhattan_even_chr_color blue --manhattan_points_cex 1.5 --manhattan_cex_axis 2 --manhattan_cex_lab 2 --generate_snp_indel_qq_plot --qqlines --qq_points_bg black --qq_lambda 

gzip /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir2"/rbc.all.1000G_p3.allchr.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.011_rsq_gt_0.80_noindels.table 
rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir2"/rbc.all.1000G_p3.${dir2}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table 


done 
