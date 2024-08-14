#!/bin/sh 
cd /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/
for dir1 in *; do
#dir1="kynurenine_13k_day42_narm"

cd /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/ 

for CHR in {1..23}
do 
cat /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats | awk '$6>0.01' | awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")' > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels

awk '{ print $1,$2,$3,$4,$5,$8,$9,$11,$12,$15 }' /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns

awk 'NR>0{$(NF+1)="snp"}1' /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp 

awk -v c=${CHR} 'NR>0{$(NF+1)=c}1' /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp_wchr

awk 'BEGIN {FS=" "; OFS="\t"} {print $1,$12,$2,$8,$9,$10,$3,$5,$7,$6,$11 }' /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp_wchr > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table

cat /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr${CHR}.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table >> /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table 

done


sed -e '1i\VARIANT_ID	CHR	POSITION	BETA	SE	P	ALLELE1	FREQA1	N	RSQ	TYPE' /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table > /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.allchr.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table

gzip /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats

rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp_wchr

rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns_wsnp 

rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels_reducedcolumns

rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels.table

rm /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/highpriority/"$dir1"/rbc.all.1000G_p3.chr*.${dir1}~Age+gender+bloodcenter+numDon2yrs+10evs+SNP.stats_maf_gt_0.01_rsq_gt_0.80_noindels

done
