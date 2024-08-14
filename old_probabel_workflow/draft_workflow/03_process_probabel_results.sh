
PHENO="ATP"

working_dir=/share/storage/REDS-III/RBCOmics/metabolomics
out_file=${working_dir}/processing/${PHENO}_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats


#### this is already done. The input dosage files were filtered to MAF>=1%
# FIRST, filter to MAF >= 1%
#for (( chr=1; chr<23; chr++ )); do
#  stats_file=${working_dir}/processing/${PHENO}_rbc.all.1000G_p3.chr${chr}.age_gender_hub_donate_BMI_evs_SNP.stats
#  awk '{FS=" "} {if ($5 >= 0.01 && $5 <= 0.99) print $0}' $stats_file > ${working_dir}/processing/chr${chr}_maf01
#done



# Combine all chromosomes into one file
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  stats_file=${working_dir}/processing/${PHENO}_rbc.all.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats

  # these are un-sorted. Sort them by position here
  sort -k2 $stats_file > ${stats_file}.sorted
  

  # add chromosome column and SNP type
  # only allow Rsq > 0.8
  #awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")'
  awk -v var="$chr" 'BEGIN {OFS=" "} {if ($8 > 0.8) {print $0,var}}' ${stats_file}.sorted | awk '{if ( ($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' >> $out_file 
done

# add header to combined fileprint $0,"INDEL"nprint $0,"INDEL"
#name chrom position A1 A2 Freq1 MAF Quality Rsq n Mean_predictor_allele beta_SNP_add se_beta_SNP_add chi2_SNP chi p or_95_percent_ci
echo -e "VARIANT_ID POSITION A1 A2 Freq1 MAF Quality Rsq N Mean_predictor_allele beta_SNP_add se_beta_SNP_add chr2_SNP chi P OR_95_percent_ci CHR TYPE" > $working_dir/processing/header
cat $working_dir/processing/header $out_file > tmp
mv tmp $out_file


# compress it
gzip $out_file


# remove intermediate files
#rm $working_dir/processing/*chr*
rm -r $working_dir/processing/work
