
#PHENO="ATP"

#working_dir=/share/storage/REDS-III/RBCOmics/metabolomics

out_file=${processingDir}/${METABOLITE}_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats


# Combine all chromosomes into one file
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  stats_file=${processingDir}/${METABOLITE}_rbc.all.1000G_p3.chr${chr}.age_gender_hub_donate_evs_SNP.stats

  # these are un-sorted. Sort them by position here
  sort -k2 -n $stats_file > ${stats_file}.sorted
  

  # add chromosome column and SNP type
  # only allow Rsq > 0.8
  #awk '$8>0.80' | awk '($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")'
  awk -v var="$chr" 'BEGIN {OFS=" "} {if ($8 > 0.8) {print $0,var}}' ${stats_file}.sorted | awk '{if ( ($3=="A" || $3=="C" || $3=="G" || $3=="T") && ($4=="A" || $4=="C" || $4=="G" || $4=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' >> $out_file 
done

echo -e "creating combined stats file\n"

# add header to combined fileprint $0,"INDEL"nprint $0,"INDEL"
#name chrom position A1 A2 Freq1 MAF Quality Rsq n Mean_predictor_allele beta_SNP_add se_beta_SNP_add chi2_SNP chi p or_95_percent_ci
echo -e "VARIANT_ID POSITION A1 A2 Freq1 MAF Quality Rsq N Mean_predictor_allele beta_SNP_add se_beta_SNP_add chr2_SNP chi P OR_95_percent_ci CHR TYPE" > ${processingDir}/header
cat ${processingDir}/header ${out_file} > ${processingDir}/tmp
mv ${processingDir}/tmp $out_file

echo -e "compressing\n"
# compress it
gzip $out_file

# subset for a portable version
echo -e "subsetting\n"
portable_file=${processingDir}/${METABOLITE}_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.forFUMA.gz
zcat ${out_file}.gz | awk '{FS=" "} {print $1,$17,$2,$3,$4,$6,$8,$11,$12,$15}' | tr " " "\t" | gzip > $portable_file

# calculate lambda here so it's easier to collect
#round(qchisq(median(pValues[,1]),chiDf,lower.tail=FALSE)/qchisq(0.5,chiDf,lower.tail=FALSE),3)

# sym link to the results directory so it's easier to find later
#ln -s $processingDir/$out_file /share/storage/REDS-III/RBCOmics/metabolomics/results/
ln -s $portable_file /share/storage/REDS-III/RBCOmics/metabolomics/results/

echo -e "Stats file ready: $out_file \n"


#echo "Cleaning up intermediate files..."
# remove intermediate files
#rm $working_dir/processing/*chr*
#rm -r $processingDir/work
#rm $processingDir/${METABOLITE}_rbc.all.1000G_p3.chr*.age_gender_hub_donate_evs_SNP.stats
#rm $processingDir/${METABOLITE}_rbc.all.1000G_p3.chr*.age_gender_hub_donate_evs_SNP.stats.sorted
#rm $processingDir/header

