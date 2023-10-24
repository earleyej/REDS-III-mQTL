

# check if the pngs are there and have non-zero file size

echo "Cleaning up intermediate files..."
# remove intermediate files
echo "Removing $processingDir/work/"
rm -r $processingDir/work
echo "Removing per-chrom stats files"
rm $processingDir/${METABOLITE}_rbc.all.1000G_p3.chr*.age_gender_hub_donate_evs_SNP.stats
rm $processingDir/${METABOLITE}_rbc.all.1000G_p3.chr*.age_gender_hub_donate_evs_SNP.stats.sorted
rm $processingDir/header

