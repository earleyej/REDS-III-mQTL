#### format stats file ####
# LocusZoom expects a format like METAL
# MarkerName<\t>P-value

#Find lead SNPs from the .sig file created in the last step
ANALYTE="ATP"
WORKING_DIR=/share/storage/REDS-III/RBCOmics/metabolomics/processing
gwas_stats=$WORKING_DIR/${ANALYTE}_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz
sigFile=${gwas_stats}.sig


# create METAL formatted file that has just the local region surrounding the lead SNP
$scriptDir/tmake_locus_zoom_input.R




#### plotting ####
/share/storage/REDS-III/common/locuszoom/bin/locuszoom \
--source 1000G_Nov2014 --build hg19 \
--metal /share/storage/REDS-III/almoore/metabolomics/final/rbc.all.1000G_p3.D.Glucose.6.phosphate_v09272021_rs76059814_plusminus200kb_forlocuszoom.txt \
--pop ALL \
--refsnp rs76059814 \
--flank 200kb \
--pvalcol P \
--delim '\t' \
--markercol VARIANT_ID \
--prefix D.Glucose.6.phosphate.EUR.rs76059814.locuszoom

