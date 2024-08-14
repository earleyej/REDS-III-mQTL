/share/storage/REDS-III/common/locuszoom/bin/locuszoom \
--source 1000G_Nov2014 --build hg19 \
--metal /share/storage/REDS-III/almoore/metabolomics/final/rbc.all.1000G_p3.LPS18.0_v09272021_rs73264680_plusminus200kb_forlocuszoom.txt \
--pop EUR \
--refsnp rs73264680 \
--flank 200kb \
--pvalcol P-value \
--delim '\t' \
--markercol MarkerName \
--prefix LPS18.0.EUR.rs73264680.locuszoom 