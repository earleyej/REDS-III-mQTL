# Overlap of mQTL results with prior hits

sig = read.csv(gzfile("results/280/sumstats/annotation/RBCOmics_mQTL_annotated_sig_SNPs_2024_06_25.txt.gz"), sep="\t", header=T)


df = sig[,c("rsid","Metabolite","Chr","Pos_b37","Effect","Pvalue","lead.snp","DISEASE.TRAIT", "PUBMEDID")]

#table(df$PUBMEDID)

overlap = df[!is.na(df$PUBMEDID),]

length(unique(overlap$rsid))
length(unique(sig$rsid))

overlap = unique(overlap)
dim(overlap)


write.table(overlap, "manuscript/figures/supplemental_table_4.txt", sep="\t", col.names=T, row.names=F, quote=F)



head(sort(table(overlap$DISEASE.TRAIT), decreasing=T),50)

# 36395887 - Amy's pilot study
table(overlap[overlap$PUBMEDID == 36395887,"DISEASE.TRAIT"])
