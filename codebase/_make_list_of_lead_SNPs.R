#library(data.table)
library(dplyr)

sigFile = "ProstaglandinG2_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig"

sig = read.table(sigFile, sep="\t", header=F, stringsAsFactors=F)

colnames(sig)[c(2,3,12)]=c("chr","position","pval")

# this doesn't work
#sig %>% group_by(chr) %>% mutate(start = position - 200000, end = position + 200000) %>% filter(pval == min(pval))

# this will give the lowest p-value by chrom
tmp = sig %>% group_by(chr) %>% filter(pval == min(pval))

# write out the list of rsids
write.table(tmp$V4, file="lead_snps", col.names=F, row.names=F, quote=F)


