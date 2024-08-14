# Supplemental Table 3 - List of lead SNPs + annotation
setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/")

df = read.csv("results/280/sumstats/annotation/RBCOmics_mQTL_annotated_lead_SNPs_VEP_2024_07_29.txt", sep="\t", header=T)


# this table is missing alleles and other information. Check the table of all sig. hits
sig = read.csv(gzfile("results/280/sumstats/annotation/RBCOmics_mQTL_annotated_sig_SNPs_2024_06_25.txt.gz"), sep="\t", header=T)
# reduce the size a bit
tmp = unique(sig[,c("rsid","Metabolite","Allele1","Allele2","Effect","StdErr","Pvalue",
                    "Direction","HetPVal","Consequence","SYMBOL","Gene","EXON","INTRON","Protein_position","Amino_acids","Codons")])

df = merge(df, 
           tmp, 
           by=c("rsid","Metabolite"))
df$EXON = paste0("exon_",df$EXON)
df$INTRON = paste0("intron_",df$INTRON)
write.table(df, file="manuscript/figures/supplemental_table_3.txt", sep="\t", col.names=T, row.names=F, quote=F)
