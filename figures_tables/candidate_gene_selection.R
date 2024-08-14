# This code will take lead SNPs from the mQTL project and generate a list of gene candidates
#   2 ways to identify gene candidates
#       1. GPT-4 (after: https://www.medrxiv.org/content/10.1101/2024.05.30.24308179v1.full.pdf)
#       2. proximity (after: https://www.nature.com/articles/s41586-024-07148-y#Sec11)

library(dplyr)
library(stringr)
setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/results/280/")


# load data #
snps = read.table("sumstats/annotation/RBCOmics_mQTL_annotated_sig_SNPs_2024_05_29.txt.gz", sep="\t", header=T, stringsAsFactors = F, quote="")


#### VEP approach ####
# VEP provides gene names for all annotations. Can I just use this information?
df = unique(snps[snps$lead.snp == "lead",c("rsid","Metabolite","Chr","Pos_b37","IMPACT","SYMBOL","BIOTYPE")])
dim(unique(df[,c("Metabolite","rsid")])) # 1087 lead SNPs
#write.table(out, file="sumstats/annotation/RBCOmics_mQTL_annotated_lead_SNPs_VEP_review.txt", sep="\t", col.names=T, row.names=F, quote=F)

# restrict to  protein-coding only
df = df[df$BIOTYPE == "protein_coding",]
dim(unique(df[,c("Metabolite","rsid")])) # 916 lead SNPs with protein_coding annotations
df = unique(df)
out = df %>%
  group_by(Metabolite, rsid) %>%
  mutate(genes = paste0(SYMBOL, collapse = ",")) %>%
  mutate(impacts = paste0(IMPACT, collapse = ","))
out = subset(out, select = -c(IMPACT, SYMBOL))
out = unique(out)
dim(out)

write.table(out, file="sumstats/annotation/RBCOmics_mQTL_annotated_lead_SNPs_VEP_review.txt", sep="\t", col.names=T, row.names=F, quote=F)

#




#### GPT-4 approach ####
# get the most recent gencode gene annotation
# load sig SNPs and isolate lead SNPs
# get all genes within 500kb (might use bedtools here)
# write out the file and upload to GPT-4


# Gencode annotation
# Did this in bash
# wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.annotation.gff3.gz
# zcat gencode.v46.annotation.gff3.gz | grep -v "^#" | grep -w "gene" | grep "protein_coding" | cut -f1,4,5,9 | awk '{OFS="\t"}{print $1,$2-1,$3,$4}' > gencode.v46.protein_coding.bed


# # load sig SNPs and filter for lead SNPs only
# # convert to bed and write out
# flank = 250000
# snps = snps[snps$lead.snp == "lead",]
# dim(snps)
# snps.out = snps[,c("Chr","Pos_b37","rsid","Metabolite")]
# snps.out$start = snps.out$Pos_b37 - flank
# snps.out$start = ifelse(snps.out$start < 0, 0, snps.out$start)
# snps.out$stop = snps.out$Pos_b37 + flank
# 
# snps.out = snps.out[,c("Chr","start","stop","Metabolite","rsid")]
# # uniq-ify
# snps.out = snps.out[!duplicated(snps.out),]
# write.table(snps.out, file="sumstats/annotation/lead_snps_2024_06_05.bed", sep="\t", col.names=F, row.names=F, quote=F)
# 
# 
# 
# # Bedtools to find intersection
# # had to rename the chrom column in gencode to remove "chr" prefix
# #bedtools intersect -a lead_snps_2024_06_05.bed -b gencode.v46.protein_coding.bed -wa -wb > lead_snps_overlap_gencode_500kb_2024_06_05.txt
# 
# # process the output to look like this:
# # SNP1  Metabolite1  gene1, gene2, gene3, gene4
# # SNP2  Metabolite1  gene1, gene2
# # SNP1  Metabolite2 gene1, gene2
# # SNP2  Metabolite2 gene1, gene2, gene3
# intersect = read.table("sumstats/annotation/lead_snps_overlap_gencode_500kb_2024_06_05.txt", sep="\t", header=F, stringsAsFactors = F, quote="")
# colnames(intersect) = c("Chr","Start","Stop","Metabolite", "rsID","gencode_chr","gencode_start","gencode_stop","details")
# df = intersect[,c("Metabolite","rsID","details")]
# # extract out the "gene_name" info
# tmp = gsub("^.+?;gene_name=","",df$details)
# df$gene = gsub(";.+$","",tmp)
# df = subset(df, select = -c(details))
# df = df %>%
#   group_by(Metabolite, rsID) %>%
#   mutate(genes = paste0(gene, collapse = ",")) 
# df = subset(df, select = -c(gene))
# df = df[!duplicated(df),]
# dim(df)
# head(df)
# df = df[with(df, order(Metabolite, rsID)),]
# 
# #Write this out and upload to GPT4
# write.table(df, file="sumstats/annotation/for_gpt_metabolite_pos_genes.txt", sep="\t", col.names=T, row.names=F, quote=F)
# 
# 
# #
# 




