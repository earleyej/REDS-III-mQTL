# Eric Earley
# RTI International
# Oct 2023

# This code creates a composite Manhattan plot of significant mQTL hits within the 280 Metabolite RBC project
# These hits came from a meta-analysis of REDS-III cohort (N=13K) population-stratified GWASes

# input: summary stats for all sig SNPs
# output: png manhattan plot

# Note: If I only plot the significant SNPs then the width of chromosomes is messed up.
#   I could specify the boundaries of the chroms manually... that feels like a pain
#   OR I could load in one full sumstat file for one metabolite and add to plot data, then maybe remove from plot?


setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/results/280/")

df = read.table(gzfile("sumstats/annotation/RBCOmics_mQTL_annotated_sig_SNPs_2024_05_29.txt.gz"), sep="\t", header=T, stringsAsFactors = F, quote="")
dim(df)

# save a subset of this table for plotting
plot.me = df[,c("rsid","Metabolite","Chr","Pos_b37","Pvalue","lead.snp","DISEASE.TRAIT", "PUBMEDID")]
# remove SNPs lacking an rsID
plot.me = plot.me[!is.na(plot.me$rsid),]
# plot ALL significant SNPs, but uniq-ify
plot.me = plot.me[with(plot.me, order(Chr, Pos_b37)), ]
dim(plot.me)
plot.me = plot.me[!duplicated(plot.me),]
dim(plot.me)

 
# remove all other non-plottable columns
plot.me.out = plot.me[,c("rsid","Chr","Pos_b37","Pvalue")]
# uniq-ify again
plot.me.out = plot.me.out[!duplicated(plot.me.out),]
# convert Pvalues = 0 to 10^-350
plot.me.out$Pvalue = ifelse(plot.me.out$Pvalue == 0, 1e-350, plot.me.out$Pvalue)

#### Add non-sig SNPs from a random metabolite ####
# this dummy.stats came from metal.AC_10_0_.tbl.annotated.gz
dummy = read.table("plots/dummy.stats.gz", sep="\t", header=T, stringsAsFactors = F)
colnames(dummy)[1] = "rsid"
dummy = dummy[dummy$Pvalue > 0.00000005,] # remove significant hits, since I already collected those in plot.me
plot.me.out = rbind(plot.me.out, dummy)
plot.me.out = plot.me.out[with(plot.me.out, order(Chr, Pos_b37)),]

#### write output ####
write.table(plot.me.out, file="plots/RBCOmics_mQTL_annotated_sig_SNPs_2024_05_29_FORPLOT.txt", sep="\t", col.names=T, row.names=F, quote=F)



#also write a list of rsIDs to highlight in plot. These will be any rsID with a pubmedid of 36395887
prior.hits = plot.me[grepl("36395887",plot.me$PUBMEDID),]
prior.hits.out = unique(prior.hits$rsid)
write.table(prior.hits.out, file="plots/prior.hits", col.names=F, row.names=F, quote=F)






#### run this in bash ####
Rscript ~/eearley/OneDrive\ -\ Research\ Triangle\ Institute/Documents/genomic_resources/generate_gwas_plots.v10.R \
  --in RBCOmics_mQTL_annotated_sig_SNPs_2024_05_29_FORPLOT.txt \
  --in_chromosomes autosomal_nonPAR \
  --in_header \
  --out RBCOmics_mQTL_annotated_sig_SNPs_2024_05_29.plot \
  --col_id rsid \
  --col_chromosome Chr \
  --col_position Pos_b37 \
  --col_p Pvalue \
  --generate_manhattan_plot \
  --manhattan_pch 19 \
  --manhattan_points_cex 1 \
  --highlight_list prior.hits \
  --manhattan_highlight_color yellow \
  --manhattan_odd_chr_color_sig blue \
  --manhattan_even_chr_color_sig red
####














#### OLD ####
# df$pval = NA
# df$pos = NA
# df$chr = NA

# but it also needs p-value for the lead SNPs
files=list.files("results/sumstats/", pattern="sig")
# only keep the 'ALL' cohort
files = files[grepl("ALL",files, ignore.case=T)]

# for each row in this data.frame, add on the pvalue, chrom, and position
for (i in 1:dim(df)[1]) {
  metabolite = df[i,"Metabolite"]
  snp = df[i,"SNP"]
  print(paste0(metabolite, ":", snp))
  this.file = files[grepl(metabolite,files)][1]
  tmp = read.table(gzfile(paste0("results/sumstats/",
                                 this.file)), 
                   sep="\t", header=F, stringsAsFactors = F)  
  tmp$snp = sapply(strsplit(tmp$V4,split=":"), `[`,1)
  pval = tmp[tmp$snp == snp,"V12"]
  pval = ifelse(length(pval) == 0, NA, pval)
  
  pos = ifelse(length(tmp[tmp$snp == snp,"V3"]) == 0, NA, tmp[tmp$snp == snp,"V3"])
  chr = ifelse(length(tmp[tmp$snp == snp,"V2"]) == 0, NA, tmp[tmp$snp == snp,"V2"])
  
    
  df[i,"pval"] = pval  
  df[i,"pos"] = pos
  df[i,"chr"] = chr
}


# save this file
write.table(df, file="results/lead_snps_2023_10_18.txt", sep="\t", col.names=T, row.names=F, quote=F)



# plot in a second script:
#Rscript ~/eearley/OneDrive\ -\ Research\ Triangle\ Institute/Documents/genomic_resources/generate_gwas_plots.v10.R --in lead_snps_2023_10_18.txt --in_chromosomes autosomal_nonPAR --in_header --out lead_snps_2023_10_18.plot --col_id SNP --col_chromosome chr --col_position pos --col_p pval --generate_manhattan_plot --manhattan_pch 19 --manhattan_points_cex 1.5



