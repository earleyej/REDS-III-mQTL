#### NOTES: ####
# This script will combine the .sig files (created by saigeGDS) and the annotation (created by ENSEMBL Variant Effect Predictor)

# 
setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/results/280/sumstats/ALL/")
# CHANGE THIS WHEN RUNNING A NEW ANNOTATION!
annotFile = "../annotation/vep_raw_results_20240220.txt"
outFile = "../annotation/redsiii_sigSNPs_annot_20240220.txt"





######################
#### Extract rsID ####
###################### 

# pull out the rsID column from each file, sort, and uniq-ify.
# then upload to VEP for annotation
files=list.files(pattern="sig$")
files=files[sapply(files, file.size) > 1] #remove files with no sig results

out=NULL
for (i in 1:length(files)) {
  cat(files[i],"\n")
  sig=read.table(files[i],sep="\t",header=F,stringsAsFactors = F)
  if (dim(sig)[2]>1) {
    sig$rsid=ifelse(grepl("rs",sig$V4),
                  gsub(":.*$","",sig$V4),
                  sig$V1)  
    out=c(out,sig$rsid)
  } else {
    print("failed dim check")
  }
} 

# # Fructose has a different format:
# sig = read.table("DFructose16bisphosphate_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig",
#                  sep="\t", header=F, stringsAsFactors = F)
# sig$rsid=ifelse(grepl("rs",sig$V4),
#                 gsub(":.*$","",sig$V4),
#                 sig$V4)
# out=c(out, sig$rsid)


# uniq-ify
out = unique(sort(out))
out = out[grepl("^rs",out)]
length(out)


write.table(out, file="annotate/FOR_ANNOT.txt",quote=F, row.names=F, col.names=F)



#############################
#### GWAS Catalog Lookup ####
#############################
# Run this jupyter notebook:
# C:\Users\eearley\OneDrive - Research Triangle Institute\Documents\REDS-III\metabolomics\results\glycolysis_probabel\ALL\annotate\Glycolysis_rsID_GWAS_Catalog_Lookup.ipynb

# add metabolite to this file
gwas_cat = read.table("annotate/Glycolysis_GWAS_Catalog_RESULTS_2023_12_07.txt", sep="\t", header=T, stringsAsFactors = F)





files=list.files(pattern="sig")
#FructoseBiphosphate has a different format, so process that one later
files=files[!grepl("Fructose",files)]

out=NULL
for (i in 1:length(files)) {
  stats = read.table(files[i],sep=" ",header=F, stringsAsFactors = F)
 
  #remove extraneous stuff from rsID column
  stats$V1=ifelse(grepl("rs",stats$V1),
                 gsub(":.*$","",stats$V1),
                 stats$V1)


  #subset the stats
  #name pos A1 A2 Freq1 MAF Quality Rsq n Mean_predictor/2 chrom position beta_mu beta_sex beta_age beta_SNP  sebeta_mu sebeta_sex sebeta_age sebeta_SNP sigma2 SNP_Z SNP_chi2

  #stats = stats[,c(1,2,3,4,6,8,11,12,15,17)]
  #colnames(stats) = c("ID","POS","REF","ALT","MAF","RSQ","BETA","SE","P","CHROM")

  #combine
  x = merge(stats,gwas_cat, by.x="V1",by.y="variant")
  x$Metabolite=gsub("_rbc.*$","",files[i])


  #subset and re-order
  colnames(x)[1:18]=c("rsID","POSITION","REF","ALT","FreqREF","MAF","Quality","Rsq","N","Mean_predictor_allele","beta_SNP_add","se_beta_SNP_add","chr2_SNP","chi","P","OR_95_percent","CHR","TYPE")
  colnames(x)[19:20]=c("GWAS_Catalog_Trait","GWAS_Catalog_P")
  x = x[,c("Metabolite","rsID","TYPE","CHR","POSITION",colnames(x)[c(3:16,19,20)])]
  x = x[,c("Metabolite","rsID","TYPE","CHR","POSITION","REF","ALT","MAF","beta_SNP_add","se_beta_SNP_add","P","GWAS_Catalog_Trait","GWAS_Catalog_P")]
  colnames(x)[c(9,10)]=c("BETA","SE")
  out = rbind(out, x)
}
stats = read.table("DFructose16bisphosphate_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig",sep="\t",header=F, stringsAsFactors = F)
stats = stats[,c(2,3,4,5,6,7,10,11,12)]
colnames(stats) = c("CHR","POSITION","rsID","ALT","REF","REF_AF","BETA","SE","P")
stats$rsID = ifelse(grepl("rs",stats$rsID),
                  gsub(":.*$","",stats$rsID),
                  stats$rsID)
x = merge(stats,gwas_cat, by.x="rsID",by.y="variant")
x$Metabolite="Fructose-biphosphate"
colnames(x)[c(10,11)]=c("GWAS_Catalog_Trait","GWAS_Catalog_P")
x$MAF = ifelse(x$REF_AF > 0.5, 1-x$REF_AF, x$REF_AF)
x$TYPE = ifelse(nchar(x$ALT) > 1 | nchar(x$REF) > 1, "INDEL","SNP")
x = x[,colnames(out)]

out = data.frame(rbind(out,x))

write.table(out, file="annotate/glycolysis_GWAS_Catalog_lookup_2023_12_07.txt", sep="\t", col.names=T, row.names = F, quote=F)
# 







##################################
#### Variant Effect Predictor ####
##################################
# Take the for_annot.txt above and upload to:
# https://useast.ensembl.org/Tools/VEP



##########################################
#### Combine VEF Annotation with .sig ####
##########################################

annot = read.table(annotFile,sep="\t",header=T, stringsAsFactors = F)
dim(annot)

files=list.files(pattern="sig")
#remove files with no sig results
files=files[sapply(files, file.size) > 1]
length(files)

out=NULL
for (i in 1:length(files)) {
  cat(files[i],"\n")
  stats = read.table(files[i],sep="\t",header=F, stringsAsFactors = F, fill=T)
  if (dim(stats)[2]>1) {
    #remove extraneous stuff from rsID column
    stats$rsid=ifelse(grepl("rs",stats$V4),
                    gsub(":.*$","",stats$V4),
                    stats$V4)
  
    #subset the stats
    stats = stats[,c(2,3,4,5,6,7,10,11,12,14)]
    stats$V7 = ifelse(stats$V7 > 0.5, 1-stats$V7,stats$V7)
    colnames(stats) = c("CHR","POS","ID","ALT","REF","MAF","BETA","SE","P","rsid")
    Metabolite=gsub("_rbc.*$","",files[i])
    stats$Metabolite = Metabolite
    #add population
    pop = gsub("^.*rbc.","",files[i])
    pop = gsub("\\.1000G.*$","",pop)
    
    stats$pop = pop
    stats = stats[,c("Metabolite","pop","CHR","POS","ID","ALT","REF","MAF","BETA","SE","P","rsid")]
  
    #combine
    x = merge(stats, annot, by.y="Uploaded_variation",by.x="rsid",all.x=T,all.y=F)
  
  
    #remove rows with no annotation
    x = x[!is.na(x$Location),]
  
    #make sure the alleles match
    x = x[x$Allele == x$ALT,]
  
    #update column names
    colnames(x)[4] = "POS_b37"
    colnames(x)[12] = "Location_b38"
  
    out = rbind(out, x)
  } else {
    print("failed dim test")
  }
  
}

length(unique(out$Metabolite))
#table(out$Metabolite, out$pop)

# intron/exon column will be interpreted as date in Excel.
out$INTRON = gsub("/","_of_",out$INTRON)
out$EXON = gsub("/","_of_",out$EXON)

#
write.table(out, file=outFile, sep="\t",col.names=T, row.names=F,quote=F)









# 
# 
# 
# 
# 
# 
# 
# #### Glycolysis paper ####
# 
# setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/results/sumstats/glycolysis/")
# 
# #combine all the .sig files together and write out just the rsIDs with metabolite
# files=list.files(pattern="sig$")
# 
# #FructoseBiphosphate has a different format, so process that one later
# files=files[!grepl("Fructose",files)]
# 
# 
# out=NULL
# i=1
# for (i in 1:length(files)) {
#   cat(files[i])
#   sig=read.table(files[i],sep=" ",header=F,stringsAsFactors = F)
#   sig$V1=ifelse(grepl("rs",sig$V1),
#                   gsub(":.*$","",sig$V1),
#                 sig$V1)  
#   #tmp=data.frame(rsid = sig$V1)
#   #tmp=data.frame("Metabolite"=gsub("_rbc.*$","",files[i]),
#   #               "rsID"=sig$V1)
#   out=c(out,sig$V1)
# }
# 
# # add in the FructoseBiphosphate which is formatted differently
# sig = read.table("DFructose16bisphosphate_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig", sep="\t",header=F, stringsAsFactors = )
# sig$rsid = gsub(":.*$","",sig$V4)
# sig = sig[grepl("^rs",sig$rsid),]
# out = c(out, sig$rsid)
# 
# # uniq-ify
# out = out[!duplicated(out)]
# out = out[grepl("^rs",out)]
# 
# 
# write.table(out, file="annotate/rsIDs_for_annotation_20230614.txt",sep="\t",col.names=F,row.names=F,quote=F)
# 
# # now annotate this with either osasis or the ENSEMBL variant predictor
# # https://useast.ensembl.org/Homo_sapiens/Tools/VEP?db=core;tl=yp5sowVSmZLKDl7o-9257374
# 
# 
# 
# 
# #### MERGING ANNOTATION RESULTS WITH .SIG #
# 
# oasis = read.table("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/results/oasis/oasis_output/OASIS_variants_next_9_annotation.txt", sep="\t", header=T, stringsAsFactors = F, fill=T)
# # ensembl variant predictor
# #annot = read.table("annotate/ensembl_variant_predictor_20230614.txt",sep="\t",header=T, stringsAsFactors = F)
# 
# 
# i=1
# 
# out=NULL
# for (i in 1:length(files)) {
#   #files=list.files(pattern="sig")
#   stats = read.table(files[i],sep=" ",header=F, stringsAsFactors = F)
#   
#   #remove extraneous stuff from rsID column
#   stats$V1=ifelse(grepl("rs",stats$V1),
#                   gsub(":.*$","",stats$V1),
#                   stats$V1)
#   
#   
#   #subset the stats
#   #name pos A1 A2 Freq1 MAF Quality Rsq n Mean_predictor/2 chrom position beta_mu beta_sex beta_age beta_SNP  sebeta_mu sebeta_sex sebeta_age sebeta_SNP sigma2 SNP_Z SNP_chi2
#   
#   #stats = stats[,c(1,2,3,4,6,8,11,12,15,17)]
#   #colnames(stats) = c("ID","POS","REF","ALT","MAF","RSQ","BETA","SE","P","CHROM")
#   
#   #combine
#   x = merge(stats,oasis, by.x="V1",by.y="rsNum")
#   x$Metabolite=gsub("_rbc.*$","",files[i])
#   
#   
#   #subset and re-order
#   #x = x[,c("Metabolite","rsNum","SNPname","Chr37","Pos37","Ref","Alt","V6","V8","V11","V12","V15","Gene","Dist","Type","Function","AAchg","RglmDB","Dnase","Reg","SIFT","PP2_HDIV","LRT","MT","MA","FATHMM","metaSVM","REVEL","GERP.2","CADD","ClinVar","islets","adipose","liver","sklmus")]
#   #colnames(x)[c(7:11)]=c("MAF","Rsq","BETA","SE_BETA","P")  
#   colnames(x)[1:18]=c("VARIANT_ID","POSITION","A1","A2","Freq1","MAF","Quality","Rsq","N","Mean_predictor_allele","beta_SNP_add","se_beta_SNP_add","chr2_SNP","chi","P","OR_95_percent","CHR","TYPE")
#   
#   out = rbind(out, x)
# }
# # checking for na rows
# sum(is.na(out$CHROM))
# # none!
# 
# 
# #
# write.table(out, file="annotate/Glycolysis_redsiii_sigSNPs_annot_2023_0615.txt", sep="\t",col.names=T, row.names=F,quote=F)
# 
# 
# stats = read.table("DFructose16bisphosphate_rbc.all.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig",sep="\t",header=F, stringsAsFactors = F)
# stats = stats[,c(2,3,4,5,6,7,10,11,12)]
# colnames(stats) = c("CHROM","POS","ID","ALT","REF","REF_AF","BETA","SE","P")
# stats$ID = ifelse(grepl("rs",stats$ID),
#                   gsub(":.*$","",stats$ID),
#                   stats$ID)
# x = merge(annot, stats, by.x="Uploaded_variation",by.y="ID")
# write.table(x, file="annotate/FructoseBiphosphate_sigSNPs_annot_2023_0615.txt", sep="\t", col.names=T, row.names = F, quote=F)
# 
