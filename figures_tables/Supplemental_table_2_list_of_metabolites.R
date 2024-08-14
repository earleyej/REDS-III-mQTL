# Supplemental Table 2 - List of metabolites and annotate pathway
setwd("C:/Users/eearley/OneDrive - Research Triangle Institute/Documents/REDS-III/metabolomics/")

#%%%%%%%%%%%%%%%%%%%%%%%%%
#### Annotate Pathway ####
#%%%%%%%%%%%%%%%%%%%%%%%%%

# open the list of metabolites with hits - these have annotated pathways
# hits = read.table("data/metabolite_list_with_hits.txt", sep="\t", header=T, stringsAsFactors = F)
# 
# # list of the 373 metabolites
# complete = read.table("data/list_of_metabolite_full_373.txt")
# 
# # merge them
# df = merge(complete, hits, by.x="V1", by.y="Metabolite",all=T)
# 
# # open exogenous list
# exo = read.table("data/exogenous_compounds.txt")
# df$Pathway = ifelse(df$V1 %in% exo$V1, "Exogenous", df$Pathway)
# sort(table(df$Pathway))
# 
# colnames(df)[1] = "Metabolite"
# 
# # all of the metabolites with no sig. hits are lacking a pathway name
# df[is.na(df$Pathway),]
# dim(df[is.na(df$Pathway),])
# # AC_XX are Carnitine
# df[grepl("^AC",df$Metabolite),"Pathway"] = "Carnitine"
# # FA_XX are Fatty acids
# df[grepl("^FA_",df$Metabolite),"Pathway"] = "Fatty acids"
# 
# # This will be faster by hand... I'll ask Angelo to help complete the pathway IDs

# done - there are some missing, but here is a start:
df = read.csv("data/mQTL_metabolite_list_with_pathways.txt", sep="\t", header=T)

#%%%%%%%%%%%%%%%%%%%
#### Add Lambda ####
#%%%%%%%%%%%%%%%%%%%

# I calculated lambda on BDCat
lambda = read.table("results/280/sumstats/lambda_table.txt", sep="\t", header=T, stringsAsFactors = F)
# merge with df
df = merge(df, lambda, by = "Metabolite", all = T)



#%%%%%%%%%%%%%%%%%%%%%%%%%
#### Add No. Sig hits ####
#%%%%%%%%%%%%%%%%%%%%%%%%%
# I downloaded all 280 .sig files. If they have o bytes then there were no hits
df$hits = ifelse(is.na(df$Lambda), NA, 0)

path = "results/280/sumstats/sig/"
sig.files = list.files(path, pattern="sig$")
length(sig.files)
for (file in sig.files) {
  metabolite = gsub("metal\\.","",file)
  metabolite = gsub("\\.tbl.+$","", metabolite)
  print(metabolite)
  if (file.size(paste0(path,file)) == 0) {
    df[df$Metabolite == metabolite,"hits"] = 0
  } else {
    tmp = read.csv(paste0(path,file), sep="\t", header=F)
    df[df$Metabolite == metabolite,"hits"] = dim(tmp)[1]
  }
}


#%%%%%%%%%%%%%%%%%%%%%%%%%%
#### Add No. lead SNPs ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%


lead = read.table("results/280/sumstats/annotation/RBCOmics_mQTL_annotated_lead_SNPs_VEP_2024_07_16.txt", sep="\t", header=T, stringsAsFactors = F)
lead = unique(lead[,c("rsid","Metabolite")])
lead.tab = data.frame(table(lead$Metabolite))
colnames(lead.tab) = c("Metabolite","lead")

df = merge(df, lead.tab, by="Metabolite", all=T)

df$lead = ifelse(df$hits == 0,0,df$lead)


#



write.table(df, file="manuscript/figures/supplemental_table_2.txt", col.names=T, row.names=F, quote=F, sep="\t")
