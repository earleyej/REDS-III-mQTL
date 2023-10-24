library(data.table) 
library(impute)
library(pcaMethods)
library(imputeLCMD) #this loads package 'norm' which is superceded by 'norm2' -- this may be a problem


#written by Eric J Earley
# adapted from original code by Amy Moore
# RTI International 2023

# CHANGE LOG
# April 27, 2023
#	updating the Additive Solution breakdown. Only HUB BSRI used AS3; the other three HUBs used AS1




#### LOAD RAW DATA ####
setwd("/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes")
raw = read.csv("/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/phenotype_prep/221109_REDSIII_norm_values_renamed.csv")
dim(raw)
# 13091 375


# merge the metabolite data with some other metadata for these participants. 
# 'phenofull' contains just BMI, sex, age at enrollment, and BSI Subj ID, as well as sampleID
phenofull=read.csv("phenotype_prep/rbcomics_metadata_for_metabolomics.csv")
raw=merge(raw, phenofull, by.x="sample_id", by.y="BSI.Subj.ID", all.x=F, all.y=F)

#set 0 to NA
raw[raw == 0] <- NA





#### CHECK DRUGS/EXOGENOUS METABOLITES"
exogenous_metabolites <- c("X4Hydroxywarfarin", "X6Methylprednisone", "acetaminophen", "alprazolam", "amlodipine", "amoxicillin", "AMPHETAMINE", "atorvastatin", "bupropion", "cannabidiol", "carboxyibuprofen", "carvedilol", "clonazepam", "cocaine", "cotinine", "CotinineNoxide", "desmethylranitidine", "diazepam", "doxycycline", "enalaprilat", "Escitalopram", "famotidine", "Ferricgluconate", "fexofenadine", "fluconazole", "fluoxetine", "gabapentin", "glipizide", "Hydroxycotinine", "HYDROCODONE", "ketamine", "lidocaine_208", "lidocaine_264", "lisinopril", "lorazepam", "losartan", "meloxicam", "methocarbamol", "metoprolol", "metronidazole", "NAPROXEN", "Nicotinamide", "nicotinamideriboside", "Nicotinateribonucleotide", "nicotine", "norcotinine", "omeprazole", "Oxyresveratrol", "pioglitazone", "PSILOCIN", "ranitidine", "RanitidineNoxide", "Retinol", "risperidone", "sertraline", "simvastatin", "Trazadone", "valsartan", "verapamil", "zolpidem")
drugs <- subset(raw, select=exogenous_metabolites)

# Calculate some basic statistics
sumstats <- function(x) {
	Unobserved <- apply(x, 2, function(i) sum(is.na(i)))
	Mean <- apply(x, 2, function(i) mean(i, na.rm=T))
	Median <- apply(x, 2, function(i) median(i, na.rm=T))
	SD <- apply(x, 2, function(i) sd(i, na.rm=T))
	CV <- apply(x, 2, function(i) sd(i, na.rm=T)/mean(i, na.rm=T))
	result <- data.frame(Unobserved ,Mean, Median, SD, CV)
	return(result)
}

drugstats <- (sumstats(drugs))
# EJE: I'm keeping this file name in case Amy uses it later in another script
#write.table(drugstats, file="221109_REDSIII_norm_values_renamed_drugstats_20230125.txt", sep="\t", quote=F, row.names=T)


# remove exogenous drugs. These won't be associated with genetics
# and things like cocaine, ketamine, etc. are likely mislabeled anyway
baddrugs <- rownames(drugstats[which(drugstats$Median > 0),]) 
raw_nobaddrugs <- raw[,!(names(raw) %in% baddrugs)]
#drugs_nobaddrugs <- drugs[,!(names(drugs) %in% baddrugs)]



#### PROCESS METABOLITES ####
# These data were generated using "Additive solutions" which bias results.
# Stratify this QC process by AS1 and AS3 (there is no AS2)
# ARC  BCW BSRI ITxM
#3174 3272 3564 3019
raw_nobaddrugs$additive <- ifelse(raw_nobaddrugs$HUB %in% c("BSRI"), "AS3", "AS1")

# remove missingness HERE before stratifying
# put aside the non-metabolite columns, filter, then add them back
admin_vars = c("sample_id","subject_id","Gender","age","numDonWithinTwoYear","HUB","Subject.ID","additive")
data_vars = raw_nobaddrugs[,!colnames(raw_nobaddrugs) %in% admin_vars]
data_vars = data_vars[, which(colMeans(!is.na(data_vars)) >= 0.70)]
raw_nobaddrugs = cbind(raw_nobaddrugs[,admin_vars],data_vars)

data_1 <- subset(raw_nobaddrugs, additive=="AS1") #
data_3 <- subset(raw_nobaddrugs, additive=="AS3") #

# keep the columns in common between the two AS groups
#data_1_colnames <- colnames(data_1)
#data_3_colnames <- colnames(data_3)
#joint_colnames <- intersect(data_1_colnames, data_3_colnames) #
#unique_colnames <- setdiff(data_1_colnames, data_3_colnames)
#data_1 <- data_1[,joint_colnames] 
#data_3 <- data_3[,joint_colnames]  


# LOG TRANSFORM
data1rel = apply(data_1[,!colnames(data_1) %in% admin_vars],2,log)
data3rel = apply(data_3[,!colnames(data_3) %in% admin_vars],2,log)

#### IMPUTE MISSING DATA ####
data1rel_imputed <- impute.QRILC(data1rel, tune.sigma=1)
data3rel_imputed <- impute.QRILC(data3rel, tune.sigma=1)

imputed.data1rel_imputed <- data1rel_imputed[[1]] 
imputed.data3rel_imputed <- data3rel_imputed[[1]] 

# only keep the subject IDs
imputed.data1rel_imputed <- cbind(data_1[,c(1:2)], imputed.data1rel_imputed)
imputed.data3rel_imputed <- cbind(data_3[,c(1:2)], imputed.data3rel_imputed)

imputed.final <- rbind(imputed.data1rel_imputed, imputed.data3rel_imputed)

# Inverse normal transformation
# There are still long tails for some of these metabolites. This INT will adjust so the distributions look more normal
# this requires the GenABEL package which is not easily installed on HPC.
# perform the INT on a local computer after stripping away the sample IDs
#write.table(imputed.final[,!colnames(imputed.final) %in% admin_vars], file="phenotype_prep/redsiii_toBeINT.txt",sep="\t",col.names=T,row.names=F, quote=F)

x.int = read.table("phenotype_prep/redsIII_metabolites_INT_2023_05_05.txt",sep="\t",header=T,stringsAsFactors=F)

out.final = cbind(imputed.final[,c(1,2)],x.int)

write.table(out.final, file="phenotype_prep/redsiii_norm_impute_int_20230505.txt", sep="\t", quote=F, row.names=F, col.names=T)

