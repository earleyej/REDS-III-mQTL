#!/share/apps/R/bin/R

# Eric J. Earley
# RTI International
# April 2023


# this code will take processed metabolite data, pre-calculated eigenvectors, and other clinical data and create a phenotype file that will be used in ProbABEL

ANALYTE = "ATP"
OUTFILE = paste0("phenotype_",ANALYTE,"_N12433_20230502.txt")
ANALYTE_FILE = "phenotype_prep/redsiii_norm_impute_int_20230501.txt"


# columns should be:
# BSI.Subj.ID invlnkynurenine age male BCW BSRI ITxM BMI EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10
# 	note: BCW, BSRI, ITxM are blood donation sites; Probabel doesn't accept characters, so these will be encoded as 1|0 and each site will get a covariate.
setwd("/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/")

# QC'd metabolite data
analytes=read.table(ANALYTE_FILE, sep="\t", header=T, stringsAsFactors=F)




# demographic info from the original clinical data
demographic = read.table("phenotype_prep/rbcomics_metadata_for_metabolomics.csv",sep=",",header=T,stringsAsFactors=F)
# re-code Gender
demographic$male = ifelse(demographic$Gender == "M",0,1)
# expand and re-code HUB
demographic$BCW = ifelse(demographic$HUB == "BCW",1,0)
demographic$BSRI = ifelse(demographic$HUB == "BSRI",1,0)
demographic$ITxM = ifelse(demographic$HUB == "ITxM",1,0)
# remove Subject.ID and Gender
demographic = subset(demographic, select= -c(Gender, Subject.ID, HUB))


# Amy Moore already calculated PCs for this cohort. Just pull them from a phenotype file she created for phase1
# this file is missing numDonWithinTwoYear and there are other issuesr
#amy= read.table("phenotype_prep/kynurenine_13k_day42_narm.txt",sep=" ", header=T, stringsAsFactors=F)

# I also calculated these. Which one to use? Mine
eric=read.table("phenotype_prep/PCA/combined.allChr.pca.eigenvec",sep=" ",header=F,stringsAsFactors=F)
colnames(eric) = c("BSI.Subj.ID","FID",paste0("PC",1:20))
pcs = eric[,c("BSI.Subj.ID",paste0("PC",1:10))]
demographic = merge(demographic, pcs, by="BSI.Subj.ID")


# merge by subject with analyte
ANALYTE = "ATP"
idx = grep(ANALYTE,colnames(analytes))


pheno=merge(analytes[,c(1,idx)], 
  demographic,
  by.x="sample_id", 
  by.y="BSI.Subj.ID")


write.table(pheno, file=OUTFILE, sep=" ", col.names=T, row.names=F, quote=F)

