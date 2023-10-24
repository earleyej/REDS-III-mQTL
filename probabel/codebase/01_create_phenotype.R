#!/share/apps/R/bin/R

# Eric J. Earley
# RTI International
# April 2023


# this code will take processed metabolite data, pre-calculated eigenvectors, and other clinical data and create a phenotype file that will be used in ProbABEL
args <- commandArgs(TRUE)
#if (hasArg(args)) {
#	stop("Usage:\nRscript 01_create_phenotype.R --metaboliteFile <.txt> --pc <.eigenvec> --metabolite <'ATP'> --demographic <.txt> --out <.txt>\n\n")
#}
loop = TRUE
while (loop) {

        if (args[1] == "--metaboliteFile") {
                ANALYTE_FILE = args[2]
                hasAnalyteFile = TRUE
		cat(paste0("\n\nMetabolite file: ",ANALYTE_FILE, "\n"))
        }
	if (args[1] == "--out") {
		OUTFILE = args[2]
		hasOutFile = TRUE
	}
	if (args[1] == "--pc") {
		PCFILE = args[2]
		hasPCFile = TRUE
		cat(paste0("PC file: ",PCFILE,"\n"))
	}
	if (args[1] == "--metabolite") {
		ANALYTE = args[2]
		hasAnalyte = TRUE
		cat(paste0("Metabolite: ",ANALYTE,"\n"))
	}
	if (args[1] == "--demographic") {
		demographicFile = args[2]
		hasDemographicFile = TRUE
		cat(paste0("Demographic file: ",demographicFile,"\n"))
	}
        if (length(args) > 1) {
                args = args[2:length(args)]
        } else {
                loop=FALSE
        }
}



# columns should be:
# BSI.Subj.ID invlnkynurenine age male BCW BSRI ITxM BMI EV1 EV2 EV3 EV4 EV5 EV6 EV7 EV8 EV9 EV10
# 	note: BCW, BSRI, ITxM are blood donation sites; Probabel doesn't accept characters, so these will be encoded as 1|0 and each site will get a covariate.
#setwd("/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/")

cat("\nCombining PCs, demographics, and metabolite data...\n\n")
# QC'd metabolite data
analytes=read.table(ANALYTE_FILE, sep="\t", header=T, stringsAsFactors=F)





# demographic info from the original clinical data
demographic = read.table(demographicFile,sep=",",header=T,stringsAsFactors=F)
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
eric=read.table(PCFILE,sep=" ",header=F,stringsAsFactors=F)
colnames(eric) = c("BSI.Subj.ID","FID",paste0("PC",1:20))
pcs = eric[,c("BSI.Subj.ID",paste0("PC",1:10))]
demographic = merge(demographic, pcs, by="BSI.Subj.ID")


# merge by subject with analyte
#ANALYTE = "ATP"
idx = grep(ANALYTE,colnames(analytes))


pheno=merge(analytes[,c(1,idx)], 
  demographic,
  by.x="sample_id", 
  by.y="BSI.Subj.ID")


write.table(pheno, file=OUTFILE, sep=" ", col.names=T, row.names=F, quote=F)
cat(paste0("Created phenotype file: ", OUTFILE,"\n\n"))
