#Eric Earley
# RTI International, 2023
# this script will take as arguments:
#		rsID, chrom number, phenoFile
# and return boxplots of metabolite vs. genotype
# as well as per-ancestry plots


# input = rsID
# input = chrom number
# input = phenotype file

rsid = "rs144986629"
chrom = "chr17"
phenoFile = "/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/phenotype_DHexosephosphate.txt"
workingDir = "/share/storage/REDS-III/RBCOmics/metabolomics/processing/DHexosephosphate/"
imputeDir = paste0("/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01/",chrom)
famFileDir = "/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/populations"


# step 1
# search for the correct chrom chunk from the imputation directory using the mlinfo files
cmd = paste0("grep -n ",rsid," ",imputeDir,"/*mlinfo")
stdout = system(cmd, intern=TRUE)
parse.me=strsplit(stdout,split="\t")[[1]][1]
mlinfoFile=strsplit(parse.me,split=":")[[1]][1]
mlinfoLine=as.numeric(strsplit(parse.me,split=":")[[1]][2])

# step 2 
# pull out the genotype from the mldose file and add it to the phenotype
mldoseFile=gsub("mlinfo","mldose.gz",mlinfoFile)
mldose = read.table(mldoseFile,sep=" ",header=F,stringsAsFactors=F)
# columns are loci, rows are samples
column = mlinfoLine + 2
add.me = data.frame("SampleID" = mldose[,1],
		   "genotype" = mldose[,column])



# step 3 
# gather the per ancestry .fam file and add population to the phenotype
pheno =read.table(phenoFile, sep=" ",header=T,stringsAsFactors=F)
pheno = merge(pheno, add.me, by.x="sample_id", by.y="SampleID")

# gather sample IDs by population
afr = read.table(paste0(famFileDir,"/rbc.AFRAMRCN.1000G_p3.chr22.0.fam"), sep="\t",header=F)
asn = read.table(paste0(famFileDir,"/rbc.ASIAN.1000G_p3.chr22.0.fam"), sep="\t",header=F)
eur = read.table(paste0(famFileDir,"/rbc.CAUCASIAN.1000G_p3.chr22.0.fam"), sep="\t",header=F)
his = read.table(paste0(famFileDir,"/rbc.HISPANIC.1000G_p3.chr22.0.fam"), sep="\t",header=F)
oth = read.table(paste0(famFileDir,"/rbc.OTHER.1000G_p3.chr22.0.fam"), sep="\t", header=F)

# redundant with code below
#pheno$pop = NA
#pheno$pop = ifelse(pheno$sample_id %in% fam_afr$V1, "AFR", pheno$pop)
#pheno$pop = ifelse(pheno$sample_id %in% fam_asn$V1, "ASN", pheno$pop)
#pheno$pop = ifelse(pheno$sample_id %in% fam_eur$V1, "EUR", pheno$pop)
#pheno$pop = ifelse(pheno$sample_id %in% fam_his$V1, "HIS", pheno$pop)
#pheno$pop = ifelse(pheno$sample_id %in% fam_oth$V1, "OTH", pheno$pop)


# now we have: metabolite, genotype, population

# step 4
# plot

plot.me = function(pheno, title, pop) {
  tmp = pheno[pheno$sample_id %in% pop[,1],]
  title=paste(title,", N=",dim(tmp)[1])
  metabolite=colnames(pheno)[2]
  # reverse order of the genotype groups
  tmp$genotype = ifelse(tmp$genotype == 0, "Alt/Alt",
                        ifelse(tmp$genotype == 1, "Het","Ref/Ref"))
  tmp$genotype = factor(tmp$genotype, levels = c("Ref/Ref","Het","Alt/Alt"))
  p2 = ggplot(data=tmp, aes(x=genotype,y=.data[[metabolite]])) +
    geom_boxplot() +
    theme_classic() +
    ylab("") +
    xlab("") +
    ggtitle(title)
  return(p2)  
}
p.all=plot.me(pheno, title="ALL",pop=pheno)
p.afr=plot.me(pheno,title="AFR",pop=afr)
p.asn=plot.me(pheno,title="ASN",pop=asn)
p.eur=plot.me(pheno,title="EUR",pop=eur)
p.his=plot.me(pheno,title="HIS",pop=his)
p.oth=plot.me(pheno,title="OTH",pop=oth)
plot_list=list(p.all,p.afr,p.asn,p.eur,p.his,p.oth)

png(file=paste0(workingDir,"/DHexosephosphate_vs_genotype.png"),res=600,units="cm",height=20,width=25)
do.call("grid.arrange", c(plot_list, ncol = 3, bottom="", left="Hexosephosphate") )  
dev.off()





