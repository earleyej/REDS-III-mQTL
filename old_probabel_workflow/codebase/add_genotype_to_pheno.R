rsid = "rs144986629"
chrom = "chr17"
metabolite = "DHexosephosphate"
phenoFile = paste0("/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/phenotype_",metabolite,".txt")
workingDir = paste0("/share/storage/REDS-III/RBCOmics/metabolomics/processing/",metabolite,"/")
imputeDir = paste0("/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01/",chrom)
#famFileDir = "/share/storage/REDS-III/RBCOmics/metabolomics/phenotypes/populations"
outFile = paste0(workingDir,"phenotype_",metabolite,"_",rsid,".txt")

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

write.table(pheno, outFile, sep="\t", col.names=T, row.names=F, quote=F)

