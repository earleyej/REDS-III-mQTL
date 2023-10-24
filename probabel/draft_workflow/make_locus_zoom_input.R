#library(data.table)
library(dplyr)

# Eric J. Earley
# RTI International
# May 4, 2023


## WORK IN PROGRESS ###
sigFile="/share/storage/REDS-III/RBCOmics/metabolomics/processing/ATP/ATP_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz.sig"
statsFile="/share/storage/REDS-III/RBCOmics/metabolomics/processing/ATP/ATP_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz"

sig = read.table(sigFile, sep=" ", header=F, stringsAsFactors=F)
colnames(sig)[c(1,2,15,17)]=c("VARIANT_ID","POSITION","P","CHR")

# This only works if there is 1 lead SNP per chrom
# change this to be lowest P within a 500kb region
tmp = data.frame(sig %>% 
	group_by(CHR) %>% 
	slice_min(order_by = P))


# from the full stats file, pull out 200kb flanking each lead SNP
#cat(paste0("Reading stats file...", statsFile, "\n"))
#gwas = read.table(gzfile("ATP_rbc.all.1000G_p3.age_gender_hub_donate_evs_SNP.stats.gz"),sep=" ", header=T, stringsAsFactors=F)
# THIS TAKES FOREVER... Maybe just cut out a section at a time


i=1
for (i in 1:length(tmp$VARIANT_ID)) {
  cmd=paste0("zcat ",statsFile," | grep -C 50000 ",tmp$VARIANT_ID[i]," > tmp",i)
  print(paste0("Creating metal-formatted region file: ",tmp$VARIANT_ID[i]))
  system(cmd)
  region=read.table(paste0("tmp",i), sep=" ", header=F, stringsAsFactors=F)
  region.sub=region[ (region$V2 > (tmp$POSITION[i] - 100000) & region$V2 < (tmp$POSITION[i] + 100000)),]
  out=region.sub[,c("V1","V15")]
  colnames(out)=c("MarkerName","P-value")
  outFile=paste0("locuszoom_",tmp$CHR[i],"_",tmp$POSITION[i],"_",tmp$VARIANT_ID[i],".metal")
  write.table(out, file=outFile, sep="\t",col.names=T,row.names=F,quote=F)
  print(paste0("Created file: ",outFile))
}
