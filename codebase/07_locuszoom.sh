
scriptDir=/share/storage/REDS-III/RBCOmics/metabolomics/workflow/codebase


###########################################
#### step 1 - make a list of lead SNPs ####
###########################################
# start with the .sig
# pull out the rsID with the lowest pval in a 200kb region
Rscript $scriptDir/_make_list_of_lead_SNPs.R



###################################################################
#### step 2 - calculate LD for each lead SNP in a 200kb window ####
###################################################################
# plink
# To obtain all LD values for a set of SNPs versus one specific SNP, use the --ld-snp command in conjunction with --r2. For example, to get a list of all values for every SNP within 1Mb of rs12345, use the command
#plink --file mydata 
#          --r2 
#          --ld-snp rs12345 #the tag snp
#          --ld-window-kb 1000 #only calculate LD for SNPs in this window
#          --ld-window 99999 #only consider SNPs that are 99,999 SNPs apart from tag (set this arbitrarily high to get all SNPs)
#          --ld-window-r2 0 #minimum LD to report; 0 means that all results will be reported
rsid_list=`cat lead_snps`

for rsid in $rsid_list; do
  echo $rsid
#  plink --bfile chr1 \
#    --r2 \
#    --ld-snp rs12345 \
#    --ld-window-kb 200 \
#    --ld-window 99999 \
#    --ld-window-r2 0 \
#    --out rs12345
done

# step 3 - run locuszoom for each lead SNP

