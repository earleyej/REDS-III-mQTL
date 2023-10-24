# Select multiple associated SNPs through a stepwise selection procedure
#gcta64  --bfile test  --chr 1 --maf 0.01 --cojo-file test.ma --cojo-slct --out test_chr1


#test.ma format:
#SNP A1 A2 freq b se p N 
#rs1001 A G 0.8493 0.0024 0.0055 0.6653 129850 
#rs1002 C G 0.0306 0.0034 0.0115 0.7659 129799 
#rs1003 A C 0.5128 0.0045 0.0038 0.2319 129830
#...

# A1 is effect allele; A2 is the other allele
# must provide the complete GWAS summary stats
# 'test' bfile is for the whole chromosome; it is used to calculate LD

[ -z "$ancestry" ] && export ancestry=ALL


#### first, make the input.ma ####
zcat $processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz | awk '{OFS="\t"} {print $4,$6,$5,$7,$10,$11,$12,$9}' | sed 's/\tX\t/\t23\t/' > $processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.cojo.ma


# make a list of chromosomes that contain at least 1 significant SNP
export chroms=`cat ${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.gz.sig | cut -f2 | sort | uniq | sed 's/X/23/'`


#### Run cojo-slct by chrom ###
for chr in $chroms; do
  /home/eearley/apps/qsub_job.sh \
    --job_name cojo_${chr}_$METABOLITE \
    --script_prefix ${processingDir}/cojo_${chr} \
    --mem 8 \
    --priority 0 \
    --cpu 1 \
    --program /home/eearley/apps/gcta-1.94.1-linux-kernel-3-x86_64/gcta-1.94.1 \
      --bfile /share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01/chr${chr}/chr${chr} \
      --chr ${chr} \
      --maf 0.01 \
      --cojo-file $processingDir/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.cojo.ma \
      --cojo-slct \
      --out $processingDir/cojo_chr${chr}
done
#


