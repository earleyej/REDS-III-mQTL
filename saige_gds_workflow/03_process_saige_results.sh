# config file sets these environment variables:
# processingDir
# METABOLITE
# scriptDir

# if $ancestry is empty, set as "ALL"
[ -z "$ancestry" ] && export ancestry=ALL 


export outFile=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats
echo -e "id\tchr\tpos\trs.id\talt\tref\tmaf\tmac\tnum\tbeta\tSE\tpval\ttype" > $outFile


# saigeGDS results are seprated into chunks (~1.3k). This code will concatenate into one outFile
for (( chr=1; chr<24; chr++ )); do
  echo $chr
  for chunk in $(grep -P "^$chr\s+" /share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_by_race/final_chunks | perl -lane 'print $F[1];'); do
    baseName=rbc.${ancestry}.1000G_p3.chr${chr}.${chunk}.assoc.txt #rbc.ALL.1000G_p3.chr5.28.assoc.txt

    # add a column for SNP type: "SNP" or "INDEL"
    # replace chromosome "X" with "23" for locuszoom
    # if pval column ($12) is == 0, replace it with a non-zero pvalue (1.0e-350)
    # convert AF to MAF
    tail -n +2 ${processingDir}/${baseName} | awk 'BEGIN {OFS="\t"} {if ( ($5=="A" || $5=="C" || $5=="G" || $5=="T") && ($6=="A" || $6=="C" || $6=="G" || $6=="T")) {print $0,"SNP"} else {print $0,"INDEL"} }' | sed 's/chr//' | sed 's/\tX/\t23\t/' | awk '{OFS="\t"} { $12 = ($12 == 0?"1e-350":$12) } 1' | awk '{OFS="\t"} { $7 = ($7 > 0.5 ? 1-$7 :$7) } 1' >> $outFile
  done
done

echo "compressing the combined stats file.."
gzip -f $outFile

#### now subset columns for easier transport ####
echo "creating gwas summary stats subset for FUMA"
portable_file=${processingDir}/${METABOLITE}_rbc.${ancestry}.1000G_p3.age_gender_hub_donate_10evs_SNP.stats.forFUMA.gz
zcat ${outFile}.gz | awk '{FS=" "} {print $4,$2,$3,$5,$6,$7,$10,$11,$12}' | tr " " "\t" | gzip > $portable_file


