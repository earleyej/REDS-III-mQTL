#!/bin/sh
cd /share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/s1/$1 

working_dir=/share/storage/REDS-III/almoore/metabolomics/output/GWAS/continuous/s1/$1
imputation_root=/share/storage/REDS-III/RBCOmics/data/imputation/v1/imputations_all_maf0.01
phenotype_root=/share/storage/REDS-III/almoore/metabolomics/input/phenotypes/continuous/s1

method=palinear 

for (( chr=1; chr<24; chr++ )); do 
	out_file=v09272021.rbc.all.1000G_p3.chr${chr}.${1}~age+gender+hub+donate+evs+SNP.stats 
	phenotype_file=$1.txt 
	geno_prefix=rbc.ALL.1000G_p3.chr
	
		/share/storage/REDS-III/common/software/nextflow/nextflow-0.25.1-all \
		/share/storage/REDS-III/common/software/pipelines/_pipeline.association.out_stats_files.v0.1.nf \
		--final_chunks $imputation_root/chunks/final_chunks.chr$chr \
		--input_pheno $phenotype_root/$phenotype_file \
		--imputation_dir $imputation_root/chr$chr \
		--example_mldose $imputation_root/chr$chr/rbc.ALL.1000G_p3.chr$chr.0.mach.mldose.gz \
		--geno_prefix $geno_prefix \
		--working_dirs $working_dir \
		--out $working_dir/$out_file \
		--method $method
	done 
	
