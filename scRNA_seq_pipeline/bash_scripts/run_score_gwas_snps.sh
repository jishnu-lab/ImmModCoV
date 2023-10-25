#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J Score_gwas_snps
#SBATCH --output= #output file name
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user= #put your email id

echo 'started running'

##args
SNP_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/SNP_dir #directory containing all snp datasets
SNP_dataset=SNP_dataset #Directory containing SNP files
SNP_dataset_name=SNP_dataset #Name of SNP dataset
module purge
module load gcc/10.2.0 r/4.2.0

##Script
Rscript --no-save --no-restore --verbose /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/Score_gwas_snps.R $SNP_dir $SNP_dataset $SNP_dataset_name > /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/test.Rout 2>&1
