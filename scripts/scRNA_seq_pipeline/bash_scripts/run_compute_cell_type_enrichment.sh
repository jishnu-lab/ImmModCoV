#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH -J Compute_per_gene_prop_cells_odds_ratio
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-3 # job array index #needed only if submitting multiple SNP datasets together
#SBATCH --output=test.out

RNA_seq_folder=`ls /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/RNA_dataset/output_dataset| head -n $SLURM_ARRAY_TASK_ID | tail -n 1`
echo $RNA_seq_folder
##args
SNP_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/SNP_dir
SNP_dataset=SNP_dataset
SNP_file=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/bash_input/SNP_file_names.txt
phenotype=COVID
RNA_seq_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/RNA_dataset/output_dataset
output_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/plots_output/enrichment_plots
snp_type=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/bash_input/SNP_file_types.txt
module purge
module load gcc/10.2.0 r/4.2.0

##Script
Rscript --no-save --no-restore --verbose /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/Compute_cell_type_enrichment.R $SNP_dir $SNP_dataset $SNP_file $phenotype $RNA_seq_dir $RNA_seq_folder $output_dir $snp_type > /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/test.Rout 2>&1
