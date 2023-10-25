#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 1-00:00 # Runtime in D-HH:MM
#SBATCH -J create_single_cell_count_matrix
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=prg65@pitt.edu

##args
RNA_seq_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/RNA_dataset/input_dataset #directory with rna seq object
RNA_seq_file=example_processed_dataset.rds #rna seq object
output_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/RNA_dataset/output_dataset #directory to store cell type specific datasets
phenotype=COVID #phenotype label used in status column

module purge
module load gcc/10.2.0 r/4.2.0

##Script
Rscript --no-save --no-restore --verbose /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/create_cell_type_specific_RNA_dataset.R $RNA_seq_dir $RNA_seq_file $output_dir $phenotype > /ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/test.Rout 2>&1
