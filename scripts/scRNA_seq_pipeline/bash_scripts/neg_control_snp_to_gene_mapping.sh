#!/bin/bash
#
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 0-01:00 # Runtime in D-HH:MM
#SBATCH -J gws_snp_bedtools
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=
#SBATCH --output=gws_snp_bedtools.out

##Script to map GWAS SNPs to genes

SNP_dir=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/example/SNP_dir/SNP_dataset #directory with snp bed file
SNP_file=SNP_dataset.neg.snps.txt #name of snp file
Output_file_prefix=SNP_dataset #name you want output file to start with
ref_abc_file=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/ref_files_for_snp_mapping/RoadmapUABC_with_anno_BLD.txt
ref_gene_body=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/ref_files_for_snp_mapping/Gene_length.bed
ref_2kb_file=/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/ref_files_for_snp_mapping/Gene_2kb_upstream.bed

#load modules
module purge
module load gcc/8.2.0 bedtools/2.30.0

cd $SNP_dir

##Overlap with abc file
bedtools intersect -a $SNP_file -b $ref_abc_file -wa -wb > $Output_file_prefix.neg.anno.abcroadmap.bed
##Overlap with 2kb file
bedtools intersect -a $SNP_file -b $ref_2kb_file -wa -wb > $Output_file_prefix.neg.anno.2kb.bed
##Overlap with gene body file
bedtools intersect -a $SNP_file -b $ref_gene_body -wa -wb > $Output_file_prefix.neg.anno.gene.body.bed
