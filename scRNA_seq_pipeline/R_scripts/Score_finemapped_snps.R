##Score GWAS SNPs##
library(stringr)
library(tidyr)
library(tidyverse)
source('/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/Functions_for_snp_scoring.R')
args <- commandArgs(trailingOnly = TRUE) 
print(args)
SNP_dir <- args[1]
SNP_dataset <- args[2] 
SNP_dataset_name <- args[3]
print(SNP_dataset_name)
output_file_name <- paste0(SNP_dataset_name, ".", "finemapped.snps.scored.txt")
snp_file_abc <-  read.table(paste0(SNP_dir,'/',SNP_dataset,'/', SNP_dataset_name, ".", "finemapped.anno.abcroadmap.bed"))
snp_file_2kb <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/', SNP_dataset_name, ".", "finemapped.anno.2kb.bed"))
snp_file_gene_body <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/', SNP_dataset_name, ".", "finemapped.anno.gene.body.bed"))
protein_coding_snps <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/', SNP_dataset_name, ".", "protein.coding.snps.txt"))
print(SNP_dataset_name)
##Process the gws snps files
snp_file_abc <- snp_file_abc[, c('V1', 'V2', 'V3', 'V7', 'V8')]
colnames(snp_file_abc) <- c('chr', 'start', 'stop', 'GENE', 'ANNO')
snp_file_2kb <- snp_file_2kb[, c('V1', 'V2', 'V3', 'V7', 'V8')]
colnames(snp_file_2kb) <- c('chr', 'start', 'stop', 'GENE', 'ANNO')
snp_file_gene_body <- snp_file_gene_body[, c('V1', 'V2', 'V3', 'V7', 'V8')]
colnames(snp_file_gene_body) <- c('chr', 'start', 'stop', 'GENE', 'ANNO')
#combine the dfs 
snp_file <- rbind(snp_file_2kb, snp_file_abc)
snp_file <- rbind(snp_file, snp_file_gene_body)
snp_file$SNP <- paste0(snp_file$chr, ':', snp_file$start)
snp_file <- snp_file[, c('SNP', 'GENE', 'ANNO')]
snp_file_wo_dup <- snp_file[!duplicated(snp_file),]
#remove protein_coding snps 
snp_file_non_coding <- snp_file_wo_dup[!(snp_file_wo_dup$SNP %in% protein_coding_snps$V1),]
protein_coding_snps_gws <- protein_coding_snps[which(protein_coding_snps$V1 %in% snp_file_wo_dup$SNP),]
#divide a score of 1 across all the annotations and genes that a snp maps to#
y <- split(snp_file_non_coding, snp_file_non_coding$SNP)
gws_snps <- data.frame("SNP" = character(), "GENE" = character(), "ANNO" = character(), "SCORE" = numeric())
for (snp_df in 1:length(y)) {
  gws_snps <- rbind(gws_snps, score_gwas_snps(y[[snp_df]]))
}
gws_snps$CHR <- str_remove(gws_snps$SNP, ":.*")
gws_snps$SNP <- str_remove(gws_snps$SNP, ".*:")
gws_snps <- gws_snps[, c("CHR", "SNP", "GENE", "SCORE")]
##Add scores for protein coding GWAS SNPs##
intreface_snps <- read.table(paste0(SNP_dir,'/', SNP_dataset,'/',SNP_dataset_name, ".", "interface.snps.scored.txt"))
not_interface_snps <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/',SNP_dataset_name, ".", "not.interface.snps.scored.txt")) 
protein_coding_snps_scored <- rbind(not_interface_snps, intreface_snps)
colnames(protein_coding_snps_scored) <- c("CHR", "SNP", "GENE", "SCORE") 
protein_coding_snps_scored$SNP_id <- paste0(protein_coding_snps_scored$CHR, ":", protein_coding_snps_scored$SNP)
protein_coding_snps_gws_scored <- protein_coding_snps_scored[protein_coding_snps_scored$SNP_id %in% protein_coding_snps_gws,]
protein_coding_snps_gws_scored <- protein_coding_snps_gws_scored[, c("CHR", "SNP", "GENE", "SCORE")]
gws_snps <- rbind(gws_snps, protein_coding_snps_gws_scored)
#Save table 
write.table(gws_snps, file = paste0(SNP_dir, '/', SNP_dataset, "/", output_file_name), row.names = F, col.names = F, sep = '\t', quote = F)