library(stringr)
library(tidyverse)
source('/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/Functions_for_snp_scoring.R')
args <- commandArgs(trailingOnly = TRUE) 
print(args)
SNP_dir <- args[1]
SNP_dataset <- args[2]
SNP_dataset_name <- args[3]
output_file_name <- paste0(SNP_dataset_name, ".", "neg.snps.scored.txt")
snp_file_abc <-  read.table(paste0(SNP_dir,'/',SNP_dataset,'/',SNP_dataset_name, ".", "neg.anno.abcroadmap.bed"))
snp_file_2kb <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/',SNP_dataset_name, ".", "neg.anno.2kb.bed"))
snp_file_gene_body <- read.table(paste0(SNP_dir,'/',SNP_dataset,'/',SNP_dataset_name, ".", "neg.anno.gene.body.bed"))
##Process the neg snps files
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
#divide a score of 1 across all the annotations and genes that a snp maps to#
y <- split(snp_file_wo_dup, snp_file_wo_dup$SNP)
gws_snps <- data.frame("SNP" = character(), "GENE" = character(), "ANNO" = character(), "SCORE" = numeric())
for (snp_df in 1:length(y)) {
  gws_snps <- rbind(gws_snps, score_gwas_snps(y[[snp_df]]))
}
gws_snps$CHR <- str_remove(gws_snps$SNP, ":.*")
gws_snps$SNP <- str_remove(gws_snps$SNP, ".*:")
gws_snps <- gws_snps[, c("CHR", "SNP", "GENE", "SCORE")]
#Save table 
write.table(gws_snps, file = paste0(SNP_dir, '/', SNP_dataset, "/", output_file_name), row.names = F, col.names = F, sep = '\t', quote = F)