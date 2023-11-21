###Compute snp type score for each cell type##
library(R.utils)
library(tidyr)
library(Matrix.utils)
library(stringr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(sparseMatrixStats)
source('/ix/djishnu/Priyamvada/RNA_seq_GWAS_integrated_analysis/R_scripts/Functions_to_source_prop_per_gene.R')
args <- commandArgs(trailingOnly = TRUE)
print(args)
SNP_dir <- args[1]
SNP_dataset <- args[2]
SNP_file <-  args[3] #file with all snp names
phenotype <- args[4]
RNA_seq_dir <- args[5]
RNA_seq_folder <- args[6]
output_dir <- args[7]
snp_type <- args[8] #file with snp types
SNP_files <- read.table(SNP_file)
SNP_files <- SNP_files$V1
SNP_types <- read.table(snp_type)
SNP_types <- SNP_types$V1
RNA_seq_files <- list.files(paste0(RNA_seq_dir, '/', RNA_seq_folder), pattern = '.rds')
cell_types <- str_remove(RNA_seq_files, '_count_matrix.rds')
for (cell in 1:length(cell_types)){
  cell_type_RNA_seq <- readRDS(paste0(RNA_seq_dir, '/', RNA_seq_folder, '/', RNA_seq_files[cell]))
  plot_df <- data.frame("snp_names" = character(), "gene_names" = character(),  "prop_cells" = numeric(), 'snp_type' = character(),  'status' = character(), 'n_genes' = numeric())
  SNP_df_for_plotting <- data.frame("CHR" = character(), "POS" = integer(), "GENE" = character(), "SCORE" = numeric(), "SNP_type" = character(), "n_genes" = numeric())
  SNP_df <- data.frame("CHR" = character(), "POS" = integer(), "GENE" = character(), "SCORE" = numeric())
  for (snp in 1:length(SNP_files)) {
    SNP_df_type <- read.table(paste0(SNP_dir, '/', SNP_dataset, '/', SNP_files[snp]))
    SNP_df_type_plotting <- SNP_df_type
    n_genes <- length(unique(SNP_df_type$V3))
    SNP_type <- SNP_types[snp]
    SNP_df_type_plotting$SNP_type <- SNP_type
    SNP_df_type_plotting$n_genes <- n_genes
    SNP_df <- rbind(SNP_df, SNP_df_type)
    SNP_df_for_plotting <- rbind(SNP_df_for_plotting, SNP_df_type_plotting)
  }
  create_snp_score_matrix(SNP_df = SNP_df)
  plot_df <- Gene_wise_enrichment_healthy_covid_sep(SNP_matrix = SNP_mat, RNA_seq = cell_type_RNA_seq, phenotype = phenotype, plot_df = plot_df, threshold = 0.95, snp_df = SNP_df_for_plotting, cell_prop_threshold = 0.95)
  plot_df_healthy <- plot_df[plot_df$status == 'healthy',]
  plot_df_phenotype <- plot_df[plot_df$status == phenotype,]
  plot_df_list <- list(plot_df_healthy, plot_df_phenotype)
  df_for_plotting <- data.frame('SNP_TYPE' = character(), "FRACTION" = numeric(), 'STATUS' = character())
  for (i in 1:length(plot_df_list)) {
    fraction_df <- plot_df_list[[i]]
    n_genes_snp_type <- fraction_df[, c("snp_type", "n_genes")] %>% unique()
    unique_genes <- fraction_df %>% group_by(snp_type) %>% summarise(n=n())
    unique_genes <- left_join(unique_genes, n_genes_snp_type, by = 'snp_type')
    unique_genes$fraction <- unique_genes$n/unique_genes$n_genes
    unique_genes$SE <- sqrt((unique_genes$fraction*(1 - unique_genes$fraction))/unique_genes$n_genes)
    unique_genes$UPPER_CI <- unique_genes$fraction + unique_genes$SE*(qnorm(0.975))
    unique_genes$LOWER_CI <- unique_genes$fraction - unique_genes$SE*(qnorm(0.975))
    unique_genes$status <- unique(fraction_df$status)
    p_val_df <- data.frame("group1" = character(), 'group2' = character(), "p_val" = numeric())
    snps_to_test <- unique(unique_genes$snp_type)
    snps_to_test <- snps_to_test[which(snps_to_test != 'INTERFACE')]
    for (snp in 1:length(snps_to_test)) {
      x_interface <- unique_genes[which(unique_genes$snp_type == 'INTERFACE'),]$n
      x_snp <- unique_genes[which(unique_genes$snp_type == snps_to_test[snp]),]$n
      successes <- c(x_interface, x_snp)
      n_interface <- unique_genes[which(unique_genes$snp_type == 'INTERFACE'),]$n_genes
      n_snp <- unique_genes[which(unique_genes$snp_type == snps_to_test[snp]),]$n_genes
      trials <- c(n_interface, n_snp)
      p_val_snp <- prop.test(successes, trials, correct = F, alternative = 'greater')$p.value
      if (is.na(p_val_snp)) {
        p_val_snp <- 1
      }
      p_val_df <- rbind(p_val_df, data.frame("group1" = "INTERFACE", "group2" = snps_to_test[snp], "p_val" = p_val_snp))
    }
    p_val_df$status <- unique(unique_genes$status)
    df_for_plotting <- rbind(df_for_plotting, unique_genes[, c('snp_type', 'fraction', 'SE', 'UPPER_CI', 'LOWER_CI', 'status')])
    p_val_df <- transform(p_val_df, stars = stars_pval_cell_type_bxplot(p_val_df$p_val))
  }
  df_for_plotting <- df_for_plotting[which(df_for_plotting$status == phenotype),]
  p_val_df <- p_val_df[which(p_val_df$status == phenotype),]
  p_val_df <- p_val_df[, c('group1', 'group2', 'p_val', 'stars')]
  if(dir.exists(paste0(output_dir, '/', RNA_seq_folder)) == FALSE) {
    dir.create(paste0(output_dir, '/', RNA_seq_folder))
    print("RNA directory does not exist")} else {
      print("RNA directory exists")
    }
  if(dir.exists(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset)) == FALSE) {
    dir.create(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset))} else {
      print("SNP directory exists")
    }
  pdf(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset, '/', cell_types[cell], '_', 'compare.pdf'), height = 4, width = 6, paper = "USr")
  print(ggbarplot(data=df_for_plotting, x="snp_type", y="fraction", fill = "snp_type") +
          geom_errorbar(aes(ymin=fraction, ymax=fraction+SE,), width=.2,
                        position=position_dodge(.9)) +
          stat_pvalue_manual(p_val_df, y.position = (max(df_for_plotting$fraction+df_for_plotting$SE) + 0.02), step.increase = 0.1, label = 'stars') +
          scale_x_discrete(limits = c("INTERFACE", "NOT_INTERFACE", "GWAS","NEGATIVE_CONTROL"), labels = c('I', 'NI', 'G', 'NC')) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
                plot.title = element_text(face = 'bold'), axis.text = element_text(face = 'bold', size = 14, color = 'black'), 
                axis.title = element_text(face = 'bold', size = 12, color = 'black'), legend.position = 'none') +
          scale_fill_manual(values = c("INTERFACE" = '#5975A3', "NOT_INTERFACE" = '#CC8963', "GWAS" = '#5F9E6E', "NEGATIVE_CONTROL" = '#B55D60')) + 
          labs(title = paste0("Enrichment for", " ", cell_types[cell]), x = "SNP TYPE", y = "FRACTION OF IMPORTANT GENES"))
  dev.off()
}