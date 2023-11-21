##Compute difference in odds ratio and plot both delta OR and negative log-p-val##
library(R.utils)
library(tidyr)
library(Matrix.utils)
library(stringr)
library(dplyr)
library(rstatix)
library(ggplot2)
library(ggpubr)
library(sparseMatrixStats)
library(epitools)
library(ggprism)
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
odds_ratio_table_plotting <- data.frame('DELTA_OR' = numeric(), 'P_VAL' = numeric(), 'SNP_TYPE' = character(), 'CELL_TYPE' = character())
for (cell in 1:length(cell_types)){
  print(cell_types[cell])
  cell_type_RNA_seq <- readRDS(paste0(RNA_seq_dir, '/', RNA_seq_folder, '/', RNA_seq_files[cell]))
  plot_df <- data.frame("snp_names" = character(), "gene_names" = character(),  "prop_cells" = numeric(), 'cell_type' = character(),  'status' = character(), 'n_genes' = numeric(), snp_type = character())
  SNP_df_for_plotting <- data.frame("CHR" = character(), "POS" = integer(), "GENE" = character(), "SCORE" = numeric(), "SNP_type" = character(), "n_genes" = numeric())
  SNP_df <- data.frame("CHR" = character(), "POS" = integer(), "GENE" = character(), "SCORE" = numeric())
  odds_ratio_df <- data.frame("ODDS_RATIO" = numeric(),"SE" = numeric(),'CELL_TYPE' = character(), 'SNP_TYPE' = character(), 'SAMPLE' = character())
  for (snp in 1:length(SNP_files)) {
    SNP_df_type <- read.table(paste0(SNP_dir, '/', SNP_dataset, '/', SNP_files[snp]))
    SNP_df_type_plotting <- SNP_df_type
    n_genes <- length(unique(SNP_df_type$V3))
    SNP_type <- SNP_types[snp]
    SNP_df_type_plotting$SNP_type <- SNP_type
    SNP_df_type_plotting$n_genes <- n_genes
    SNP_df_type_plotting$CELL_TYPE <- cell_types[cell]
    SNP_df <- rbind(SNP_df, SNP_df_type)
    SNP_df_for_plotting <- rbind(SNP_df_for_plotting, SNP_df_type_plotting)
  } 
  create_snp_score_matrix(SNP_df = SNP_df)
  plot_df <- Gene_wise_enrichment_healthy_covid_sep_for_odds_ratio(SNP_matrix = SNP_mat, RNA_seq = cell_type_RNA_seq, phenotype = phenotype, threshold = 0.95, snp_df = SNP_df_for_plotting, cell_prop_threshold = 0.95, plot_df = plot_df)
  unique_genes <- plot_df %>% group_by(snp_type, status) %>% summarise(n=n()) 
  samples <- c('healthy', phenotype)
  for (sample in 1:length(samples)) {
    unique_genes_sample <- unique_genes[which(unique_genes$status == samples[sample]),]
    unique_genes_snp_types <- unique_genes_sample$snp_type
    unique_genes_list <- split(unique_genes_sample, unique_genes_sample$snp_type)
    n_genes_per_type <- plot_df[, c("snp_type", "n_genes")] %>% unique()
    snp_types <- unique(n_genes_per_type$snp_type)
    snp_types <- snp_types[which(snp_types != 'NEGATIVE_CONTROL')]
    for (snp in 1:length(snp_types)) {
      if (snp_types[snp] %in% unique_genes_snp_types) {
        unique_genes_df <- unique_genes_list[[snp_types[snp]]]
        unique_genes_df <- rbind(unique_genes_df, unique_genes_list[["NEGATIVE_CONTROL"]])
        n_genes_snp <- n_genes_per_type[n_genes_per_type$snp_type == snp_types[snp], 2]
        n_genes_neg_control <- n_genes_per_type[n_genes_per_type$snp_type == "NEGATIVE_CONTROL", 2]
        n_genes <- c(n_genes_snp, n_genes_neg_control)
        unique_genes_df$all_genes <- n_genes
        unique_genes_df$not_crossed <- (unique_genes_df$all_genes - unique_genes_df$n)
        unique_genes_df <- unique_genes_df[, colnames(unique_genes_df) != 'all_genes']
        odds_ratio_df <- odds_ratio_calc(unique_genes_df = unique_genes_df, phenotype = samples[sample], snp_type = snp_types[snp], odds_ratio_df = odds_ratio_df, cell_type = cell_types[cell])
      } else {
        print(paste0("No significant genes for", " ", snp_types[snp]))
      }
    }
  }
  odds_ratio_table_plotting <- diff_odds_ratio_test(odds_ratio_df = odds_ratio_df, cell_type = cell_types[cell], phenotype = phenotype, odds_ratio_table_plotting = odds_ratio_table_plotting)
}

if(dir.exists(paste0(output_dir, '/', RNA_seq_folder)) == FALSE) {
  dir.create(paste0(output_dir, '/', RNA_seq_folder))
  print("RNA directory does not exist")} else {
    print("RNA directory exists")
  }
if(dir.exists(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset)) == FALSE) {
  dir.create(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset))} else {
    print("SNP directory exists")
  }
odds_ratio_snp_type_list <- split(odds_ratio_table_plotting, odds_ratio_table_plotting$SNP_TYPE)
colors <- data.frame('snp_type' = c('INTERFACE', 'GWAS', 'NOT_INTERFACE'), p_val_color = c('#107010', '#2020DF', '#e7861b'))
for (snp in 1:length(snp_types)) {
  print(snp_types)
  odds_ratio_snp_type <- odds_ratio_snp_type_list[[snp]]
  odds_ratio_snp_type <- odds_ratio_snp_type[!is.infinite(odds_ratio_snp_type$DELTA_OR),]
  odds_ratio_snp_type <- odds_ratio_snp_type[!is.na(odds_ratio_snp_type$DELTA_OR),]
  SNP <- unique(odds_ratio_snp_type$SNP_TYPE)
  if(length(cell_types) > length(unique(odds_ratio_snp_type$CELL_TYPE))) {
    missing_cell_types <- cell_types[!cell_types %in% unique(odds_ratio_snp_type$CELL_TYPE)]
    missing_cell_types_df <- data.frame('DELTA_OR' = rep(0, length(missing_cell_types)), 'P_VAL' = rep(1, length(missing_cell_types)), 'SNP_TYPE' = rep(SNP, length(missing_cell_types)), 'CELL_TYPE' = missing_cell_types)
    odds_ratio_snp_type <- rbind(odds_ratio_snp_type, missing_cell_types_df)
  } else {
    print(paste0("All cell types present for ", SNP))
  }
  colors_df <- colors[colors$snp_type %in% SNP, ]
  odds_ratio_snp_type <- odds_ratio_snp_type[order(odds_ratio_snp_type$DELTA_OR, decreasing = TRUE),]
  odds_ratio_snp_type <- transform(odds_ratio_snp_type, stars = stars_pval(odds_ratio_snp_type$`P_VAL`))
  #rownames(odds_ratio_snp_type) <- odds_ratio_snp_type$CELL_TYPE
  #odds_ratio_snp_type_t <- data.frame(t(odds_ratio_snp_type))
  cols_to_keep <- c('DELTA_OR', 'P_VAL', 'CELL_TYPE')
  df_for_plot <- odds_ratio_snp_type[ ,colnames(odds_ratio_snp_type) %in% cols_to_keep]
  plot_stars <- data.frame('group1' = odds_ratio_snp_type$CELL_TYPE, 'group2' = odds_ratio_snp_type$CELL_TYPE, 'p.adj' = odds_ratio_snp_type$stars, y.position = odds_ratio_snp_type$DELTA_OR)
  pdf(paste0(output_dir, '/', RNA_seq_folder, '/', SNP_dataset, '/', SNP, '_', 'diff_odds_ratio.pdf'), height = 10, width = 15)
  print(ggplot(data = df_for_plot, aes(x = CELL_TYPE, y = DELTA_OR)) +
          geom_bar(stat="identity", position = position_dodge(), fill = colors_df$p_val_color, color = 'black', size = 0.75) +
          geom_hline(yintercept=0, 
                     color = "black") +
          scale_x_discrete(guide = guide_axis(angle = 90), limits = odds_ratio_snp_type$CELL_TYPE) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black", size = 0.75), 
                plot.title = element_text(face = 'bold'), axis.text = element_text(face = 'bold', size = 16, color = 'black'), 
                axis.title = element_text(face = 'bold', size = 16, color = 'black')) +
          labs(title = paste0("CELL TYPE ENRICHMENT FOR", " ", SNP), x = "CELL TYPE", y = "DELTA OR", size = 16) +
          add_pvalue(plot_stars, bracket.size = 0, fontface = 'bold', label.size = 8, tip.length = 0)
  )
  dev.off()
}