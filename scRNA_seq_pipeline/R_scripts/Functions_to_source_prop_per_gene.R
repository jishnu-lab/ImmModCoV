###############FUNCTIONS TO USE FOR ESTIMATING CELL TYPE ENRICHMENT AT PER GENE LEVEL USING ALL SNPS###################
##########################################
##Create SNP score matrix##
##########################################
create_snp_score_matrix <- function(SNP_df) {
  colnames(SNP_df) <- c("CHR", "SNP", "GENE", "SCORE")
  #remove rows with no gene mapped
  SNP_df <- SNP_df[which(SNP_df$GENE != "NA"),]
  SNP_df$SNP_id <- paste0(SNP_df$CHR, ":", SNP_df$SNP)
  SNP_df <- SNP_df[, c("SNP_id", "GENE", "SCORE")]
  SNP_matrix <- matrix(ncol = length(unique(SNP_df$GENE)), nrow = length(unique(SNP_df$SNP_id)))
  colnames(SNP_matrix) <- unique(SNP_df$GENE)
  rownames(SNP_matrix) <- unique(SNP_df$SNP_id)
  for (i in 1:nrow(SNP_df)) {
    SNP_matrix[rownames(SNP_matrix) == SNP_df$SNP_id[i], 
               colnames(SNP_matrix) == SNP_df$GENE[i]] <- SNP_df$SCORE[i]
  }
  SNP_matrix[is.na(SNP_matrix)] <- as.numeric(0) 
  #replace all NAs with 1
  assign("SNP_mat", SNP_matrix, envir = .GlobalEnv)
}

########################################################################################################################################################################
#######PROP OF CELLS ABOVE THRESHOLD PERCENTILE METHOD WITHOUT EXPECTED COLLAPSED TO CELL TYPE PROCESS ALL SNP TYPES TOGETHER WITH COVID AND HEALTHY SEPARATELY#########
########################################################################################################################################################################

Gene_wise_enrichment_healthy_covid_sep <- function(SNP_matrix, RNA_seq, phenotype, plot_df, threshold, snp_df, cell_prop_threshold) {
  colnames(snp_df) <- c("CHR", "POS", "GENE", "SCORE", "SNP_type", "n_genes")
  snp_df <- snp_df[which(snp_df$GENE != 'NA'),]
  snp_df$snp_names <- paste0(snp_df$CHR, ':', snp_df$POS)
  index_order <- match(colnames(SNP_matrix), rownames(RNA_seq))
  #remove NA
  index_order <- index_order[which(index_order != 'NA')]
  RNA_seq <- RNA_seq[index_order,]
  index_order_snps <- match(rownames(RNA_seq), colnames(SNP_matrix))
  SNP_matrix <- SNP_matrix[,index_order_snps]
  observed_score <- (SNP_matrix %*% RNA_seq)
  ##SUBSET TO HEALTHY AND COVID CELLS##
  healthy_cell_names <- grep('Healthy', colnames(observed_score), ignore.case = T)
  phenotype_cell_names <- grep(phenotype, colnames(observed_score), ignore.case = T)
  observed_score_healthy <- observed_score[, healthy_cell_names]
  observed_score_phenotype <- observed_score[, phenotype_cell_names]
  observed_score_list <- c(observed_score_healthy, observed_score_phenotype)
  status <- c('healthy', phenotype)
  for (i in 1:length(status)) {
    score_mat <-  observed_score_list[[i]]
    percentile_cut_off <- quantile(score_mat, probs = c(threshold))
    prop_cells_per_snp <- data.frame('prop_cells' = apply(score_mat, 1, function(x) length(which(x >= percentile_cut_off))/length(x)))
    prop_cells_per_snp$snp_names <- rownames(prop_cells_per_snp)
    prop_cells_per_snp <- left_join(x = prop_cells_per_snp, y = snp_df, by = 'snp_names')
    prop_cells_per_snp <- prop_cells_per_snp[, c('snp_names', 'prop_cells', 'GENE', 'SNP_type', 'n_genes')]
    prop_cells_per_snp <- prop_cells_per_snp %>% group_by(GENE, SNP_type) %>% slice(which.max(prop_cells))
    prop_cells_per_snp$Status <- status[i]
    cell_prop_percentile_cut_off <- quantile(prop_cells_per_snp$prop_cells, probs = c(cell_prop_threshold))
    prop_cells_per_snp <- prop_cells_per_snp[prop_cells_per_snp$prop_cells >= cell_prop_percentile_cut_off,]
    plot_df <- rbind(plot_df, data.frame("snp_names" = prop_cells_per_snp$snp_names, "gene_names" = prop_cells_per_snp$GENE, "prop_cells" = prop_cells_per_snp$prop_cells, 'snp_type' = prop_cells_per_snp$SNP_type, 'status' = prop_cells_per_snp$Status, 'n_genes' = prop_cells_per_snp$n_genes))
  }
  return(plot_df)
}


########################################################################################################################################################################
#######PROP OF CELLS ABOVE THRESHOLD PERCENTILE METHOD WITHOUT EXPECTED COLLAPSED TO CELL TYPE PROCESS ALL CELL TYPES TOGETHER FOR ODDS RATIO#########
########################################################################################################################################################################

Gene_wise_enrichment_healthy_covid_sep_for_odds_ratio <- function(SNP_matrix, RNA_seq, phenotype, plot_df, threshold, snp_df, cell_prop_threshold) {
  colnames(snp_df) <- c("CHR", "POS", "GENE", "SCORE", "SNP_TYPE", "n_genes", "CELL_TYPE")
  snp_df <- snp_df[which(snp_df$GENE != 'NA'),]
  snp_df$snp_names <- paste0(snp_df$CHR, ':', snp_df$POS)
  index_order <- match(colnames(SNP_matrix), rownames(RNA_seq))
  #remove NA
  index_order <- index_order[which(index_order != 'NA')]
  RNA_seq <- RNA_seq[index_order,]
  index_order_snps <- match(rownames(RNA_seq), colnames(SNP_matrix))
  SNP_matrix <- SNP_matrix[,index_order_snps]
  observed_score <- (SNP_matrix %*% RNA_seq)
  ##SUBSET TO HEALTHY AND COVID CELLS##
  healthy_cell_names <- grep('Healthy', colnames(observed_score), ignore.case = T)
  phenotype_cell_names <- grep(phenotype, colnames(observed_score), ignore.case = T)
  observed_score_healthy <- observed_score[, healthy_cell_names]
  observed_score_phenotype <- observed_score[, phenotype_cell_names]
  observed_score_list <- c(observed_score_healthy, observed_score_phenotype)
  status <- c('healthy', phenotype)
  for (i in 1:length(status)) {
    score_mat <-  observed_score_list[[i]]
    percentile_cut_off <- quantile(score_mat, probs = c(threshold))
    prop_cells_per_snp <- data.frame('prop_cells' = apply(score_mat, 1, function(x) length(which(x >= percentile_cut_off))/length(x)))
    prop_cells_per_snp$snp_names <- rownames(prop_cells_per_snp)
    prop_cells_per_snp <- merge(x = prop_cells_per_snp, y = snp_df, by = 'snp_names')
    prop_cells_per_snp <- prop_cells_per_snp[, c('snp_names', 'prop_cells', 'GENE', 'CELL_TYPE', 'n_genes', 'SNP_TYPE')]
    prop_cells_per_snp <- prop_cells_per_snp %>% group_by(GENE, SNP_TYPE) %>% slice(which.max(prop_cells))
    prop_cells_per_snp$Status <- status[i]
    cell_prop_percentile_cut_off <- quantile(prop_cells_per_snp$prop_cells, probs = c(cell_prop_threshold))
    prop_cells_per_snp <- prop_cells_per_snp[prop_cells_per_snp$prop_cells >= cell_prop_percentile_cut_off,]
    plot_df <- rbind(plot_df, data.frame("snp_names" = prop_cells_per_snp$snp_names, "gene_names" = prop_cells_per_snp$GENE, "prop_cells" = prop_cells_per_snp$prop_cells, 'cell_type' = prop_cells_per_snp$CELL_TYPE, 'status' = prop_cells_per_snp$Status, 'n_genes' = prop_cells_per_snp$n_genes, 'snp_type' = prop_cells_per_snp$SNP_TYPE))
  }
  return(plot_df)
}

########calculate odds ratio using negative control####
odds_ratio_calc_neg_control <- function(unique_genes_df, phenotype, snp_type, odds_ratio_df, cell_type) {
  row_healthy <- which(unique_genes_df$snp_type == 'NEGATIVE_CONTROL')
  row_phenotype <- which(unique_genes_df$snp_type == snp_type)
  or <- (unique_genes_df$n[row_phenotype]*unique_genes_df$not_crossed[row_healthy])/(unique_genes_df$not_crossed[row_phenotype]*unique_genes_df$n[row_healthy])
  se<-sqrt((1/unique_genes_df$n[row_phenotype])+(1/unique_genes_df$not_crossed[row_healthy])+(1/unique_genes_df$not_crossed[row_phenotype])+(1/unique_genes_df$n[row_healthy]))
  upper_ci <- exp(log(or)+1.96*se)
  lower_ci <- exp(log(or)-1.96*se)
  z_stat <-(log(or)/se)
  p_val <- (exp(-0.717*z_stat-0.416*z_stat*z_stat))
  odds_ratio <- data.frame("ODDS_RATIO" = log(or),"LOWER_CI" = log(lower_ci), "UPPER_CI" = log(upper_ci), "P_VALUE" = p_val, 'CELL_TYPE' = cell_type, 'SNP_TYPE' = snp_type)
  odds_ratio_df <- rbind(odds_ratio_df, odds_ratio)
  return(odds_ratio_df)
}
########calculate odds ratio using negative control to compute difference in odds####
odds_ratio_calc <- function(unique_genes_df, phenotype, snp_type, odds_ratio_df, cell_type) {
  row_healthy <- which(unique_genes_df$snp_type == 'NEGATIVE_CONTROL')
  row_phenotype <- which(unique_genes_df$snp_type == snp_type)
  or <- (unique_genes_df$n[row_phenotype]*unique_genes_df$not_crossed[row_healthy])/(unique_genes_df$not_crossed[row_phenotype]*unique_genes_df$n[row_healthy])
  se<-sqrt((1/unique_genes_df$n[row_phenotype])+(1/unique_genes_df$not_crossed[row_healthy])+(1/unique_genes_df$not_crossed[row_phenotype])+(1/unique_genes_df$n[row_healthy]))
  odds_ratio <- data.frame("ODDS_RATIO" = or,"SE" = se,'CELL_TYPE' = cell_type, 'SNP_TYPE' = snp_type, 'SAMPLE' = phenotype)
  odds_ratio_df <- rbind(odds_ratio_df, odds_ratio)
  return(odds_ratio_df)
}
########Test the difference in odds ratio#######
diff_odds_ratio_test <- function(odds_ratio_df, cell_type, phenotype, odds_ratio_table_plotting) {
  snp_freq <-  odds_ratio_df %>% group_by(SNP_TYPE) %>% filter(n()>=2) 
  if(nrow(snp_freq) > 1) { #added if loop to avoid error when none of the snp types passed the threshold
    snp_types_in_odds <- unique(snp_freq$SNP_TYPE)
    odds_ratio_df <- odds_ratio_df[odds_ratio_df$SNP_TYPE %in% snp_types_in_odds,]
    odds_ratio_df_list <- split(odds_ratio_df, odds_ratio_df$SNP_TYPE)
    snp_types <- unique(odds_ratio_df$SNP_TYPE)
    snp_table <- data.frame('DELTA_OR' = numeric(), 'P_VAL' = numeric(), 'SNP_TYPE' = character(), 'CELL_TYPE' = character())
    for (snp in 1:length(snp_types)) {
      odds_ratio_table <- odds_ratio_df_list[[snp_types[snp]]]
      #testing if OR diff is diff from zero##
      delta_OR <- odds_ratio_table[odds_ratio_table$SAMPLE == phenotype,]$ODDS_RATIO - odds_ratio_table[odds_ratio_table$SAMPLE == 'healthy',]$ODDS_RATIO   
      delta_se <- sqrt((odds_ratio_table[odds_ratio_table$SAMPLE == phenotype,]$SE)^2 + (odds_ratio_table[odds_ratio_table$SAMPLE == 'healthy',]$SE)^2)
      stat <- delta_OR/delta_se
      p_val <- pnorm(stat, lower.tail = F)
      snp_table_type <- data.frame('DELTA_OR' = delta_OR, 'P_VAL' = p_val, 'SNP_TYPE' = snp_types[snp], 'CELL_TYPE' = cell_type)
      snp_table <- rbind(snp_table, snp_table_type)
    }
    odds_ratio_table_plotting <- rbind(odds_ratio_table_plotting, snp_table)
  } else {
    odds_ratio_table_plotting <- odds_ratio_table_plotting
  }
}
######Stars for p-val########
stars_pval <- function(x){
  stars <- c("***", "**", "")
  var <- c(0, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}
stars_pval_cell_type_bxplot <- function(x){
  stars <- c("***", "**", "NS")
  var <- c(0, 0.01, 0.05, 1)
  i <- findInterval(x, var, left.open = T, rightmost.closed = T)
  stars[i]
}

