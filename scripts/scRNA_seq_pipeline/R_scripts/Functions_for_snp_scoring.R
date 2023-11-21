#Functions for scoring SNP files#
##########################################
#Scoring GWAS SNPs#
##########################################
score_gwas_snps <- function(snp_df) {
  snp_df <- snp_df %>% mutate(Rank = case_when(ANNO == 'promoter' | ANNO == '2kb' ~ 1, ANNO == 'Roadmap_enhancer' | ANNO == 'intergenic' | ANNO == 'genic' ~ 2, ANNO == 'gene_body' ~ 3))
  snp_df <- snp_df[order(snp_df$Rank, decreasing = TRUE), ]
  snp_df <- snp_df[!duplicated(snp_df$GENE), ] 
  df <- data.frame("SNP" = character(), "GENE" = character(), "ANNO" = character(), "SCORE" = numeric())
  SNP <- unique(snp_df$SNP)
  promoter_n <- length(which(snp_df$ANNO == 'promoter' | snp_df$ANNO == '2kb'))
  enhancer_n <- length(which(snp_df$ANNO == 'Roadmap_enhancer' | snp_df$ANNO == 'intergenic' | snp_df$ANNO == 'genic'))
  gene_body_n <- length(which(snp_df$ANNO == 'gene_body'))
  score_split <- 1/((3*promoter_n) + (2*enhancer_n) + (1*gene_body_n))
  for (i in 1:nrow(snp_df)){
    gene <- snp_df$GENE[i]
    anno <- snp_df$ANNO[i]
    if(anno == 'Roadmap_enhancer' | anno == 'intergenic' | anno == 'genic') {
      score <- 2*score_split
    } else if(anno == 'promoter' | anno == '2kb') {
      score <- 3*score_split
    } else {
      score <- 1*score_split
    }
    df_snp <- data.frame("SNP" = SNP, "GENE" = gene, "ANNO" = anno, "SCORE" = score)
    df <- rbind(df, df_snp)
  }
  snp_specific_df <- return(df)
}