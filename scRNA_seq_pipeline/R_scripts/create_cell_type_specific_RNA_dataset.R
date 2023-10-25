##Process single cell RNA-seq matrices##
##Normalize per cell##
##Label the names by statusand split into cell-type objects to make analysis easier##
##Works only if status is set to COVID and Healthy##
library(R.utils)
library(tidyr)
library(Seurat)
library(SeuratDisk)
library(edgeR)
library(Matrix.utils)
library(stringr)
library(foreach)
library(doParallel)
args <- commandArgs(trailingOnly = TRUE) 
print(args)
RNA_seq_dir <- args[1] #directory with rna seq object
RNA_seq_file <- args[2] #rna seq object
output_dir <- args[3]
phenotype <- args[4]
output_sub_dir <- dir.create(path=paste0(output_dir, '/', str_remove(RNA_seq_file, '.rds')))
##Note##
##RNA seq dataset needs to be in rds, cell type annotation column to be used should be named cell.type, sample label column should be orig.ident, and covid status column to be Status##
RNA_seq <- readRDS(paste0(RNA_seq_dir, '/', RNA_seq_file))
cell_types <- unique(RNA_seq@meta.data$cell.type)
for (cell in 1:length(cell_types)) {
  print(paste0("Processing", " ", cell_types[cell]))
  #subset the data set to only cells of this type 
  rna_seq_cell_specific <- subset(RNA_seq, subset = (cell.type == cell_types[cell]))
  healthy_cells <- length(grep(rna_seq_cell_specific@meta.data$Status, pattern = 'healthy', ignore.case = T))
  phenotype_cells <- length(grep(rna_seq_cell_specific@meta.data$Status, pattern = phenotype, ignore.case = T))
  if (nrow(rna_seq_cell_specific@meta.data) >= 1000 & healthy_cells >= 100 & phenotype_cells >= 100) { 
    count_matrix <- rna_seq_cell_specific@assays$RNA@counts
    metadata <- rna_seq_cell_specific@meta.data
    groups <- metadata[, c("orig.ident", "Status")]
    groups$identity <- paste0(rownames(groups), "_", cell_types[cell], "_", groups$Status)
    colnames(count_matrix) <- groups$identity
    #CPM normalize matrix 
    count_matrix <- NormalizeData(count_matrix)
    #save the dataset
    saveRDS(count_matrix, file = paste0(output_dir, '/', str_remove(RNA_seq_file, '.rds'), '/', cell_types[cell], '_', 'count_matrix.rds'))
  } else {
    print(paste0("Not enough cells for ", cell_types[cell], " in ", RNA_seq_file))
  }
}