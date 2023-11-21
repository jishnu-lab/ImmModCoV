## perform and test the scavenge pipeline
## written by: Prabal Chhibbar 
## This script is for getting the scores for each cell type in the dataset,outputs ameta file with information for each cell in the dataset 

library(SCAVENGE)
library(chromVAR)
library(gchromVAR)
library(BuenColors)
library(SummarizedExperiment)
library(data.table)
library(BiocParallel)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(igraph)
library(tidyverse)

set.seed(42) # setting the random seed 

args = commandArgs(trailingOnly = T)

inpath = args[1] # the path of the working directory 
trait_file = args[2] # the name of the bed file with the SNPs associated with the trait 
outfile = args[3] # the name of the output file where the meta information for each cell needs to be stored 

setwd(inpath) # the working directory where the analysis output is saved

# load the Rdata files with the knn graphs between the cells 

SE_pbmc5k = readRDS("scATAC_Heme_All_SummarizedExperiment.final.rds")

# subset based on the PBMCs and the healthy samples

cells_subset = c() 
names = c()


# remove the cells from the bone marrow and CD34 bone marrow donors 


for (i in 1:length(colnames(SE_pbmc5k))){
  name = strsplit(SE_pbmc5k$Internal_Name[[i]],"_")[[1]][[2]]
  if (!(name == "BM" | name == "CD34BM")){
    cells_subset = c(cells_subset,SE_pbmc5k$Internal_Name[[i]])
    
  }
}

SE_pbmc5k = SE_pbmc5k[,SE_pbmc5k$Internal_Name %in% cells_subset]

# remove the early progenitors cell types and unknown types from the summarised experiment 
remove_cell_types = c("PBMC_Rep3","PBMC_Rep2","BM_pDC","CMP","MEP","CD34_Progenitors_Rep2","CD34_Progenitors_Rep1","GMP","MPP","PBMC_Rep4","HSC","PBMC_Rep1","LMPP",
                      "CLP")

SE_pbmc5k  = SE_pbmc5k[,!(SE_pbmc5k$Group %in% remove_cell_types)]

# collapse the replicates into the same cell type 

cell_types = SE_pbmc5k$Group

cd4_rep1 = which(cell_types=="Memory_CD4_T_Cells_Rep2")
cd4_rep2 = which(cell_types=="Memory_CD4_T_Cells_Rep1")
cell_types[cd4_rep1] = "Memory_CD4_T_Cells"
cell_types[cd4_rep2] = "Memory_CD4_T_Cells"

naive_t1 = which(cell_types=="Naive_CD4_T_Cells_Rep2")
naive_t2 =  which(cell_types=="Naive_CD4_T_Cells_Rep1")
cell_types[naive_t1] = "Naive_CD4_T_Cells"
cell_types[naive_t2] = "Naive_CD4_T_Cells"

SE_pbmc5k$Group = cell_types

# gchromVAR analysis 

SE_pbmc5k = addGCBias(SE_pbmc5k,genome=BSgenome.Hsapiens.UCSC.hg19)
assayNames(SE_pbmc5k) = "counts"

# filtering the peaks with insufficient fragments

SE_pbmc5k <- filterSamples(SE_pbmc5k,
                                 min_in_peaks = 0.15, shiny = FALSE)
SE_pbmc5k <- filterPeaks(SE_pbmc5k, non_overlapping = TRUE)

SE_pbmc5k_bg = getBackgroundPeaks(SE_pbmc5k,niterations=200)
trait_import = importBedScore(rowRanges(SE_pbmc5k),trait_file,colidx=5)
SE_pbmc5k_DEV = computeWeightedDeviations(SE_pbmc5k,trait_import,background_peaks = SE_pbmc5k_bg )

# reformat the results 

z_score_mat = data.frame(colData(SE_pbmc5k), z_score = t(assays(SE_pbmc5k_DEV)[["z"]]) %>% c)

## temp remove the rows with Nan Z-score values 
#z_score_mat = z_score_mat %>% select(z_score_mat, z_score != "NaN")
z_score_mat = na.omit(z_score_mat)

z_score_mat <- z_score_mat[is.finite(z_score_mat$z_score),]

##

seed_idx = seedindex(z_score_mat$z_score, 0.01)
# plot the bar plot of the proportions of cell chosen from the total number of the cells of that type 

seed_numbers = data.frame(table(z_score_mat[seed_idx,]$Group))
total_numbers = data.frame(table(z_score_mat$Group))

colnames(seed_numbers) = c("celltype","freq")
colnames(total_numbers) = c("celltype","freq1")

df = merge(seed_numbers,total_numbers) 
df[["prop"]] = df$freq/df$freq1
df = df[order(df$prop,decreasing = T),] # sort the data frame based on the proportion values

# calculate scale factor 

scale_factor  = cal_scalefactor(z_score=z_score_mat$z_score,0.01)
scale_factor 

# scale factor is calculating from most enriched 1 percent of cells 

# construct m-knn graph 
# calculate tfidf-mat 

peak_by_cell_mat <- assay(SE_pbmc5k)
tfidf_mat  = tfidf(bmat = peak_by_cell_mat, mat_binary = T,TF=T,log_TF = T)

# svd analysis of TF-IDF matrix 

lsi_mat = do_lsi(tfidf_mat,dims=30)

mutualknn3 = getmutualknn(lsi_mat,30) # calculate the m-knn graph 

np_score = randomWalk_sparse(intM=mutualknn3,rownames(mutualknn3)[seed_idx],gamma=0.05)

# trait relevant score with scaled and normalized 

omit_idx = np_score == 0
sum(omit_idx)

mutualknn3 = mutualknn3[!omit_idx,!omit_idx]
np_score = np_score[!omit_idx]
TRS = np_score %>% capOutlierQuantile(.,0.95) %>% max_min_scale
TRS = TRS * scale_factor 

mono_mat = data.frame(z_score_mat[!omit_idx,], seed_idx[!omit_idx],np_score,TRS)
head(mono_mat)

print("writing out the meta files for further investigation!")
write.csv(mono_mat,outfile,row.names=F,col.names=F,sep="\t")

print("Done!")
