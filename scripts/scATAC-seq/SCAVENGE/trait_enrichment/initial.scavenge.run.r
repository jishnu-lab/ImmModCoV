# written by Prabal Chhibbar 
# Date: 13th October 2022 
# The code to perform the initial processing and reference before the permutation testing with different randomly selected seed cells 

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

set.seed(42)

args = commandArgs(trailingOnly = T)
# input args: args[1] - the trait file; args[2] - the filename for the network object 
# and args[3]: the filename for the seed cells object 
# and args[4]: the outfile name for the np scores 
# and args[5]: the outfile name for the z score matrix 

setwd("/ix/djishnu/Prabal/Aim2/covid/data/approach1b/scavenge") # the working directory where the analysis output is saved


# step1: preprocessing the summarised experiment object and saving it for the permutation analysis. this will be done for each trait file 

# the path to the trait file
load("peak_by_cell_matrix-summarizedexperiment.rda")
load("For_Prabal_covidpbmc-lsi-meta.rda")
#SE_pbmc5k = readRDS("scATAC_Heme_All_SummarizedExperiment.final.rds")
SE_pbmc5k = proj_PeakMatrix
#SE_pbmc5k = readRDS("scATAC-Healthy-Hematopoiesis-191120.rds")
trait_file = args[1]
#SE_pbmc5k = obj # remove this after the analysis
# gchromVAR analysis
assayNames(SE_pbmc5k) = 'counts'
#SE_pbmc5k = SE_pbmc5k[,SE_pbmc5k$Health_state=="D"] # subset the disease state
remove_cell_type = c("Undefined")

SE_pbmc5k = SE_pbmc5k[,!(SE_pbmc5k$BioCluster %in% remove_cell_type)]

# add the column with the total reads to be used as the sequencing depth

temp = data.frame(colData(SE_pbmc5k))
SE_pbmc5k$TotalReads = rowSums(temp[c("ReadsInBlacklist","ReadsInPromoter","ReadsInPromoter","ReadsInPeaks")]) # adding the reads in different regions together in order to get total reads


# filtering the peaks with not enough fragments
SE_pbmc5k$depth = SE_pbmc5k$TotalReads

SE_pbmc5k <- filterSamples(SE_pbmc5k,
                           min_in_peaks = 0.15, shiny = FALSE)
SE_pbmc5k <- filterPeaks(SE_pbmc5k, non_overlapping = TRUE)


SE_pbmc5k = addGCBias(SE_pbmc5k,genome=BSgenome.Hsapiens.UCSC.hg19)
#assayNames(SE_pbmc5k) = "counts"
SE_pbmc5k_bg = getBackgroundPeaks(SE_pbmc5k,niterations=200)
trait_import = importBedScore(rowRanges(SE_pbmc5k),trait_file,colidx=5)
SE_pbmc5k_DEV = computeWeightedDeviations(SE_pbmc5k,trait_import,background_peaks = SE_pbmc5k_bg )

# reformat the results

z_score_mat = data.frame(colData(SE_pbmc5k), z_score = t(assays(SE_pbmc5k_DEV)[["z"]]) %>% c)

## temp remove the rows with Nan Z-score values
#z_score_mat = z_score_mat %>% select(z_score_mat, z_score != "NaN")
z_score_mat1 = na.omit(z_score_mat)

z_score_mat1 <- z_score_mat[is.finite(z_score_mat$z_score),]

##
seed_idx = seedindex(z_score_mat1$z_score, 0.05)

# plot the bar plot of the proportions of cell chosen from the total number of the cells of that type

seed_numbers = data.frame(table(z_score_mat1[seed_idx,]$BioCluster))
total_numbers = data.frame(table(z_score_mat1$BioCluster))

colnames(seed_numbers) = c("celltype","freq")
colnames(total_numbers) = c("celltype","freq1")

df = merge(seed_numbers,total_numbers)
df[["prop"]] = df$freq/df$freq1
df = df[order(df$prop,decreasing = T),] # sort the data frame based on the proportion values

df$celltype = factor(df$celltype,levels = df$celltype) # done in order to prevent lexical sorting of the variable by ggplot
# calculate scale factor

scale_factor  = cal_scalefactor(z_score=z_score_mat1$z_score,0.01)
scale_factor
# scale factor is calculating from most enriched 1 percent of cells

# construct m-knn graph
# calculate tfidf-mat

peak_by_cell_mat <- assay(SE_pbmc5k)
tfidf_mat  = tfidf(bmat = peak_by_cell_mat, mat_binary = T,TF=T,log_TF = T)

# svd analysis of TF-IDF matrix

lsi_mat = do_lsi(tfidf_mat,dims=30)
mutualknn3 = getmutualknn(lsi_mat,30) # calculate the m-knn graph

# the knn network needs to be saved to be used at the later time 

saveRDS(mutualknn3,file=args[2]) 

# get the seed cells and save them too!

seed_cells = rownames(mutualknn3)[seed_idx]

# save the seed cells

saveRDS(seed_cells,file=args[3])

# perform the reference network propagation and save the scores to be used during permutations 

np_score = randomWalk_sparse(intM=mutualknn3,rownames(mutualknn3)[seed_idx],gamma=0.05)

saveRDS(np_score,file=args[4])

# save the z-score matrix 

saveRDS(z_score_mat,file=args[5])

print("Done!")
