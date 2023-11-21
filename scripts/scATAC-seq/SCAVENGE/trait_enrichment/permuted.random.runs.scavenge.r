##written by: Prabal Chhibbar 
##Date: october 13th 2022 
## perform the permutation test for quality of the network propagation results from scavenge by using randomly selected degree-matched seed cells multiple times 

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

#set.seed(args[5]) # seed is passed as an argument to ensure complete randomness and minimal overlaps in the parallel processes

args = commandArgs(trailingOnly = T)

# the input args: args[1]: the mutualknn3 network object filename, args[2]: the seed cells object ; args[3]: the original network propagation scores; args[4]: the output file to save the enrichment list object and args[5]: the seed for the script!
set.seed(args[5])
setwd("/ix/djishnu/Prabal/Aim2/covid/data/approach1b/scavenge")

# step 1: read in the network and the seed cells objects, the z-score matrix was also saved in the init script but is not needed to be read in this script 

mutualknn3 = readRDS(args[1])
seed_cells = readRDS(args[2])

# step2: make the degree based bins for the different nodes in the network for degree matched random sampling

# calculate the degree of each node in the cell-cell network

node_degrees = rowSums(mutualknn3)

# create the degree bins for degree matched random sampling for the permutation test!
bin1 = c() #degrees between 0 and 7
bin2 = c() # degrees between 7 and 15
bin3 = c() # degrees between 15 and 20
bin4 = c() # degrees between 20 and 25
bin5 = c() # degrees between 25 and 30

# sort all the nodes based on the degrees
for (name in names(node_degrees))
{
  deg = node_degrees[[name]]
  if (between(deg,0,7)) {bin1 = c(bin1,name)}
  else if(between(deg,7,15)) {bin2 = c(bin2,name)}
  else if(between(deg,15,20)) {bin3 = c(bin3,name)}
  else if(between(deg,20,25)) {bin4 = c(bin4,name)}
  else if(between(deg,25,30)) {bin5 = c(bin5,name)}
}
# remove all the informative seed cells from the sampling process

bin1 = setdiff(x=bin1,y=seed_cells)
bin2 = setdiff(x=bin2,y=seed_cells)
bin3 = setdiff(x=bin3,y=seed_cells)
bin4 = setdiff(x=bin4,y=seed_cells)
bin5 = setdiff(x=bin5,y=seed_cells)

# step3: perform the reference network propagation to be done to compare the results!

#np_score = randomWalk_sparse(intM=mutualknn3,rownames(mutualknn3)[seed_idx],gamma=0.05)
np_score = readRDS(args[3]) # the network propagation score from the reference run in the init step 
# step4: perform the permutation a 100 times per process - 10 processes in total, 10*100=1000

# initialize the list with zeros against the nodes in the vector
enrichment_list = vector("list",length(names(node_degrees)))
names(enrichment_list) = names(node_degrees)

for (name in names(enrichment_list)) {enrichment_list[[name]] = 0}

# perform the degree aware seed cell list permutations, n =10

for (i in 1:1000)
{
  if(i%%10==0){print(i)}
  random_seed_cells = c()
  #sample the degree matched random seed cells
  for (cell in seed_cells)
  {
    if (between(node_degrees[[cell]],0,7)){random_seed_cells = c(random_seed_cells,sample(bin1,1))}
    else if (between(node_degrees[[cell]],7,15)) {random_seed_cells = c(random_seed_cells, sample(bin2,1))}
    else if (between(node_degrees[[cell]],15,20)) {random_seed_cells = c(random_seed_cells,sample(bin3,1))}
    else if (between(node_degrees[[cell]],20,25)) {random_seed_cells = c(random_seed_cells,sample(bin4,1))}
    else if (between(node_degrees[[cell]],25,30)) {random_seed_cells =c(random_seed_cells,sample(bin5,1))}
  }
random_np_scores = randomWalk_sparse(intM=mutualknn3,random_seed_cells,gamma=0.05)
  for (name in names(random_np_scores)){print(name);if (random_np_scores[[name]]>=np_score[[name]]) {enrichment_list[[name]] = enrichment_list[[name]] +1  }}



}

# save the enrichment list to be collated in the end 

saveRDS(enrichment_list,file=args[4])

print("done!")

