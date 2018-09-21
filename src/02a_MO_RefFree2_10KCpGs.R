# Reference-free deconvolution of macrophage EPIC data set using 20K CpGs
# Script author: David Chen
# Date: 09/25/17
# Notes:
# -- This script is for Linux Cluster computation

rm(list=ls());
library(matrixStats);
library(RefFreeEWAS);
library(doParallel); registerDoParallel(detectCores() - 1);
load("090817_Armstrong_MO_betas.RData");

## Select 10K with top sample variance:
set.seed(123);
n.CpG <- 10000;
Klist <- 2:10; 
iters <- 25; 

## DO NOT RUN: local testing
# RefFree_Array_MO <- RefFreeCellMixArray(Y_shortened, Klist=2, iters=1);
# pheatmap::pheatmap(RefFree_Array_MO[[1]]$Omega, border_color=NA, main="Matrix Mu");
# RefFree_Array_MO2 <- RefFreeCellMixArray(betas_lungMOs, Klist=2, iters=1, Yfinal=Y_shortened);
# RefFree_MO_Boots  <- RefFreeCellMixArrayDevianceBoots(RefFree_Array_MO2, Y_shortened, R=1, bootstrapIterations=1);
# RefFree_MO_Boots
## END

## Step 0: Define putative cell types, iterations, num. bootstraps:
sele <- order(rowVars(betas_lungMOs), decreasing=TRUE)[1:n.CpG]; 
Y_shortened <- betas_lungMOs[sele, ]; 

## Steps 1-2: Alternate fixing matrices Mu & Omega by iterating from 2-10 cell types:
RefFree_Array_MO <- RefFreeCellMixArray(
  Y_shortened,
  Klist = Klist,
  iters = iters
); 

print( sapply(RefFree_Array_MO, deviance, Y=Y_shortened) );

## Use the full beta matrix plus the k most variable probes to infer the cell types:
RefFree_Array_MO2 <- RefFreeCellMixArray(
  betas_lungMOs,
  Klist = Klist,
  iters = iters,
  Yfinal = Y_shortened
); 
print( sapply(RefFree_Array_MO2, deviance, Y=Y_shortened) );

## Export results:
rm(betas_lungMOs, target_mo, Klist, iters, n.CpG, sele);
save(
  list = ls(),
  file = "092417_RefFree_on_10K_Macrophage_EPIC.RData", 
  compress = TRUE
);
