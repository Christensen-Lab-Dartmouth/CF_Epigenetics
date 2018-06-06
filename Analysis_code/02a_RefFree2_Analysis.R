####################################################################################
# Reference-free deconvolution of EPIC data set
# Script author: David Chen
# Date: 09/19/17
# Notes:
####################################################################################

rm(list=ls())

library(matrixStats)
library(RefFreeEWAS)
library(doParallel); registerDoParallel(detectCores() - 1)

load("090817_Armstrong_MO_betas.RData")
# load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/090817_Armstrong_MO_betas.RData") #local testing
target_mo$Sample_Name <- paste("Sample", target_mo$Sample_Name, sep="_");
target_mo$Complete_Barcode <- paste(target_mo$Slide, target_mo$Array,sep="_");
if(all(target_mo$Complete_Barcode == colnames(betas_lungMOs))){ #checkpoint
  print("Samples already matched, so proceed to name switching...");
  colnames(betas_lungMOs) <- target_mo$Sample_Name; 
} else {
  print("Matching sample names...");
  betas_lungMOs <- betas_lungMOs[ , match(target_mo$Complete_Barcode, colnames(betas_lungMOs))]; 
  print("Now ready to name switching...");
  colnames(betas_lungMOs) <- target_mo$Sample_Name; 
}

## Select 10,000 CpGs with top sample variance:
set.seed(123);
n.CpG <- 10000;
sele <- order(rowVars(betas_lungMOs), decreasing=TRUE)[1:n.CpG]; #most variable across all samples

## Due to error in matrix inversion: separately analyze each lobe (i.e., unique patients per analysis)
rll.samps <- target_mo$Sample_Name[target_mo$lobe=="RLL"];
rul.samps <- target_mo$Sample_Name[target_mo$lobe=="RUL"];

## Shortened:
Y_shortened.RLL <- betas_lungMOs[sele, colnames(betas_lungMOs) %in% rll.samps ]; 
Y_shortened.RUL <- betas_lungMOs[sele, colnames(betas_lungMOs) %in% rul.samps]; 

## Full subsets:
Y_RLL <- betas_lungMOs[, colnames(betas_lungMOs) %in% rll.samps ]; 
Y_RUL <- betas_lungMOs[, colnames(betas_lungMOs) %in% rul.samps]; 

## Define parameters:
Klist <- 2:8; #update 
iters <- 25;
bootstrapIterations <- 10; 
R <- 1500; #more bootstraps may be necessary

#----------------------------------------RefFreeEWAS 2.0: RLL----------------------------------------
## Steps 1-2: Alternate fixing matrices Mu & Omega by iterating from 2-10 cell types:
RefFree_Array_RLL <- RefFreeCellMixArray(
  Y_shortened.RLL,
  Klist = Klist,
  iters = iters
); 
print( sapply(RefFree_Array_RLL, deviance, Y=Y_shortened.RLL) )

## Use the full beta matrix plus the k most variable probes to infer the cell types:
RefFree_Array_RLL2 <- RefFreeCellMixArray(
  Y_RLL,
  Klist = Klist,
  iters = iters,
  Yfinal = Y_shortened.RLL
); 
print( sapply(RefFree_Array_RLL2, deviance, Y=Y_shortened.RLL) )

## Step 3: Determine the optimal number of putative cell types K by bootstrapping:
RefFree_RLL_Boots <- RefFreeCellMixArrayDevianceBoots(
  RefFree_Array_RLL2,
  Y_shortened.RLL,
  R = R, #more bootstraps may be necessary
  bootstrapIterations = bootstrapIterations
)

## Step 4: Export
save(
  list = c("Y_shortened.RLL","RefFree_Array_RLL","RefFree_Array_RLL2","RefFree_RLL_Boots"),
  file = "092917a_RefFree_on_EPIC_RLL.RData", 
  compress = TRUE
);

print("RLL complete and exported!")

## Remove saved objects for RLL to avoid confusion/easy mistake:
rm(Y_shortened.RLL, RefFree_Array_RLL, RefFree_Array_RLL2, RefFree_RLL_Boots); 

#----------------------------------------RefFreeEWAS 2.0: RUL----------------------------------------
## Steps 1-2: Alternate fixing matrices Mu & Omega by iterating from 2-10 cell types:
RefFree_Array_RUL <- RefFreeCellMixArray(
  Y_shortened.RUL, #double check
  Klist = Klist,
  iters = iters
); 
print( sapply(RefFree_Array_RUL, deviance, Y=Y_shortened.RUL) )

## Use the full beta matrix plus the k most variable probes to infer the cell types:
RefFree_Array_RUL2 <- RefFreeCellMixArray(
  Y_RUL,
  Klist = Klist,
  iters = iters,
  Yfinal = Y_shortened.RUL
); 
print( sapply(RefFree_Array_RUL2, deviance, Y=Y_shortened.RUL) )

## Step 3: Determine the optimal number of putative cell types K by bootstrapping:
RefFree_RUL_Boots <- RefFreeCellMixArrayDevianceBoots(
  RefFree_Array_RUL2,
  Y_shortened.RUL, 
  R = R, 
  bootstrapIterations = bootstrapIterations
);

save(
  list = list("Y_shortened.RUL","RefFree_Array_RUL","RefFree_Array_RUL2","RefFree_RUL_Boots"),
  file = "092917b_RefFree_on_EPIC_RUL.RData", 
  compress = TRUE
);

print("RUL complete and exported!")
print("Analysis complete!")
