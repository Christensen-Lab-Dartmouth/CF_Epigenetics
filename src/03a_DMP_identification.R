####################################################################################
# Probe-wise differential methylation analysis
# Script author: David Chen
# Date: 03/01/18
# Notes:
# 1. Do NOT adjust for cell types.
# 2. Account for paired nature.
####################################################################################

rm(list=ls())

library(limma)
library(matrixStats)
library(doParallel); registerDoParallel(detectCores() - 1)

#--------------------------------------Load data & switch out barcode with sample IDs--------------------------------------
## Annotation package:
library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);
annot.850kb3 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);

## Methylation data:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/090817_Armstrong_MO_betas.RData");

## Master covariate file:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/112017_MasterCovariates_BAL.RData");

#--------------------------------------Select most variable CpGs--------------------------------------
meth_vars <- sort(rowVars(betas_lungMOs), decreasing=TRUE);
k <- sum(meth_vars > 0.01);
sele <- order(rowVars(betas_lungMOs), decreasing=TRUE)[1:k];
betas_lungMOs <- betas_lungMOs[sele, ];
dim(betas_lungMOs)

#--------------------------------------Construct Design Matrix--------------------------------------
masterCovar$Status <- factor(masterCovar$Status, levels=c("H","CF"));
masterCovar$Subject <- as.factor(masterCovar$Subject);
masterCovar$Complete_Barcode <- as.character(masterCovar$Complete_Barcode); 

## Build design matrix and check
if( ! identical(colnames(betas_lungMOs), masterCovar$Complete_Barcode) ){
  print("Proceed to sample name matching...")
  betas_lungMOs <- betas_lungMOs[ , match(masterCovar$Complete_Barcode, colnames(betas_lungMOs))];
  mVals <- minfi::logit2(betas_lungMOs); 
  
  print("M-value matrix calculated! Proceed to model-matrix & subject-block setup...")
  myDesign <- model.matrix( ~ Status + Sex + Age, data=masterCovar);
  myBlock <- masterCovar$Subject; #block for calculating duplicate correlation for limma
  if(
    identical(masterCovar$Complete_Barcode, colnames(mVals)) &
    identical(masterCovar$Subject, myBlock)
  ){
    print("Checkpoint passed! OK to continue differential analysis")
  } else {
    stop("Something went wrong...")
  }
}

masterCovar[ , "Status"]
## H  CF H  H  H  H  CF CF H  H  H  CF CF CF CF CF
myDesign[ , "StatusCF"] #1=CF
## 0  1  0  0  0  0  1  1  0  0  0  1  1  1  1  1 

dupCor <- duplicateCorrelation(betas_lungMOs, myDesign, block=myBlock);
dupCor$consensus.correlation #larger DF, and more powerful
## [1] 0.7622167

#--------------------------------------Execute LIMMA--------------------------------------
stopifnot(identical( colnames(betas_lungMOs), colnames(mVals) ) ) #important checkpoint

## Execute statistical tests:
fit <- lmFit(
  mVals,
  design = myDesign, 
  correlation = dupCor$consensus.correlation, 
  block = myBlock
);
fit2 <- eBayes(fit);

## Subset annotation & identify:
ann850KSub <- annot.850kb3[match(rownames(mVals),annot.850kb3$Name), ];
DMPs <- topTable(
  fit2,
  number = Inf, 
  coef = "StatusCF", 
  genelist = ann850KSub, 
  adjust.method = "fdr", 
  sort.by = "p"
);

#--------------------------------------Gene List--------------------------------------
thresh <- 0.05;
mean(DMPs$adj.P.Val < thresh)
sum(DMPs$adj.P.Val < thresh)
DMPs$isSignifAt0.05 <- DMPs$adj.P.Val < thresh;
DMPs$isSignifAt0.10 <- DMPs$adj.P.Val < thresh*2;

#--------------------------------------Export--------------------------------------
## Export relevant objects as a workspace:
setwd("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/")
save(
  list = c("DMPs","dupCor","fit","fit2","myDesign","ann850KSub"),
  file = "030118_DMP_mostVariable26K_NoCellType.RData",
  compress = TRUE
)
