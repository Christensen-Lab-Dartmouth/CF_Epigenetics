# Recursively Partitioned Mixture Model (RPMM)
# Script author: David Chen
# Date: 09/25/2017
# Notes:

rm(list=ls());
library(matrixStats);
library(RPMM);
library(doParallel); registerDoParallel(detectCores() - 1);
OUT_FILE_NAME <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/092517_RPMM/RPMM_cluster_membership_final.csv"; 

load_Yinv <- function(n.CpG=10000) {
  #'@description Load, subset, and transform data for RPMM
  load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/090817_Armstrong_MO_betas.RData");
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
  sele <- order(rowVars(betas_lungMOs), decreasing=TRUE)[1:n.CpG]; 
  Y_inv <- t(betas_lungMOs[sele, ]); 
  assign("Y_inv", Y_inv, envir=.GlobalEnv);
}

getRPMMClustLabels <- function(rpmmObject, Y_inv=NULL) {
  #'@description Extracts RPMM hard cluster labels
  #'@param rpmmObject RPMM object
  #'@param Y_inv Optional. Input matrix for RPMM computation. If provided, the sample name will be updated
  
  hardLabels <- blcTreeLeafClasses(rpmmObject);
  hardLabels <- as.data.frame(hardLabels);
  colnames(hardLabels) <- "RPMMClusters";
  if(! is.null(Y_inv)) rownames(hardLabels) <- hardLabels$Sample_Name <- rownames(Y_inv);
  return(hardLabels);
}

getRPMMSampOrder <- function(rpmmClusters, Y_inv) {
  #'@description Retrieves sample orders fo heat map visualization
  #'@describeIn Hinoue et al. online tutorial
  #'@param rpmmClusters data.frame with row.names = sample ID & 1 column named RPMM with cluster assignments
  #'@param Y_inv Input matrix for RPMM computation
  
  sampOrder <- c();
  for(r in names(table(rpmmClusters$RPMM))) {
    samps <- rownames(rpmmClusters)[rpmmClusters$RPMM == r];
    clu <- t(Y_inv[samps, ]);
    s_i <- seriation::seriate(clu, margin=2);
    so_i <- seriation::get_order(s_i);
    sampOrder <- c(sampOrder, samps[so_i]);
  }
  sampOrder <- data.frame(
    Sample_Name = sampOrder,
    RPMMSampleOrder = 1:length(sampOrder)
  );
  return(sampOrder);
}

main <- function() {
  ## RPMM with max_level of 2:
  load_Yinv();
  sup_rpmm <- blcTree(Y_inv, verbose=1, maxlevel=2-1); 
  print(sup_rpmm);
  
  ## Extract RPMM cluster labels:
  rpmmClusters <- getRPMMClustLabels(sup_rpmm, Y_inv);
  
  ## Retrieve RPMM sample order:
  sampOrders <- getRPMMSampOrder(rpmmClusters, Y_inv);  
  
  ## Export:
  dat_sup <- merge(rpmmClusters, sampOrders, by="Sample_Name");
  write.csv(dat_sup, file=OUT_FILE_NAME, row.names=FALSE, quote=FALSE);
}

main(); 
