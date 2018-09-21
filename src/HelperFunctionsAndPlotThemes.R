# Helper Functions and Plot Themes
# Script author: David Chen
# Date: (ongoing)
# Notes:
# -- This script is to be imported/sourced for multiple analyses.

#-------------------------------------------------- Data/Annotation Loading Methods --------------------------------------------------
loadEPICandSampleAnnotation <- function(dnam_path="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/090817_Armstrong_MO_betas.RData",
                                        rpmm_path="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/092517_RPMM/RPMM_cluster_membership_final.csv",
                                        bySubject=FALSE) {
  #'@description Loads EPIC data and update annotation, then send each to global environment
  #'@param dnam_path RData path for EPIC data and annotation
  #'@param rpmm_path CSV path for RPMM reuslt
  #'@param bySubject Should Sample_Name be named based on subject (TRUE) or sample (FALSE)?
  
  load(dnam_path);
  dat_sup <- read.csv(rpmm_path, stringsAsFactors=FALSE);
  
  ## Match data & annotation:
  target_mo$Sample_Name <- paste("Sample", target_mo$Sample_Name, sep="_");
  target_mo$Complete_Barcode <- paste(target_mo$Slide, target_mo$Array,sep="_");
  target_mo <- merge(target_mo, dat_sup, by="Sample_Name");
  
  if(bySubject) target_mo$Sample_Name <- paste0(target_mo$status, target_mo$matched_pair, "_", target_mo$lobe);
  if(all(target_mo$Complete_Barcode == colnames(betas_lungMOs))){ #checkpoint
    print("Samples already matched, so proceed to name switching...");
    colnames(betas_lungMOs) <- target_mo$Sample_Name; 
  } else {
    print("Matching sample names...");
    betas_lungMOs <- betas_lungMOs[ , match(target_mo$Complete_Barcode, colnames(betas_lungMOs))]; 
    print("Now ready to name switching...");
    colnames(betas_lungMOs) <- target_mo$Sample_Name; 
  }
  
  ## Return both EPIC data & sample annotation to global environment:
  assign("betas_lungMOs", betas_lungMOs, envir=.GlobalEnv);
  assign("target_mo", target_mo, envir=.GlobalEnv);
}

loadEPICannot <- function(forHeatmap=TRUE) {
  #'@description Load EPIC annotation for heatmap
  #'@param forHeatmap If TRUE, a data.frame specific for `pheatmap` will be returned. Defaults to TRUE
  
  require(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);
  annot.850kb3 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);
  if(! forHeatmap) return(annot.850kb3);
  
  row_annot <- data.frame(
    row.names = annot.850kb3$Name, 
    Island = ifelse(annot.850kb3$Relation_to_Island=="Island", "Yes", "No"),
    Promoter = ifelse(grepl("TSS",annot.850kb3$UCSC_RefGene_Group), "Yes", "No"),
    Enhancer = ifelse(annot.850kb3$X450k_Enhancer=="TRUE" | annot.850kb3$Phantom4_Enhancers != "" | annot.850kb3$Phantom5_Enhancers != "", "Yes", "No")
  );
  return(row_annot);
}


#-------------------------------------------------- Exploratory & Analytic Methods --------------------------------------------------
getMostVariableSites <- function(data, k) {
  #'@description Get most variable CpGs
  #'@param data DNA methylation data with row, column = CpGs, samples
  #'@param k Number of most variable CpGs to select
  
  sele <- order(rowVars(data), decreasing=TRUE)[1:k];
  mat <- data[sele, ];
  return(mat);
}


#--------------------------------------------------Plot Themes--------------------------------------------------
require(ggplot2);
myVolcanoTheme <- theme_classic() +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_text(size=21, color="black"),
        legend.title=element_blank(), legend.text=element_blank(), legend.position="none"); 

myBoxplotTheme <- theme_classic() +
  theme(axis.text.x=element_text(size=20,color="black"), axis.title.x=element_blank(),
        axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
        legend.text=element_text(size=20,color="black",face="bold"), legend.title=element_blank(), legend.position="top");
