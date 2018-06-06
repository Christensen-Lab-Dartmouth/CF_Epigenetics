###############################################################################################################
# Macrophage (MO) Methylation Heat Maps
# Script author: David Chen
# Date: 09/07/2017
###############################################################################################################

rm(list=ls())

library(ggplot2)
library(matrixStats)
library(pheatmap)
library(reshape2)
library(doParallel); registerDoParallel(detectCores() - 1)

library(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);
annot.850kb3 <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b3.hg19);

#-------------------------------------Data loading & annotation-------------------------------------
## Switch out barcode with sample IDs:
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

## Load & merge in RPMM class membership:
dat_sup <- read.csv("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/092517_RPMM/092517_RPMM_class_membership.csv")
colnames(dat_sup) <- gsub("class", "RPMM", colnames(dat_sup)); 
target_mo <- merge(target_mo, dat_sup[,c("Sample_Name","RPMM")], by="Sample_Name");

## Heat map annotations:
heat_annot <- data.frame(
  row.names = target_mo$Sample_Name,
  Status = target_mo$status,
  Lobe = target_mo$lobe,
  RPMM = target_mo$RPMM
);

row_annot <- data.frame(
  row.names = annot.850kb3$Name, 
  Island = ifelse(annot.850kb3$Relation_to_Island=="Island", "Yes", "No"),
  Promoter = ifelse(grepl("TSS",annot.850kb3$UCSC_RefGene_Group), "Yes", "No"),
  Enhancer = ifelse(annot.850kb3$X450k_Enhancer=="TRUE" | annot.850kb3$Phantom4_Enhancers != "" | annot.850kb3$Phantom5_Enhancers != "", "Yes", "No")
);

ann_colors <- list(
  Lobe = c(RUL="black", RLL="lightgray"),
  Status = c(CF="black",H="lightgray"),
  RPMM = c(rR="black",rL="lightgray"),
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  Island = c(Yes="black", No="lightgray")
);

#-------------------------------------Unsupervised heat map-------------------------------------
k <- 1e4; #most variable k CpGs:
sele <- order(rowVars(betas_lungMOs), decreasing=TRUE)[1:k];
mat <- betas_lungMOs[sele, ];
png("~/Downloads/Figure1A.png", res=300, units="in", height=8.27, width=11.69);
pheatmap(
  mat,
  show_rownames = FALSE, #CpGs
  show_colnames = TRUE, #samples
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors =  ann_colors,
  color = colorRampPalette(c("yellow", "blue"))(2048),
  border_color = NA, 
  fontsize = 13
);
dev.off();
rm(sele, mat, k);

#-------------------------------------DMPs between cases & controls-------------------------------------
## Model:
DMPs <- read.csv("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/030118_CF_limma_without_CellType_adjustment.csv", stringsAsFactors=F);
myCpGs <- DMPs$Name[DMPs$adj.P.Val < 0.05 & abs(DMPs$logFC) >= 3.50];
  
## Heatmap:
mat <- betas_lungMOs[rownames(betas_lungMOs) %in% myCpGs , ];
pheatmap(
  mat,
  show_rownames = FALSE, #CpGs
  show_colnames = TRUE, #samples
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors =  ann_colors,
  color = colorRampPalette(c("yellow", "blue"))(1024),
  border_color = NA,
  fontsize = 12
)
