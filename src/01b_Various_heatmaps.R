# Macrophage (MO) Methylation Heat Maps
# Script author: David Chen
# Date: (ongoing)
# Notes:

rm(list=ls());
library(ggplot2);
library(matrixStats);
library(pheatmap);
library(reshape2);
library(doParallel); registerDoParallel(detectCores() - 1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsAndPlotThemes.R");

nCpGs <- 1e4; #most variable k CpGs
MVAL_THRESH <- 3.50; #asked by reviewer
FDR_THRESH <- 0.05;
HEAT_COLORS <- colorRampPalette(c("yellow","black","blue"))(1024); 

loadEPICandSampleAnnotation(bySubject=TRUE);

#-------------------------------------Heat map: unsupervised clustering-------------------------------------
heat_annot <- data.frame(
  row.names = target_mo$Sample_Name,
  Status = target_mo$status,
  Lobe = target_mo$lobe,
  RPMM = target_mo$RPMMClusters
);

row_annot <- loadEPICannot();

ann_colors <- list(
  Lobe = c(RUL="black", RLL="lightgray"),
  Status = c(CF="black", H="lightgray"),
  RPMM = c(rR="black", rL="lightgray"),
  Promoter = c(Yes="black", No="lightgray"),
  Enhancer = c(Yes="black", No="lightgray"),
  Island = c(Yes="black", No="lightgray")
);

mat <- getMostVariableSites(betas_lungMOs, k=nCpGs);
sampOrders <- target_mo$Sample_Name[order(target_mo$RPMMSampleOrder, decreasing=TRUE)];
png("~/Downloads/Figure1A.png", res=300, units="in", height=8.27, width=11.69);
pheatmap(
  mat[ , sampOrders],
  show_rownames = FALSE, #CpGs
  show_colnames = TRUE, #samples
  cluster_cols = FALSE, 
  gaps_col = 5,
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  color = HEAT_COLORS,
  border_color = NA, 
  fontsize = 13
);
dev.off();

## Fisher's exact test:
contTab <- table(target_mo$RPMMClusters, target_mo$status);
contTab <- contTab[c(2,1), ];
print(contTab);
fisher.test(contTab);

#-------------------------------------Heat map for DMPs between cases & controls-------------------------------------
## Model output:
DMP_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/030118_CF_limma_without_CellType_adjustment.csv"; 
DMPs <- read.csv(DMP_PATH, stringsAsFactors=FALSE);
myCpGs <- DMPs$Name[DMPs$adj.P.Val < FDR_THRESH & abs(DMPs$logFC) >= MVAL_THRESH];

## Heatmap:
mat <- betas_lungMOs[rownames(betas_lungMOs) %in% myCpGs, ];
png("~/Downloads/Figure_DMPs_log2FC_over3.png", res=300, units="in", height=8.27, width=11.69);
pheatmap(
  mat,
  show_rownames = FALSE, #CpGs
  show_colnames = TRUE, #samples
  annotation_col = heat_annot,
  annotation_row = row_annot,
  annotation_colors = ann_colors,
  color = HEAT_COLORS,
  border_color = NA,
  fontsize = 12
);
dev.off();
