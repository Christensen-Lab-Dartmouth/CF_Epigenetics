# MethylationEPIC RefFree Continuation and Data Visualization
# Script author: David Chen
# Date: 09/24/2017
# Notes:  
# -- This is continuation of script tagged "02a_", but run locally instead of Linux

rm(list=ls());
library(ggplot2);
library(matrixStats);
library(pheatmap);
library(RefFreeEWAS);
library(reshape2);
library(doParallel); registerDoParallel(detectCores() - 1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsAndPlotThemes.R");
DIR <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/";
setwd(DIR);

## Decide Optimal K:
load("091917_MO_EPIC/092417_RefFree_on_10K_Macrophage_EPIC.RData");
RefFree_MO_Boots <- RefFreeCellMixArrayDevianceBoots(RefFree_Array_MO2, Y_shortened); #may need several tries
mean_Ks <- apply(RefFree_MO_Boots[-1, ], 2, mean, trim=0.25);
optimalK <- names(mean_Ks)[which.min(mean_Ks)];
optimalK

## Retrieve Cell Type Proportions:
myOmega <- RefFree_Array_MO[[optimalK]]$Omega;
loadEPICandSampleAnnotation(bySubject=TRUE);
if(! identical(target_mo$Complete_Barcode, rownames(myOmega))) {
  myOmega <- myOmega[match(target_mo$Complete_Barcode, rownames(myOmega)), ];
  stopifnot(identical(target_mo$Complete_Barcode, rownames(myOmega))); #checkpoint
  rownames(myOmega) <- target_mo$Sample_Name;
}
rownames(myOmega) <- target_mo$Sample_Name;

## Heat map annotations:
heat_annot <- data.frame(
  row.names = target_mo$Sample_Name,
  Status = target_mo$status,
  Lobe = target_mo$lobe,
  RPMM = target_mo$RPMMClusters
);

ann_colors <- list(
  Lobe = c(RUL="black", RLL="lightgray"),
  Status = c(CF="red",H="blue"),
  RPMM = c(rR="black",rL="lightgray")
);

fixInNamespace("draw_colnames","pheatmap"); #vjust = 1, hjust = 0.5, rot = 0
png("~/Downloads/Figure3A.png", res=300, units="in", height=8.27, width=5.84);
pheatmap(
  myOmega,
  labels_col = paste("Cell type", colnames(myOmega)),
  color = colorRampPalette(c("white", "darkorchid"))(1024),
  annotation_row = heat_annot, #samples
  annotation_colors = ann_colors,
  fontsize = 12,
  border_color = NA
)
dev.off();

## Box plots:
plt.myOmega <- myOmega; #copy
colnames(plt.myOmega) <- paste0("CellType",colnames(plt.myOmega)); 
plt.myOmega <- melt(plt.myOmega, varnames=c("Sample_Name","CellType"));
plt.myOmega <- merge(plt.myOmega, target_mo[ , c("Sample_Name","status","lobe")], all.x=TRUE); #merge in values
plt.myOmega$status <- gsub("CF", "Cystic fibrosis (CF)", plt.myOmega$status);
plt.myOmega$status <- gsub("H", "Healthy (H)", plt.myOmega$status);
plt.myOmega$lobe <- gsub("RUL", "Right upper lobe (RUL)", plt.myOmega$lobe);
plt.myOmega$lobe <- gsub("RLL", "Right upper lobe (RLL)", plt.myOmega$lobe);
plt.myOmega$CellType <- gsub("CellType","Cell type", plt.myOmega$CellType);

png("~/Downloads/Figure3B.png", res=300, units="in", height=8.27, width=5.84);
ggplot(plt.myOmega, aes(x=CellType, y=value, color=status)) +
  geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA) +
  geom_point(aes(color=status), position=position_jitterdodge(jitter.width=0.2)) + # add jitter
  scale_color_manual(values=c("red","blue")) +
  scale_y_continuous(limits=c(0,1.05)) +
  labs(y="Proportion") +
  myBoxplotTheme
dev.off();
