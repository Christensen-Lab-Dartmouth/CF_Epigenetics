# Volcano & Density Plots for DMPs
# Script author: David Chen
# Date: 03/01/2018
# Revision date: 09/20/2018
# Notes:

rm(list=ls());
library(ggrepel);
library(matrixStats);
library(RnBeads);
library(doParallel); registerDoParallel(detectCores() - 1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsAndPlotThemes.R");

FDR_THRESH <- 0.05;
MVAL_THRESH <- 3.50; #originally set to 3.50; #in M-value
DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/"; 

drawLabeledVolcanoPlot <- function(DMPs, pThresh, mThresh, xScaleUp=1.3, xScaleDown=1.15) {
  ## Direction of change:
  DMPs$dir[DMPs$logFC > mThresh & DMPs$adj.P.Val < pThresh] <- "Positive";
  DMPs$dir[DMPs$logFC < -mThresh & DMPs$adj.P.Val < pThresh] <- "Negative";
  DMPs$dir[abs(DMPs$logFC) < mThresh] <- ""; #special selection (removal) method

  nPos <- sum(DMPs$dir == "Positive", na.rm=TRUE);
  nNeg <- sum(DMPs$dir == "Negative", na.rm=TRUE);
  
  ## Gene names:
  DMP.genes <- as.character(DMPs$UCSC_RefGene_Name); 
  DMP.genes <- strsplit(DMP.genes, split=";");
  DMPs$Gene <- rep(NA, length(DMP.genes) ); 
  for(k in 1:length(DMP.genes)){
    DMP.genes[[k]] <- unique(DMP.genes[[k]]); 
    if(length(DMP.genes[[k]]) > 1) {
      DMPs$Gene[k] <- paste(DMP.genes[[k]], collapse="; "); 
    } else if(length(DMP.genes[[k]]) == 1) {
      DMPs$Gene[k] <- DMP.genes[[k]]; 
    } else if(length(DMP.genes[[k]]) == 0) {
      DMPs$Gene[k] <- NA;
    }
  }
  DMPs$Label <- DMPs$Gene; #initialize
  DMPs$Label[(DMPs$logFC <= xScaleUp*mThresh & DMPs$logFC >= -xScaleDown*mThresh) | DMPs$adj.P.Val >= pThresh] <- NA; #higher threshold for gene labeling
  
  ## Plot:
  plt <- ggplot(DMPs, aes(x=logFC, y=-log10(adj.P.Val), color=dir)) +
    geom_point(aes(size=dir, alpha=dir)) + #override
    scale_x_continuous(expand=c(0,0), limits=c(-7,7), breaks=seq(-5,5,by=1)) +
    scale_color_manual(values = c("lightgray","royalblue","olivedrab3")) +
    scale_size_manual(values=c(1, 1.1, 1.1)) +
    scale_alpha_manual(values=c(0.9, 1, 1)) +
    geom_text_repel(aes(label=Label), color="black", size=4) +
    geom_hline(yintercept=-log10(pThresh), linetype="dashed") +
    geom_vline(xintercept=c(-mThresh,mThresh), linetype="dashed") +
    labs(x="log2FC(M-value)", y="-log10(FDR)") +
    annotate("text",5.2,3.5, label=paste(nPos,"HYPERmethylated"),size=6.5,color="olivedrab3") +
    annotate("text",-5.2,3.5,label=paste(nNeg,"HYPOmethylated"),size=6.5,color="royalblue") +
    myVolcanoTheme
  print(plt);
}

main <- function() {
  ## Load differential analysis results & summary stats:
  load(paste0(DATA_PATH,"030118_DMP_mostVariable26K_NoCellType.RData"));
  load(paste0(DATA_PATH,"030218_log2BetaFCs_for_visualization_purposes.RData"));
  DMPs <- merge(DMPs, log2BetaFCs, by="Name");
  
  ## Export in response to reviewer:
  dmpSub <- subset(DMPs, isSignifAt0.05); 
  dmpSub$isConsideredDifferential <- abs(dmpSub$logFC) >= MVAL_THRESH;
  write.csv(dmpSub, file="~/Downloads/full_CpGs_at_FDR0.05.csv", row.names=FALSE, quote=FALSE);
  
  ## Volcano plot with x-axis being log2FC(M):
  png("~/Downloads/Figure2A.png", res=300, units="in", height=8.27, width=11.69);
  drawLabeledVolcanoPlot(DMPs, pThresh=FDR_THRESH, mThresh=MVAL_THRESH);
  dev.off();
}

main(); 
