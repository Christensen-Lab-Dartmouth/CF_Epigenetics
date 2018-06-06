###############################################################################################################
# Volcano & Density Plots for DMPs
# Script author: David Chen
# Date: 03/01/2018
# Notes:
###############################################################################################################

rm(list=ls())

library(ggrepel)
library(matrixStats)
library(RnBeads)
library(doParallel); registerDoParallel(detectCores() - 1)

## Load differential analysis results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/030118_DMP_mostVariable26K_NoCellType.RData"); 

## Load beta-value summary statistics:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/113017_DMPs/030218_log2BetaFCs_for_visualization_purposes.RData")

## Merge:
DMPs <- merge(DMPs, log2BetaFCs, by="Name");

## Define a flexible set of thresholds:
pThresh <- 0.05;
mThresh <- 3.50; #M-value

## Code in unique genes as labels:
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

## Volcano plot:
DMPs$dir[DMPs$logFC > mThresh & DMPs$adj.P.Val < pThresh] <- "Positive";
DMPs$dir[DMPs$logFC < -mThresh & DMPs$adj.P.Val < pThresh] <- "Negative";
DMPs$dir[abs(DMPs$logFC) < mThresh] <- ""; #special selection (removal) method
sum(DMPs$dir != "", na.rm=T)

## Count information for adding to volcano plot:
sum(DMPs$dir == "Positive", na.rm=T)
sum(DMPs$dir == "Negative", na.rm=T)

## Add labels:
DMPs$Label <- DMPs$Gene; #initialize
DMPs$Label[(DMPs$logFC < 1.05*mThresh & DMPs$logFC > -1.15*mThresh) | DMPs$adj.P.Val >= pThresh] <- NA; #higher threshold for gene labeling
sum(! is.na(DMPs$Label))

## Volcano plot: M-value log2FC
png("~/Downloads/Figure2A.png", res=300, units="in", height=8.27, width=11.69);
ggplot(DMPs, aes(x=logFC, y=-log10(adj.P.Val), color=dir)) +
  geom_point(aes(size=dir, alpha=dir)) + #override
  scale_x_continuous(expand=c(0,0), limits=c(-7,7)) +
  scale_color_manual(values = c("lightgray","royalblue","olivedrab3")) +
  scale_size_manual(values=c(1, 1.05, 1.05)) +
  scale_alpha_manual(values=c(0.9, 1, 1)) +
  geom_text_repel(aes(label=Label), color="black", size=4) +
  geom_hline(yintercept=-log10(pThresh), linetype="dashed") +
  geom_vline(xintercept=c(-mThresh,mThresh), linetype="dashed") +
  theme_classic() +
  labs(x="log2 M-value fold change", y="-log10(FDR)") +
  theme(axis.text.x=element_text(color="black",size=16),axis.title.x=element_text(size=21, color="black"), 
        axis.text.y=element_text(color="black",size=16),axis.title.y=element_text(size=21, color="black"),
        legend.title=element_blank(), legend.text=element_blank(), legend.position="none") +
  annotate("text",5.2,3.5, label="51 HYPERmethylated",size=6.5,color="olivedrab3") +
  annotate("text",-5.2,3.5,label="58 HYPOmethylated",size=6.5,color="royalblue")
dev.off();

