####################################################################################
# MethylationEPIC RefFree2.0 Result Visualization
# Script author: David Chen
# Date: 09/24/2017
# Notes:  
####################################################################################

rm(list=ls())

library(gdata)
library(ggplot2)
library(matrixStats)
library(pheatmap)
library(RColorBrewer); 
library(RefFreeEWAS)
library(reshape2)
library(doParallel); registerDoParallel(detectCores() - 1)

#----------------------------------------------Decide optimal K---------------------------------------------
RefFree_MO_Boots <- RefFreeCellMixArrayDevianceBoots(RefFree_Array_MO2, Y_shortened);
# save(RefFree_MO_Boots, file="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/091917_MO_EPIC/092517_local_solution_for_bootstraps.RData", compress=TRUE)

## Load cluster-computed results:
load("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/091917_MO_EPIC/092417_RefFree_on_10K_Macrophage_EPIC.RData")

mean_Ks <- apply(RefFree_MO_Boots[-1, ], 2, mean, trim=0.25);
optimalK <- names(mean_Ks)[which.min(mean_Ks)];
optimalK

## Plot distribution of matrices Mu:
par(mfrow=c(3,3), oma=c(2,0,2,0))
for(x in 1:length(RefFree_Array_MO)){
  hist(RefFree_Array_MO[[x]]$Mu, main=paste0("K=",names(RefFree_Array_MO)[x]), xlab="", col="royalblue", border=F);
}
mtext("Distribution of matrix Mu for 9 presumptive K's", outer=TRUE);

## Check constrain the range to unit interval and check distribution of matrices Omega:
for(x in 1:length(RefFree_Array_MO)){
  RefFree_Array_MO[[x]]$Omega[RefFree_Array_MO[[x]]$Omega < 0] <- 0.0;
  RefFree_Array_MO[[x]]$Omega[RefFree_Array_MO[[x]]$Omega > 1] <- 1.0;
  hist(RefFree_Array_MO[[x]]$Omega, main=paste0("K=",names(RefFree_Array_MO)[x]), xlab="", col="orange", border=F);
}
mtext("Distribution of matrix Omega for 9 presumptive K's", outer=TRUE);

#---------------------------------------------Heat map visualization---------------------------------------------
## Load data & switch out barcode with sample ID:
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

## (Added) Load & merge in RPMM class membership:
dat_sup <- read.csv("~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/092517_RPMM/092517_RPMM_class_membership.csv")
colnames(dat_sup) <- gsub("class", "RPMM", colnames(dat_sup)); 

## Heat map annotations:
heat_annot <- data.frame(
  Sample_Name = target_mo$Sample_Name,
  Status = target_mo$status,
  Lobe = target_mo$lobe
);
heat_annot <- merge(heat_annot, dat_sup[ , c("Sample_Name","RPMM")], by="Sample_Name");
rownames(heat_annot) <- heat_annot$Sample_Name;
heat_annot$Sample_Name <- NULL;

ann_colors <- list(
  Lobe = c(RUL="black", RLL="lightgray"),
  Status = c(CF="red",H="blue"),
  RPMM = c(rR="black",rL="lightgray")
);

## Visualize Omega_optimalK:
myOmega <- RefFree_Array_MO[[optimalK]]$Omega;
if( all(match(target_mo$Complete_Barcode, rownames(myOmega) ) == 1:16) ) {
  print("Already in order, no further action necessary!");
  rownames(myOmega) <- target_mo$Sample_Name;
} else {
  myOmega <- myOmega[ , match(target_mo$Complete_Barcode, rownames(myOmega)) ];
  rownames(myOmega) <- target_mo$Sample_Name;
  print("Match & name swithing complete!")
}

fixInNamespace("draw_colnames","pheatmap") #vjust = 1, hjust = 0.5, rot = 0
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

rbind(paste0("cell type",colnames(myOmega)), colVars(myOmega)) #check variability

# write.csv(myOmega, file="~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/091917_MO_EPIC/092517_RefFree2_OmegaK=2.csv",row.names=T, quote=F)

#---------------------------------------------Box plot visualization---------------------------------------------
## Data frame assembly:
plt.myOmega <- myOmega; #make a copy
colnames(plt.myOmega) <- paste0("CellType",colnames(plt.myOmega)); 
plt.myOmega <- melt(plt.myOmega, varnames=c("Sample_Name","CellType"));
plt.myOmega <- merge(plt.myOmega, target_mo[ , c("Sample_Name","status","lobe")], all.x=TRUE); #merge in values

plt.myOmega$status <- gsub("CF", "Cystic fibrosis (CF)", plt.myOmega$status);
plt.myOmega$status <- gsub("H", "Healthy (H)", plt.myOmega$status);

plt.myOmega$lobe <- gsub("RUL", "Right upper lobe (RUL)", plt.myOmega$lobe);
plt.myOmega$lobe <- gsub("RLL", "Right upper lobe (RLL)", plt.myOmega$lobe);

plt.myOmega$CellType <- gsub("CellType","Cell type", plt.myOmega$CellType);

## Data visualization
png("~/Downloads/Figure3B.png", res=300, units="in", height=8.27, width=5.84);
ggplot(plt.myOmega, aes(x=CellType, y=value, color=status)) +
  geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA) +
  geom_point(aes(color=status), position=position_jitterdodge(jitter.width=0.2)) + # add jitter
  scale_color_manual(values=c("red","blue")) +
  scale_y_continuous(limits=c(0,1.05)) +
  labs(y="Proportion") +
  theme_classic() +
  theme(
    axis.text.x=element_text(size=20,color="black"), axis.title.x=element_blank(),
    axis.text.y=element_text(size=20,color="black"), axis.title.y=element_text(size=20,color="black"),
    legend.text=element_text(size=20,color="black",face="bold"), legend.title=element_blank(), legend.position="top"
  )
dev.off();
