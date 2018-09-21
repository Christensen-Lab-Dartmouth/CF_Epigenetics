# Comparison of Cell-Type Proportions by Statistical Tests
# Script author: David Chen
# Date: 09/24/17; 04/13/18
# Notes:

rm(list=ls());
library(ggplot2);
library(lme4); library(lmerTest);
library(matrixStats);
library(reshape2);
library(doParallel); registerDoParallel(detectCores() - 1);
setwd(dirname(rstudioapi::getActiveDocumentContext()$path));
source("HelperFunctionsAndPlotThemes.R");
DATA_PATH <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/";

#-----------------------------------------------Statistical tests-----------------------------------------------
load(paste0(DATA_PATH, "112017_MasterCovariates_BAL.RData"));
masterCovar$Status <- factor(masterCovar$Status, levels=c("H","CF")); #reorder factor level

fit.CT1 <- lmer(
  CellType1 ~ (1 | Subject) + Status + Age + Sex,
  masterCovar
);
summary(fit.CT1)
## Random effects:
##   Groups   Name        Variance Std.Dev.
## Subject  (Intercept) 0.01797  0.1341  
## Residual             0.01227  0.1108  
## Number of obs: 16, groups:  Subject, 8
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)   
## (Intercept)   1.55763    0.33785  4.00000   4.610  0.00995 **
##  StatusCF    -0.54285    0.12901  4.00000  -4.208  0.01361 * 
##  Age         -0.02790    0.01295  4.00000  -2.154  0.09751 . 
##  SexM         0.27962    0.12307  4.00000   2.272  0.08553 . 

## Statistical test 2:
fit.CT2 <- lmer(
  CellType2 ~ (1 | Subject) + Status + Age + Sex,
  masterCovar
);
summary(fit.CT2)
## Random effects:
##   Groups   Name        Variance Std.Dev.
## Subject  (Intercept) 0.02071  0.1439  
## Residual             0.01265  0.1125  
## Number of obs: 16, groups:  Subject, 8
## 
## Fixed effects:
##             Estimate Std. Error       df t value Pr(>|t|)  
## (Intercept) -0.55708    0.35781  3.99999  -1.557   0.1945  
## StatusCF     0.52739    0.13663  3.99999   3.860   0.0181 *
## Age          0.02730    0.01372  3.99999   1.991   0.1174  
## SexM        -0.27046    0.13034  3.99999  -2.075   0.1066  

#-----------------------------------------------Data visualization-----------------------------------------------
plt.myOmega <- masterCovar[ , c("Sample_Name","CellType1","CellType2","Status","Lobe")];
plt.myOmega <- melt(plt.myOmega);

plt.myOmega$Status <- gsub("CF", "Cystic fibrosis (CF)", plt.myOmega$Status);
plt.myOmega$Status <- gsub("H", "Healthy (H)", plt.myOmega$Status);

plt.myOmega$Lobe <- gsub("RUL", "Right upper lobe (RUL)", plt.myOmega$Lobe);
plt.myOmega$Lobe <- gsub("RLL", "Right upper lobe (RLL)", plt.myOmega$Lobe);

png("~/Downloads/Figure2A_with_p_values.png", res=300, units="in", height=8.27, width=5.84);
ggplot(plt.myOmega, aes(x=variable, y=value, color=Status)) +
  geom_boxplot(outlier.colour=NA, outlier.fill=NA, outlier.shape=NA) +
  geom_point(aes(color=Status), position=position_jitterdodge(jitter.width=0.2)) + # add jitter
  scale_color_manual(values=c("red","blue")) +
  scale_y_continuous(limits=c(0,1.06)) +
  labs(y="Proportion") +
  myBoxplotTheme +
  
  geom_segment(aes(x=0.81, y=1.0, xend=1.19, yend=1.0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=0.81, y=0.95, xend=0.81, yend=1.0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1.19, y=0.99, xend=1.19, yend=1.0), size=0.3, inherit.aes=F) +
  annotate("text", x=1, y=1.03, label="P = 0.014", size=6) +
  
  geom_segment(aes(x=1.81, y=1.0, xend=2.19, yend=1.0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=1.81, y=0.97, xend=1.81, yend=1.0), size=0.3, inherit.aes=F) +
  geom_segment(aes(x=2.19, y=0.25, xend=2.19, yend=1.0), size=0.3, inherit.aes=F) +
  annotate("text", x=2, y=1.03, label="P = 0.018", size=6)
dev.off();
