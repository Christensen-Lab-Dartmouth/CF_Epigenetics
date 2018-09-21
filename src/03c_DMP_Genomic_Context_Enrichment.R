# Test enriched genomic contexts, stratified by direction of change
# Script author: David Chen
# Date: 03/01/2018
# Notes:

rm(list=ls());
library(ggplot2);
library(gridExtra);
library(matrixStats);
library(reshape2);
library(doParallel); registerDoParallel(detectCores() - 1);

FDR_THRESH <- 0.05;
DIR <- "~/Dropbox (Christensen Lab)/Christensen Lab - 2017/Armstrong_CF_project/"; 

labelDifferentialSet <- function(DMPs, pThresh, esThresh) {
  #'@description Label in additional variables based on existing variables
  #'@param DMPs Limma output data.frame merged with Illumina EPIC annotation
  #'@param pThresh,esThresh Threshold for significance & effect size, respectively
  
  DMPs$isConsideredDifferential[DMPs$logFC > esThresh & DMPs$adj.P.Val < pThresh] <- "HYPER";
  DMPs$isConsideredDifferential[DMPs$logFC < -esThresh & DMPs$adj.P.Val < pThresh] <- "HYPO";
  DMPs$isConsideredDifferential[is.na(DMPs$isConsideredDifferential)] <- FALSE; 
  
  ## Add 0/1 indicator columns for transcriptional context:
  DMPs$Promoter <- ifelse(grepl("TSS",DMPs$UCSC_RefGene_Group), 1, 0);
  DMPs$Enhancer <- ifelse(DMPs$Phantom4_Enhancers != "" |  DMPs$Phantom5_Enhancers != "" | DMPs$X450k_Enhancer != "", 1, 0);
  DMPs$Gene.body <- ifelse(grepl("Body",DMPs$UCSC_RefGene_Group), 1, 0);
  DMPs$DHS <- ifelse(DMPs$DNase_Hypersensitivity_NAME != "", 1, 0);
  
  ## Add 0/1 indicator columns for CGI (each category is mutually exclusive of others):
  DMPs$Island <- ifelse(DMPs$Relation_to_Island=="Island", 1, 0);
  DMPs$N_Shelf <- ifelse(DMPs$Relation_to_Island=="N_Shelf", 1, 0);
  DMPs$N_Shore <- ifelse(DMPs$Relation_to_Island=="N_Shore", 1, 0);
  DMPs$OpenSea <- ifelse(DMPs$Relation_to_Island=="OpenSea", 1, 0);
  DMPs$S_Shelf <- ifelse(DMPs$Relation_to_Island=="S_Shelf", 1, 0);
  DMPs$S_Shore <- ifelse(DMPs$Relation_to_Island=="S_Shore", 1, 0);
  
  return(DMPs); 
}

drawGenomicEnrichmentPlots <- function(DMPs, runTests=FALSE, combineSets=FALSE, confLev=0.95) {
  #'@description Explore and/or formally test genomic context enrichment
  #'@param DMPs Limma output data.frame merged with Illumina EPIC annotation
  #'@param runTests Should Fisher's exact test be run?
  #'@param combineSets Should hyper- and hypo-methylated DMP sets be combined into a "differentially methyalted" set?
  #'@param confLev Confidence level. Used only when runTests=TRUE
  
  ## Generate titles for iterative computation:
  titles <- table(DMPs$isConsideredDifferential); 
  titles <- titles[match(c("HYPO","HYPER"), names(titles))]; 
  if(combineSets) titles <- c(titles, differentially=sum(titles));
  titles <- paste(titles, names(titles));
  titles <- paste0(titles, "methylated");
  print(titles);
  pltCols <- length(titles);
  
  ## Iterative computation (Fisher tests & forest plots):
  forestList <- barplotList <- tableList <- list();
  for(directionInCF in titles) {
    ## Subsetting:
    if(grepl("HYPER", directionInCF, ignore.case=FALSE)){
      DMPs.sig <- subset(DMPs, isConsideredDifferential=="HYPER");
      myColor <- "olivedrab3";
    } else if(grepl("HYPO", directionInCF, ignore.case=FALSE)) {
      DMPs.sig <- subset(DMPs, isConsideredDifferential=="HYPO");
      myColor <- "royalblue";
    } else if(grepl("differentially", directionInCF, ignore.case=FALSE)) {
      DMPs.sig <- subset(DMPs, isConsideredDifferential %in% c("HYPER","HYPO")); 
      myColor <- "black";
    }
    
    ## Fisher's tests:
    contexts <- c("Promoter", "Enhancer", "DHS", "Island", "OpenSea");
    summFisher <- data.frame(Context=contexts, propSig=NA, propUniv=NA, OR=NA, CI.lower=NA, CI.upper=NA, P=NA);
    for(gc in contexts) {
      contTab <- rbind(
        Differential = table(DMPs.sig[ , gc]),
        Input = table(DMPs[ , gc])
      );
      contTab <- contTab[ , c(2,1)];
      summFisher$propSig[summFisher$Context==gc] <- contTab[1,1] / sum(contTab[1,1:2]);
      summFisher$propUniv[summFisher$Context==gc] <- contTab[2,1] / sum(contTab[2,1:2]);
      
      if(runTests) {
        fT <- fisher.test(contTab, conf.level=confLev);
        summFisher$OR[summFisher$Context==gc] <- fT$estimate;
        summFisher$CI.lower[summFisher$Context==gc] <- fT$conf.int[[1]]; 
        summFisher$CI.upper[summFisher$Context==gc] <- fT$conf.int[[2]]; 
        summFisher$P[summFisher$Context==gc] <- fT$p.value; 
      }
    }
    summFisher$isSignif <- summFisher$P < 0.05; 
    tableList[[directionInCF]] <- summFisher;
    
    ## Visualization of proportions:
    plt.sideBySideProps <- summFisher[ , c("Context","propSig","propUniv")];
    plt.sideBySideProps <- reshape2::melt(plt.sideBySideProps);
    plt.sideBySideProps$Context <- factor(plt.sideBySideProps$Context, levels=contexts);
    plt.sideBySideProps$variable <- gsub("propSig", "Input", plt.sideBySideProps$variable);
    plt.sideBySideProps$variable <- gsub("propUniv", "Universe", plt.sideBySideProps$variable);
    barplotList[[directionInCF]] <- ggplot(plt.sideBySideProps, aes(Context, value, fill=variable)) +
      geom_bar(stat="identity", position="dodge") +
      scale_y_continuous(limits=c(0,1)) + #empirically determined
      labs(title=directionInCF, y="Proportion") +
      theme_classic() +
      theme(axis.text=element_text(size=20,color="black"), axis.title.x=element_blank(),
            axis.title.y=element_text(size=20,color="black"), title=element_text(size=20,color="black"),
            legend.text=element_text(size=20,color="black",face="bold"), legend.title=element_blank(), legend.position="top");
    
    ## Visualization of Fisher's test results:
    if(runTests) {
      plt.summFisher <- summFisher; #copy
      plt.summFisher$Category <- factor(plt.summFisher$Context, levels=rev(c("Promoter","Enhancer","DHS","Island","OpenSea")) );
      plt.summFisher$P <- signif(plt.summFisher$P, 3);
      
      plt <- ggplot(plt.summFisher, aes(Category, OR, ymin=CI.lower, ymax=CI.upper)) +
        geom_pointrange(size=1, color=myColor) +
        coord_flip() +
        geom_hline(yintercept=1, size=0.2, linetype="dashed", color=myColor) +
        scale_y_continuous(limits=c(-0.5,6)) +
        labs(y="Odds Ratio (OR)", title=paste(directionInCF,"in CF")) +
        theme_classic() +
        theme(axis.line=element_line(color=myColor), axis.ticks=element_line(color=myColor),
              axis.text=element_text(size=20,color=myColor), title=element_text(size=20,face="bold",color=myColor),
              axis.title.x=element_blank(), axis.title.y=element_blank());

      labPos <- nrow(plt.summFisher):1 + 0.15; 
      for(i in 1:nrow(plt.summFisher)) {
        plt <- plt + annotate("text", labPos[i], plt.summFisher$OR[i], label=paste("P =", plt.summFisher$P[i]), size=7, color=myColor);
      }
      
      forestList[[directionInCF]] <- plt;
    }
  }
  
  print(gridExtra::grid.arrange(grobs=barplotList, nrow=pltCols));
  if(runTests) print(gridExtra::grid.arrange(grobs=forestList, ncol=pltCols));
  return(tableList);
}

main <- function() {
  load(paste0(DIR,"113017_DMPs/030118_DMP_mostVariable26K_NoCellType.RData"));
  
  DMPs.original <- labelDifferentialSet(DMPs, pThresh=FDR_THRESH, esThresh=3.50);
  drawGenomicEnrichmentPlots(DMPs.original, runTests=TRUE);
  
  #png("~/Downloads/Figure2B.png", res=300, units="in", height=8.27, width=11.69);
  DMPs.revised <- labelDifferentialSet(DMPs, pThresh=FDR_THRESH, esThresh=2.95);
  drawGenomicEnrichmentPlots(DMPs.revised, runTests=TRUE);
  #dev.off();
}

main(); 

