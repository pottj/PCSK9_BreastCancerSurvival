#' ---
#' title: "MR: PCSK9 on outcomes"
#' subtitle: "MR of PCSK9 on Breast Cancer Survival"
#' author: "Janne Pott"
#' date: "Last compiled on `r format(Sys.time(), '%d %B, %Y')`"
#' output:
#'   html_document:
#'     toc: true
#'     number_sections: true
#'     toc_float: true
#'     code_folding: show
#' ---
#'
#' # Introduction ####
#' ***
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/03_Exposure_for_MR_pruned.RData")
load("../results/02_Outcome_for_MR.RData")
OutcomeData = OutcomeData[rsID %in% ExposureData$rsID,]
OutcomeData[,chr := as.numeric(chr)]

#' # Filter for relevant phenotype ####
#' ***
ExposureData[,trait := paste(phenotype,setting,sep = " - ")]
myTraits = unique(ExposureData$trait)
myTraits2 = gsub("LDL-C levels - ","LDLC_",myTraits)
myTraits2 = gsub("PCSK9 gene expression - ","GE_",myTraits2)
myTraits2 = gsub("PCSK9 protein levels - ","PE_",myTraits2)
myTraits2 = gsub("[()]","",myTraits2)
myTraits2 = gsub(" -","",myTraits2)
myTraits2 = gsub(" ","_",myTraits2)

OutcomeData[,trait := paste(phenotype,setting,sep = " - ")]

dumTab1 = foreach(i = 1:length(myTraits))%do%{
  #i=17
  exposure = copy(ExposureData)
  exposure = exposure[trait == myTraits[i],]
  
  outcome2 = copy(OutcomeData)
  outcome2 = outcome2[rsID %in% exposure$rsID]
  myOutcomes = unique(outcome2$trait)
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    #j=1
    outcome = copy(outcome2)
    outcome = outcome[trait == myOutcomes[j],]
    
    setorder(exposure,chr,pos_b37)
    setorder(outcome,chr,pos_b37)
    stopifnot(exposure$rsID == outcome$rsID)
    
    if(dim(outcome)[1] == 1){
      test = MRfunction_jp(betaX = exposure$beta, seX = exposure$se, 
                           betaY = outcome$beta, seY = outcome$se)
      
      res = data.table(exposure = myTraits[i],
                       outcome = myOutcomes[j],
                       nSNPs = 1,
                       method = "MR-ratio",
                       beta = test$beta_IV,
                       se = test$se_IV2,
                       pval = test$p_IV2,
                       FStat = test$Fstat)
      res
    }else{
      nSNPs = dim(exposure)[1]
      mrob = mr_input(bx = exposure$beta, bxse = exposure$se,
                      by = outcome$beta, byse = outcome$se,
                      snps = paste0(exposure$rsID,"\n"), 
                      exposure = myTraits[i], outcome = myOutcomes[j])
      
      mod1= mr_ivw(mrob)
      
      mr_plot(mrob, interactive = F,labels = T)
      filename = paste0("../results/04_ScatterPlot_IVW/",myTraits2[i],
                        "___",myOutcomes[j],".png") 
      filename = gsub(" - ","_",filename)
      filename = gsub(" + ","_",filename)
      filename = gsub(" ","_",filename)
      
      ggsave(filename,
             height = 7, width = 7)
      
      res = data.table(exposure = mod1@Exposure,
                       outcome = mod1@Outcome,
                       nSNPs = nSNPs,
                       method = "MR-IVW",
                       beta = mod1@Estimate,
                       se = mod1@StdError,
                       pval = mod1@Pvalue,
                       FStat = mod1@Fstat, 
                       Q = mod1@Heter.Stat[1],
                       pval_Q = mod1@Heter.Stat[2])
      res
      
    }
  }
  res1 = rbindlist(dumTab2,fill=T)
  res1
}
MRTab = rbindlist(dumTab1,fill = T)
MRTab[pval <0.05, .N,by=c("outcome","method")]

#' # Save data ####
#' ***
save(MRTab, file = "../results/04_MR.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

