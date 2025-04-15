#' ---
#' title: "MR-IVW: PCSK9 on BC outcome"
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
#' In this script, I test all exposure - outcome combinations using MR-IVW. 
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
load("../results/04_IVs_pruned.RData")
load("../results/03_outcome_harmonized.RData")

IVData = IVData[pval<5e-8,]

#' # Filter for relevant phenotype ####
#' ***
myTraits = unique(IVData$phenotype)
myOutcomes = unique(outcomeData$phenotype)

dumTab1 = foreach(i = 1:length(myTraits))%do%{
  #i=1
  exposure = copy(IVData)
  exposure = exposure[phenotype == myTraits[i],]
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    #j=1
    outcome = copy(outcomeData)
    outcome = outcome[phenotype == myOutcomes[j],]
    outcome = outcome[rsID %in% exposure$rsID,]
    
    if(dim(outcome)[1] == 0){
      res = data.table(exposure = myTraits[i],
                       outcome = myOutcomes[j],
                       nSNPs = 0)
      res
    }else{
      exposure2 = exposure[rsID %in% outcome$rsID]
      nSNPs = dim(exposure2)[1]
      mrob = mr_input(bx = exposure2$beta, bxse = exposure2$se,
                      by = outcome$beta, byse = outcome$se,
                      snps = paste0(exposure2$rsID,"\n"), 
                      exposure = myTraits[i], outcome = myOutcomes[j])
      
      mod1= mr_ivw(mrob)
      
      mr_plot(mrob, interactive = F,labels = T)
      ggsave(paste0("../results/05_ScatterPlot_IVW/",myTraits[i],"_",myOutcomes[j],".png"),
             height = 7, width = 7)
      
      res = data.table(exposure = mod1@Exposure,
                       outcome = mod1@Outcome,
                       nSNPs = nSNPs,
                       beta_IVW = mod1@Estimate,
                       se_IVW = mod1@StdError,
                       pval_IVW = mod1@Pvalue,
                       FStat = mod1@Fstat, 
                       Q = mod1@Heter.Stat[1],
                       pval_Q = mod1@Heter.Stat[2])
      res
      
    }
  }
  res1 = rbindlist(dumTab2,fill=T)
  res1
}
MR_IVW = rbindlist(dumTab1)
MR_IVW[pval_IVW <0.05]

#' # Save data ####
#' ***
save(MR_IVW, file = "../results/05_MR_IVW.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

