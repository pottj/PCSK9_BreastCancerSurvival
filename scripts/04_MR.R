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
load("../results/03b_Exposure_for_MR_filtered.RData")
load("../results/03b_Outcome_for_MR_filtered.RData")

#' # Filter for relevant phenotype ####
#' ***
#' I also want a PCSK9_HMGCR combination for LDL-C.
#' 
dummy1 = copy(ExposureData)
dummy1 = dummy1[phenotype == "LDL-C levels",]
dummy1 = dummy1[grepl("PCSK9",instrumentLocus) | grepl("HMGCR",instrumentLocus),]
dummy1 = dummy1[MRapproach == "QTL"]
dummy1[,instrumentLocus := "PCSK9 and HMGCR"]
ExposureData = rbind(ExposureData,dummy1)

ExposureData[,trait := paste(phenotype,setting,instrumentLocus,MRapproach,sep=" - ")]
myTraits = unique(ExposureData$trait)
myTraits2 = gsub("LDL-C levels","LDLC",myTraits)
myTraits2 = gsub("PCSK9 GE levels","PCSK9",myTraits2)
myTraits2 = gsub("PCSK9 PE levels","PCSK9",myTraits2)
myTraits2 = gsub(" - ","_",myTraits2)
myTraits2 = gsub(" ","",myTraits2)
myTraits2 = gsub("-","",myTraits2)
myTraits2

OutcomeData[,trait := paste(phenotype,setting,source,geneticModel,sep = " - ")]
myYTraits = unique(OutcomeData$trait)
myYTraits2 = gsub(" et al.","",myYTraits)
myYTraits2 = gsub(" - ","_",myYTraits2)
myYTraits2 = gsub(" & ","",myYTraits2)
myYTraits2 = gsub("-","",myYTraits2)
myYTraits2

dumTab1 = foreach(i = 1:length(myTraits))%do%{
  #i=1
  exposure2 = copy(ExposureData)
  exposure2 = exposure2[trait == myTraits[i],]
  
  outcome2 = copy(OutcomeData)
  outcome2 = outcome2[rsID %in% exposure2$rsID]
  myOutcomes = unique(outcome2$trait)
  myOutcomes2 = gsub(" et al.","",myOutcomes)
  myOutcomes2 = gsub(" - ","_",myOutcomes2)
  myOutcomes2 = gsub(" & ","",myOutcomes2)
  myOutcomes2 = gsub("-","",myOutcomes2)
  myOutcomes2
  
  dumTab2 = foreach(j = 1:length(myOutcomes))%do%{
    #j=1
    outcome = copy(outcome2)
    outcome = outcome[trait == myOutcomes[j],]
    
    setorder(exposure2,chr,pos_b37)
    setorder(outcome,chr,pos_b37)
    exposure = copy(exposure2)
    exposure = exposure[rsID %in% outcome$rsID]
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
                        "___",myOutcomes2[j],".png")
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

