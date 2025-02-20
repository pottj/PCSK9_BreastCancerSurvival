#' ---
#' title: "Get supplemental table"
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
#' Sup Tab 1: Instruments for the MR-IVW
#' 
#' Sup Tab 2: Instruments for the MR ratio estimate
#' 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Sup Tab 1 ####
#' ***
load("../results/04_IVs_pruned.RData")
load("../results/03_outcome_harmonized.RData")

IVData = IVData[phenotype != "PCSK9_females_free",]
IVData[,flag := NULL]
names(IVData)[8:12] = paste0("PCSK9_",names(IVData)[8:12])

#' add breast cancer survival summary statistics from BCAC
out1 = copy(outcomeData)
out1 = out1[phenotype == "BCS_meta",]
matched1 = match(IVData$rsID,out1$rsID)
out1 = out1[matched1,]

IVData[,BCS_EAF := out1$EAF]
IVData[,BCS_nSamples := out1$nSamples]
IVData[,BCS_beta := out1$beta]
IVData[,BCS_se := out1$se]
IVData[,BCS_pval := out1$pval]

plot(IVData$PCSK9_EAF,IVData$BCS_EAF)
abline(0,1)

#' add breast cancer summary statistics from FinnGen and UKB meta-analysis
out2 = copy(outcomeData)
out2 = out2[phenotype == "BCP_meta",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,BC_EAF := out2$EAF]
IVData[,BC_nSamples := out2$nSamples]
IVData[,BC_beta := out2$beta]
IVData[,BC_se := out2$se]
IVData[,BC_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$BC_EAF)
abline(0,1)

stab1 = copy(IVData)

#' # Sup Tab 2 ####
#' ***
load("../results/06_MR_ratio_rs562556.RData")
MR_ratio_BCAC = copy(MR_ratio)

load("../results/07_MR_ratio_rs562556_Mei.RData")
MR_ratio_Mei = copy(MR_ratio)

MR_ratio_BCAC[,outcome_source := "Morra et al."]
MR_ratio_Mei[,outcome_source := "Mei et al."]

names(MR_ratio_BCAC)
names(MR_ratio_Mei)
names(MR_ratio_Mei)[13:15] = c("phenotype","beta","se" )

MR_ratio = rbind(MR_ratio_BCAC,MR_ratio_Mei,use.names=T,fill=T)
names(MR_ratio)[13:19] = paste0("outcome_",names(MR_ratio)[13:19])

MR_ratio = MR_ratio[outcome_phenotype %in% c("BCS_meta","EUR only")]

stab2 = MR_ratio[c(1:12,15,20:31,34),c(1:12,26,14,15,17:20,22,24,25)]
names(stab2)[8:12] = paste0("PCSK9_",names(stab2)[8:12])
names(stab2)[19] = "beta_ratio"
names(stab2)[20] = "se_ratio"
names(stab2)[21] = "pval_ratio"

#' # Save ####
#' ***
WriteXLS(x = c("stab1","stab2"), 
         ExcelFileName=paste0("../results/SupTables.xlsx"), 
         SheetNames=c("TableS1","TableS2"), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
