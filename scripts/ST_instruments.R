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
#' Sup Tab 2: all MR-IVW results
#' 
#' Sup Tab 3: Instruments for the MR ratio estimate + ratio
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

IVData[,BCSurvival_EAF := out1$EAF]
IVData[,BCSurvival_nSamples := out1$nSamples]
IVData[,BCSurvival_logHR := out1$beta]
IVData[,BCSurvival_se := out1$se]
IVData[,BCSurvival_pval := out1$pval]

plot(IVData$PCSK9_EAF,IVData$BCSurvival_EAF)
abline(0,1)

#' add breast cancer summary statistics from FinnGen and UKB meta-analysis
out2 = copy(outcomeData)
out2 = out2[phenotype == "BCP_meta",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,BreastCancer_EAF := out2$EAF]
IVData[,BreastCancer_nSamples := out2$nSamples]
IVData[,BreastCancer_logOR := out2$beta]
IVData[,BreastCancer_se := out2$se]
IVData[,BreastCancer_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$BreastCancer_EAF)
abline(0,1)

#' add parental age at death summary statistics from UKB 
out2 = copy(outcomeData)
out2 = out2[phenotype == "Parents' age at death",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,AgeAtDeath_EAF := out2$EAF]
IVData[,AgeAtDeath_nSamples := out2$nSamples]
IVData[,AgeAtDeath_beta := out2$beta]
IVData[,AgeAtDeath_se := out2$se]
IVData[,AgeAtDeath_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$AgeAtDeath_EAF)
abline(0,1)

stab1 = copy(IVData)

#' # Sup Tab 2 ####
#' ***
load("../results/05_MR_IVW.RData")
MR_IVW = MR_IVW[exposure != "PCSK9_females_free"]
MR_IVW = MR_IVW[outcome %in% c("BCP_meta","BCS_meta", "Parents' age at death")]
MR_IVW[outcome == "BCP_meta", outcome := "BreastCancer [FinnGen+UKB]"]
MR_IVW[outcome == "BCS_meta", outcome := "BCSurvival [BCAC]"]
MR_IVW[outcome == "Parents' age at death", outcome := "AgeAtDeath [UKB]"]

stab2 = copy(MR_IVW)

#' # Sup Tab 3 ####
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
MR_ratio = MR_ratio[pval<0.05,]
MR_ratio = MR_ratio[outcome_phenotype %in% c("BCP_meta","BCS_meta", "Parents' age at death","EUR only")]
MR_ratio = MR_ratio[phenotype %in% c("Brain_Cerebellum","Esophagus_Muscularis","Spleen","PCSK9_females","PCSK9_free" )]

stab3 = MR_ratio[,c(1:12,13,14,15,17:20,22,24,25)]
names(stab3)[8:12] = paste0("PCSK9_",names(stab3)[8:12])
names(stab3)[19] = "beta_ratio"
names(stab3)[20] = "se_ratio"
names(stab3)[21] = "pval_ratio"
stab3[outcome_phenotype == "BCP_meta", outcome_phenotype := "BreastCancer [FinnGen+UKB]"]
stab3[outcome_phenotype == "BCS_meta", outcome_phenotype := "BCSurvival [BCAC]"]
stab3[outcome_phenotype == "Parents' age at death", outcome_phenotype := "AgeAtDeath [UKB]"]
stab3[outcome_phenotype == "EUR only", outcome_phenotype := "BCSurvival [Mei et al.]"]
names(stab3)[7] = "exposure"
names(stab3)[13] = "outcome"

#' # Save ####
#' ***
WriteXLS(x = c("stab1","stab2","stab3"), 
         ExcelFileName=paste0("../results/SupTables.xlsx"), 
         SheetNames=c("TableS1","TableS2","TableS3"), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(stab1,stab2,stab3,file = "../results/SupTables.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
