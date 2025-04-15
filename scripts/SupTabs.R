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
#' Sup Tab 2: all MR-IVW results for PCSK9
#' 
#' Sup Tab 3: Instruments for the MR ratio estimate + ratio for PCSK9
#' 
#' Sup Tab 4: Instruments for the MR-IVW of LDLC
#' 
#' Sup Tab 5: all MR-IVW results for LDLC 
#' 
#' Sup Tab 6: Instruments for the MVMR-IVW of PCSK9 and LDLC
#' 
#' Sup Tab 7: all MVMR-IVW results for PCSK9 and LDLC
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
out1 = out1[phenotype == "BCS",]
matched1 = match(IVData$rsID,out1$rsID)
out1 = out1[matched1,]

IVData[,BCSurvival_EAF := out1$EAF]
IVData[,BCSurvival_nSamples := out1$nSamples]
IVData[,BCSurvival_nCases := out1$nCases]
IVData[,BCSurvival_logHR := out1$beta]
IVData[,BCSurvival_se := out1$se]
IVData[,BCSurvival_pval := out1$pval]

plot(IVData$PCSK9_EAF,IVData$BCSurvival_EAF)
abline(0,1)

#' add breast cancer summary statistics from FinnGen and UKB meta-analysis
out2 = copy(outcomeData)
out2 = out2[phenotype == "BC",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,BreastCancer_EAF := out2$EAF]
IVData[,BreastCancer_nSamples := out2$nSamples]
IVData[,BreastCancer_nCases := out2$nCases]
IVData[,BreastCancer_logOR := out2$beta]
IVData[,BreastCancer_se := out2$se]
IVData[,BreastCancer_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$BreastCancer_EAF)
abline(0,1)

#' add coronary atherosclerosis statistics from FinnGen and UKB meta-analysis
out2 = copy(outcomeData)
out2 = out2[phenotype == "CAD",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,CAD_EAF := out2$EAF]
IVData[,CAD_nSamples := out2$nSamples]
IVData[,CAD_nCases := out2$nCases]
IVData[,CAD_logOR := out2$beta]
IVData[,CAD_se := out2$se]
IVData[,CAD_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$CAD_EAF)
abline(0,1)

#' add parental age at death summary statistics from UKB 
out2 = copy(outcomeData)
out2 = out2[phenotype == "PAAD",]
matched2 = match(IVData$rsID,out2$rsID)
out2 = out2[matched2,]

IVData[,ParentsAgeAtDeath_EAF := out2$EAF]
IVData[,ParentsAgeAtDeath_nSamples := out2$nSamples]
IVData[,ParentsAgeAtDeath_beta := out2$beta]
IVData[,ParentsAgeAtDeath_se := out2$se]
IVData[,ParentsAgeAtDeath_pval := out2$pval]

plot(IVData$PCSK9_EAF,IVData$AgeAtDeath_EAF)
abline(0,1)

stab1 = copy(IVData)
setnames(stab1,"phenotype","exposure")

#' # Sup Tab 2 ####
#' ***
load("../results/05_MR_IVW.RData")
MR_IVW[outcome == "BC", outcome := "BreastCancer [FinnGen+UKB]"]
MR_IVW[outcome == "BCS", outcome := "BCSurvival [BCAC]"]
MR_IVW[outcome == "CAD", outcome := "CAD [FinnGen+UKB]"]
MR_IVW[outcome == "PAAD", outcome := "ParentsAgeAtDeath [UKB]"]

stab2 = copy(MR_IVW)
stab2 = stab2[exposure != "Testis",]

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
names(MR_ratio_Mei)[13:16] = paste0("outcome_",names(MR_ratio_Mei)[13:16])
MR_ratio_Mei = MR_ratio_Mei[outcome_phenotype == "EUR only"]
MR_ratio_Mei[,outcome_phenotype := "BCS_Mei"]
MR_ratio_BCAC[outcome_phenotype=="BCS",outcome_phenotype := "BCS_Morra"]

MR_ratio = rbind(MR_ratio_BCAC,MR_ratio_Mei,use.names=T,fill=T)
names(MR_ratio) 

stab3 = MR_ratio[,c(1:20,22,24,25)]
names(stab3)[7:12] = paste0("exposure_",names(stab3)[7:12])
names(stab3)[20] = "beta_ratio"
names(stab3)[21] = "se_ratio"
names(stab3)[22] = "pval_ratio"

stab3[outcome_phenotype == "BC", outcome_phenotype := "BreastCancer [FinnGen+UKB]"]
stab3[outcome_phenotype == "CAD", outcome_phenotype := "CAD [FinnGen+UKB]"]
stab3[outcome_phenotype == "BCS_Morra", outcome_phenotype := "BCSurvival [BCAC]"]
stab3[outcome_phenotype == "PAAD", outcome_phenotype := "ParentsAgeAtDeath [UKB]"]
stab3[outcome_phenotype == "BCS_Mei", outcome_phenotype := "BCSurvival [Mei et al.]"]

names(stab3)[7] = "exposure"
names(stab3)[13] = "outcome"

stab3 = stab3[exposure != "Testis",]

#' # Sup Tab 4 ####
#' ***
#' Instruments for MR-IVW for LDLC
#' 
load("../temp/09_MR_IVW_input.RData")

mySNPs_PCSK9 = c("rs2495491","rs11591147","rs693668","rs11583680")
mySNPs_HMGCR = c("rs12916","rs1525764")
mySNPs_PCSK9_HMGCR = c(mySNPs_PCSK9,mySNPs_HMGCR)

stab4 = copy(LDLC)
names(stab4)[8:12] = paste0("LDLC_",names(stab4)[8:12])

#' add breast cancer survival summary statistics from BCAC
stopifnot(stab4$rsID == outcomeData[phenotype=="BCS",rsID])
stab4[,BCSurvival_EAF := outcomeData[phenotype=="BCS",EAF]]
stab4[,BCSurvival_nSamples := outcomeData[phenotype=="BCS",nSamples]]
stab4[,BCSurvival_nCases := outcomeData[phenotype=="BCS",nCases]]
stab4[,BCSurvival_logHR := outcomeData[phenotype=="BCS",beta]]
stab4[,BCSurvival_se := outcomeData[phenotype=="BCS",se]]
stab4[,BCSurvival_pval := outcomeData[phenotype=="BCS",pval]]

plot(stab4$LDLC_EAF,stab4$BCSurvival_EAF)
abline(0,1)

#' add breast cancer summary statistics from FinnGen and UKB meta-analysis
stopifnot(stab4$rsID == outcomeData[phenotype=="BC",rsID])
stab4[,BreastCancer_EAF := outcomeData[phenotype=="BC",EAF]]
stab4[,BreastCancer_nSamples := outcomeData[phenotype=="BC",nSamples]]
stab4[,BreastCancer_nCases := outcomeData[phenotype=="BC",nCases]]
stab4[,BreastCancer_logOR := outcomeData[phenotype=="BC",beta]]
stab4[,BreastCancer_se := outcomeData[phenotype=="BC",se]]
stab4[,BreastCancer_pval := outcomeData[phenotype=="BC",pval]]

plot(stab4$LDLC_EAF,stab4$BreastCancer_EAF)
abline(0,1)

#' add CAD summary statistics from FinnGen and UKB meta-analysis
stopifnot(stab4$rsID == outcomeData[phenotype=="CAD",rsID])
stab4[,CAD_EAF := outcomeData[phenotype=="CAD",EAF]]
stab4[,CAD_nSamples := outcomeData[phenotype=="CAD",nSamples]]
stab4[,CAD_nCases := outcomeData[phenotype=="CAD",nCases]]
stab4[,CAD_logOR := outcomeData[phenotype=="CAD",beta]]
stab4[,CAD_se := outcomeData[phenotype=="CAD",se]]
stab4[,CAD_pval := outcomeData[phenotype=="CAD",pval]]

plot(stab4$LDLC_EAF,stab4$CAD_EAF)
abline(0,1)

#' add parental age at death summary statistics from UKB 
stopifnot(stab4$rsID == outcomeData[phenotype=="PAAD",rsID])
stab4[,ParentsAgeAtDeath_EAF := outcomeData[phenotype=="PAAD",EAF]]
stab4[,ParentsAgeAtDeath_nSamples := outcomeData[phenotype=="PAAD",nSamples]]
stab4[,ParentsAgeAtDeath_nCases := outcomeData[phenotype=="PAAD",nCases]]
stab4[,ParentsAgeAtDeath_beta := outcomeData[phenotype=="PAAD",beta]]
stab4[,ParentsAgeAtDeath_se := outcomeData[phenotype=="PAAD",se]]
stab4[,ParentsAgeAtDeath_pval := outcomeData[phenotype=="PAAD",pval]]

plot(stab4$LDLC_EAF,stab4$ParentsAgeAtDeath_EAF)
abline(0,1)

#' add flags to indicate usage in the MR_IVW approaches
#' 
stab4[,flag1 := F]
stab4[!is.element(rsID,c("rs2495491","rs693668")),flag1 := T]
stab4[,flag2 := F]
stab4[rsID %in% mySNPs_PCSK9,flag2 := T]
stab4[,flag3 := F]
stab4[rsID %in% mySNPs_HMGCR,flag3 := T]
stab4[,flag4 := F]
stab4[rsID %in% mySNPs_PCSK9_HMGCR,flag4 := T]

stab4[,phenotype := "LDLC_females"]
setnames(stab4,"phenotype","exposure")

#' # Sup Tab 5 ####
#' ***
#' MR-IVW results for LDLC
load("../results/09_MR_IVW_LDLC.RData")
MR_IVW[outcome == "BC", outcome := "BreastCancer [FinnGen+UKB]"]
MR_IVW[outcome == "BCS", outcome := "BCSurvival [BCAC]"]
MR_IVW[outcome == "CAD", outcome := "CAD [FinnGen+UKB]"]
MR_IVW[outcome == "PAAD", outcome := "ParentsAgeAtDeath [UKB]"]

stab5 = copy(MR_IVW)
stab5[,exposure := "LDLC_females"]

stab5[comment == "pruned",comment := "flag1"]
stab5[comment == "PCSK9",comment := "flag2"]
stab5[comment == "HMGCR",comment := "flag3"]
stab5[comment == "HMGCR_PCSK9",comment := "flag4"]
setnames(stab5,"comment","flag")

#' # Sup Tab 6 ####
#' ***
#' Instruments for MVMR-IVW for PCSK9 and LDLC
#' 
load("../temp/10_MVMR_IVW_input.RData")

stab6 = copy(LDLC)
names(stab6)[8:12] = paste0("LDLC_",names(stab6)[8:12])
stab6[,phenotype := "LDLC_females"]
setnames(stab6,"phenotype","exposure2")

#' add PCSK9 females summary statistics from Pott et al.
stab6[,exposure1 := "PCSK9_females"]
stab6[,PCSK9_EAF := PCSK9_females$EAF]
stab6[,PCSK9_nSamples := PCSK9_females$nSamples]
stab6[,PCSK9_beta := PCSK9_females$beta]
stab6[,PCSK9_se := PCSK9_females$SE]
stab6[,PCSK9_pval := PCSK9_females$pval]

stab6 = stab6[,c(1:6,13:18,7:12)]

#' add breast cancer survival summary statistics from BCAC
stopifnot(stab6$rsID == outcomeData[phenotype=="BCS",rsID])
stab6[,BCSurvival_EAF := outcomeData[phenotype=="BCS",EAF]]
stab6[,BCSurvival_nSamples := outcomeData[phenotype=="BCS",nSamples]]
stab6[,BCSurvival_nCases := outcomeData[phenotype=="BCS",nCases]]
stab6[,BCSurvival_logHR := outcomeData[phenotype=="BCS",beta]]
stab6[,BCSurvival_se := outcomeData[phenotype=="BCS",se]]
stab6[,BCSurvival_pval := outcomeData[phenotype=="BCS",pval]]

plot(stab6$LDLC_EAF,stab6$BCSurvival_EAF)
abline(0,1)

#' add breast cancer summary statistics from FinnGen and UKB meta-analysis
stopifnot(stab6$rsID == outcomeData[phenotype=="BC",rsID])
stab6[,BreastCancer_EAF := outcomeData[phenotype=="BC",EAF]]
stab6[,BreastCancer_nSamples := outcomeData[phenotype=="BC",nSamples]]
stab6[,BreastCancer_nCases := outcomeData[phenotype=="BC",nCases]]
stab6[,BreastCancer_logOR := outcomeData[phenotype=="BC",beta]]
stab6[,BreastCancer_se := outcomeData[phenotype=="BC",se]]
stab6[,BreastCancer_pval := outcomeData[phenotype=="BC",pval]]

plot(stab6$LDLC_EAF,stab6$BreastCancer_EAF)
abline(0,1)

#' add CAD summary statistics from FinnGen and UKB meta-analysis
stopifnot(stab6$rsID == outcomeData[phenotype=="CAD",rsID])
stab6[,CAD_EAF := outcomeData[phenotype=="CAD",EAF]]
stab6[,CAD_nSamples := outcomeData[phenotype=="CAD",nSamples]]
stab6[,CAD_nCases := outcomeData[phenotype=="CAD",nCases]]
stab6[,CAD_logOR := outcomeData[phenotype=="CAD",beta]]
stab6[,CAD_se := outcomeData[phenotype=="CAD",se]]
stab6[,CAD_pval := outcomeData[phenotype=="CAD",pval]]

plot(stab6$LDLC_EAF,stab6$CAD_EAF)
abline(0,1)

#' add parental age at death summary statistics from UKB 
stopifnot(stab6$rsID == outcomeData[phenotype=="PAAD",rsID])
stab6[,ParentsAgeAtDeath_EAF := outcomeData[phenotype=="PAAD",EAF]]
stab6[,ParentsAgeAtDeath_nSamples := outcomeData[phenotype=="PAAD",nSamples]]
stab6[,ParentsAgeAtDeath_nCases := outcomeData[phenotype=="PAAD",nCases]]
stab6[,ParentsAgeAtDeath_beta := outcomeData[phenotype=="PAAD",beta]]
stab6[,ParentsAgeAtDeath_se := outcomeData[phenotype=="PAAD",se]]
stab6[,ParentsAgeAtDeath_pval := outcomeData[phenotype=="PAAD",pval]]

plot(stab6$LDLC_EAF,stab6$ParentsAgeAtDeath_EAF)
abline(0,1)

#' add PCSK in statin-free individuals
#' 
stab6_2 = copy(stab6)
stab6_2[,exposure1 := "PCSK9_free"]
stab6_2[,PCSK9_EAF := PCSK9_free$EAF]
stab6_2[,PCSK9_nSamples := PCSK9_free$nSamples]
stab6_2[,PCSK9_beta := PCSK9_free$beta]
stab6_2[,PCSK9_se := PCSK9_free$SE]
stab6_2[,PCSK9_pval := PCSK9_free$pval]

stab6 = rbind(stab6,stab6_2)

#' add flag for different modes
mySNPs = c("rs2495491","rs11591147","rs693668","rs11583680","rs562556","rs12916","rs1525764")

stab6[,flag1 := F]
stab6[rsID %in% mySNPs[-5],flag1 := T]
stab6[,flag2 := F]
stab6[rsID %in% mySNPs[1:4],flag2 := T]
stab6[,flag3 := F]
stab6[rsID %in% mySNPs,flag3 := T]
stab6[,flag4 := F]
stab6[rsID %in% mySNPs[1:5],flag4 := T]

#' # Sup Tab 7 ####
#' ***
#' MVMR-IVW results for PCSK9 and LDLC
load("../results/10_MVMR_IVW_LDLC_PCSK9.RData")
MVMR_IVW[outcome == "BC", outcome := "BreastCancer [FinnGen+UKB]"]
MVMR_IVW[outcome == "BCS", outcome := "BCSurvival [BCAC]"]
MVMR_IVW[outcome == "PAAD", outcome := "ParentsAgeAtDeath [UKB]"]
MVMR_IVW[outcome == "CAD", outcome := "CAD [FinnGen+UKB]"]

stab7 = copy(MVMR_IVW)
stab7[,exposure2 := "LDLC_females"]
stab7[grepl("female",exposure1),exposure1 := "PCSK9_females"]
stab7[!grepl("female",exposure1),exposure1 := "PCSK9_free"]

stab7[comment == "HMGCR_PCSK9",comment := "flag1"]
stab7[comment == "HMGCR_PCSK9_rs562556",comment := "flag3"]
stab7[comment == "PCSK9",comment := "flag2"]
stab7[comment == "PCSK9_rs562556",comment := "flag4"]
setnames(stab7,"comment","flag")


#' # Save ####
#' ***
#' I do not want to change the names here, but change the order in which they appear in the excel sheet. 
#' 
#' First all tables with instrument information (ST1, ST3, ST4, and ST6).
#' 
#' Then the MR-IVW results (ST2 and ST5) and the MVMR-IVW results (ST7)
#' 
WriteXLS(x = c("stab1","stab3","stab4","stab6","stab2","stab5","stab7"), 
         ExcelFileName=paste0("../results/SupTables.xlsx"), 
         SheetNames=paste0("TableS",1:7), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(stab1,stab2,stab3,stab4,stab5,stab6,stab7,file = "../results/SupTables.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
