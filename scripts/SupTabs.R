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
#' Sup Tab 1: all instruments for MR approach (one line per SNP - exposure combination)
#' 
#' Sup Tab 2: all instruments for the MVMR approach (one line per SNP - exposure combination)
#' 
#' Sup Tab 3: all summary statistics for all outcomes (one line per SNP - outcome combination)
#' 
#' Sup Tab 4: MR-IVW results & MR-ratio results (no rs562556)
#' 
#' Sup Tab 5: MR-ratio results for rs562556
#' 
#' Sup Tab 6: MVMR-IVW results
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

#' # Sup Tab 1 ####
#' ***
#' All instruments for MR approach (one line per SNP - exposure combination)
#' 
load("../results/03_Exposure_for_MR_pruned.RData")
ExposureData

stab1 = copy(ExposureData)
setnames(stab1,"phenotype","exposure")
setorder(stab1,-exposure,setting,chr,pos_b37)

#' # Sup Tab 2 ####
#' ***
#' All instruments for the MVMR approach (one line per SNP - exposure combination)
#' 
load("../results/05_MVMR_input.RData")
ExposureData_MVMR

stab2 = copy(ExposureData_MVMR)

#' Add some missing columns (Z-Score, position in b38)
stab2[,exp1_absZScore := abs(exp1_beta/exp1_se)]
stab2[,exp2_absZScore := abs(exp2_beta/exp2_se)]
matched = match(stab2$rsID,stab1$rsID)
stopifnot(is.na(matched)==F)
stab2[,pos_b38 := stab1[matched,pos_b38]]

#' Add flags to indicate which SNP was used for what analysis
stab2[,flag1_PCSK9_v1 := F]
stab2[rsID %in% stab1[exposure == "PCSK9 protein levels" & !grepl("ratio",setting),rsID],flag1_PCSK9_v1 := T]
stab2[,flag2_PCSK9_v2 := F]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting == "females PCSK9 SNPs",rsID] & setting == "females",flag2_PCSK9_v2 := T]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting == "all PCSK9 SNPs",rsID] & setting == "all",flag2_PCSK9_v2 := T]
stab2[,flag3_PCSK9_HMGCR := F]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting %in% c("females PCSK9 SNPs","females HMGCR SNPs"),rsID] & setting == "females",flag3_PCSK9_HMGCR := T]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting %in% c("all PCSK9 SNPs","all HMGCR SNPs"),rsID] & setting == "all",flag3_PCSK9_HMGCR := T]
stab2[,flag4_genome_wide := F]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting %in% c("females gw SNPs"),rsID] & setting == "females",flag4_genome_wide := T]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & setting %in% c("all gw SNPs"),rsID] & setting == "all",flag4_genome_wide := T]
stab2 = stab2[flag1_PCSK9_v1==T | flag2_PCSK9_v2==T | flag3_PCSK9_HMGCR==T | flag4_genome_wide==T,]

#' Change order 
names(stab2)
stab2 = stab2[,c(1:3,21,4:12,19,13:18,20,22:25)]
setorder(stab2,-setting,chr,pos_b37)

#' # Sup Tab 3 ####
#' ***
#' All summary statistics for all outcomes (one line per SNP - outcome combination)
#' 
load("../results/02_Outcome_for_MR.RData")
OutcomeData

stab3 = copy(OutcomeData)
stab3[,chr := as.numeric(chr)]
stab3[,absZScore := abs(beta/se)]
setnames(stab3,"phenotype","outcome")
setorder(stab3,outcome,setting,chr,pos_b37)

stab3[setting == "all", setting := "Aragam et al. sex-combined"]
stab3[setting == "females", setting := "Aragam et al. females"]
stab3[setting == "UKB", setting := "Pilling et al. sex-combined"]
stab3[grepl("BC",outcome), setting := paste(setting,"females")]
stab3[,setting := gsub("BCAC","Morra et al.",setting)]

stab3[,beta_type := "linear effect"]
stab3[outcome %in% c("BC","CAD"),beta_type := "log(OR)"]
stab3[outcome %in% c("BCS"),beta_type := "log(HR)"]

stab3[outcome == "CAD",outcome := "Coronary Artery Disease"]
stab3[outcome == "BC",outcome := "Breast Cancer"]
stab3[outcome == "BCS",outcome := "Breast Cancer Survival"]
stab3[outcome == "PLD",outcome := "Parental Longevity (combined parental age at death)"]

#' Filter for SNPs in the exposure data sets
stab3 = stab3[rsID %in% stab1$rsID,]
stopifnot(stab2$rsID %in% stab3$rsID)

#' # Sup Tab 4 ####
#' ***
#' MR-IVW results & MR-ratio results (no rs562556)
#' 
load("../results/04_MR.RData")
MRTab

MRTab[outcome == "BC - FinnGen + UKB", outcome := "Breast Cancer - females"]
MRTab[outcome == "BCS - BCAC", outcome := "Breast Cancer Survival - females - Morra et al."]
MRTab[outcome == "BCS - Mei et al.", outcome := "Breast Cancer Survival - females - Mei at al."]
MRTab[outcome == "CAD - all", outcome := "Coronary Artery Disease - all"]
MRTab[outcome == "CAD - females", outcome := "Coronary Artery Disease - females"]
MRTab[outcome == "PLD - UKB", outcome := "Parental Longevity - all"]

stab4 = copy(MRTab)
stab4 = stab4[!grepl("ratio",exposure)]

#' # Sup Tab 5 ####
#' ***
#' MR-ratio results for rs562556
#' 
stab5 = copy(MRTab)
stab5 = stab5[grepl("ratio",exposure)]

stab5 = stab5[,c(1,2,5:8)]
stab5[,exposure := gsub(" ratio","",exposure)]

#' # Sup Tab 6 ####
#' ***
#' MVMR-IVW results
#' 
load("../results/05_MVMR.RData")
MVMRTab

stab6 = copy(MVMRTab)
stab6[outcome == "BC - FinnGen + UKB", outcome := "Breast Cancer - females"]
stab6[outcome == "BCS - BCAC", outcome := "Breast Cancer Survival - females - Morra et al."]
stab6[outcome == "CAD - all", outcome := "Coronary Artery Disease - all"]
stab6[outcome == "CAD - females", outcome := "Coronary Artery Disease - females"]
stab6[outcome == "PLD - UKB", outcome := "Parental Longevity - all"]

setnames(stab6,"sex","exposure_sex")
setnames(stab6,"SNPset","flag")
stab6[flag == "PCSK9_v1", flag := "flag1_PVSK9_v1"]
stab6[flag == "PCSK9_v2", flag := "flag2_PCSK9_v2"]
stab6[flag == "PCSK9_HMGCR", flag := "flag3_PCSK9_HMGCR"]
stab6[flag == "genome-wide", flag := "flag4_genome_wide"]

stab6 = stab6[,c(1,3,4,2,5:16)]
setnames(stab6,"NR_SNPs_total","nSNPs")
setnames(stab6,"SE_exp1","se_exp1")
setnames(stab6,"SE_exp2","se_exp2")
setnames(stab6,"condFstat_exp1","condFStat_exp1")
setnames(stab6,"condFstat_exp2","condFStat_exp2")
setnames(stab6,"HeteroStat","Q")
setnames(stab6,"HeteroStat_pval","pval_Q")

#' # Save ####
#' ***
#' 
WriteXLS(x = c("stab1","stab2","stab3","stab4","stab5","stab6"), 
         ExcelFileName=paste0("../results/SupTables.xlsx"), 
         SheetNames=paste0("TableS",1:6), 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

save(stab1,stab2,stab3,stab4,stab5,stab6,file = "../results/SupTables.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
