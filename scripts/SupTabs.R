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

#' # Sup Tab 0 ####
#' ***
#' Content
stab0 = data.table(Table = paste0("S",1:6),
                   Title = c("Instruments and exposure association for the MR approaches",
                             "Instruments and expoures association for the MVMR approach",
                             "Instruments and outcome associations for both MR and MVMR approaches",
                             "Results of the MR-ratio using rs562556", 
                             "Results of the MR analyses using genome-wide significant instrumets",
                             "Results of the MVMR analyses"))

#' # Sup Tab 1 ####
#' ***
#' All instruments for MR approach (one line per SNP - exposure combination)
#' 
load("../results/03b_Exposure_for_MR_filtered.RData")
ExposureData

stab1 = copy(ExposureData)
setnames(stab1,"phenotype","exposure")
setorder(stab1,-exposure,setting,chr,pos_b37)
names(stab1)
stab1 = stab1[,c(1:10,18,11:16)]
head(stab1)

#' # Sup Tab 2 ####
#' ***
#' All instruments for the MVMR approach (one line per SNP - exposure combination)
#' 
load("../results/05_MVMR_input.RData")
ExposureData_MVMR

stab2 = copy(ExposureData_MVMR)
names(stab2)

#' Add some missing columns (Z-Score, position in b38)
stab2[,exp1_absZScore := abs(exp1_beta/exp1_se)]
stab2[,exp2_absZScore := abs(exp2_beta/exp2_se)]
matched = match(stab2$rsID,stab1$rsID)
stopifnot(is.na(matched)==F)
stab2[,pos_b38 := stab1[matched,pos_b38]]
stab2[,exposure1 := gsub("protein","PE",exposure1)]

#' Add flags to indicate which SNP was used for what analysis
stab2[,flag1_PCSK9 := F]
stab2[rsID %in% stab1[exposure == "PCSK9 PE levels" & MRapproach == "QTL",rsID],flag1_PCSK9 := T]
stab2[,table(flag1_PCSK9)]

stab2[,flag2_PCSK9_HMGCR := F]
stab2[rsID %in% stab1[exposure == "LDL-C levels" & MRapproach == "QTL" & instrumentLocus != "genome-wide",rsID],flag2_PCSK9_HMGCR := T]
stab2[,table(flag2_PCSK9_HMGCR)]

#' Change order 
names(stab2)
stab2 = stab2[,c(1:3,21,4:12,19,13:18,20,22:23)]
setorder(stab2,-setting,chr,pos_b37)

#' # Sup Tab 3 ####
#' ***
#' All summary statistics for all outcomes (one line per SNP - outcome combination)
#' 
load("../results/03b_Outcome_for_MR_filtered.RData")
OutcomeData

stab3 = copy(OutcomeData)
stab3[,chr := as.numeric(chr)]
stab3[,absZScore := abs(beta/se)]
setnames(stab3,"phenotype","outcome")
setorder(stab3,outcome,setting,chr,pos_b37)

names(stab3)
stab3 = stab3[,c(1:8,10,11:16,9,18)]
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

stab3[,source := gsub("Mei et al.","meta-analysis",source)]

#' # Sup Tab 4 ####
#' ***
#' MR-IVW results 
#' 
load("../results/04_MR.RData")
MRTab

dummy1 = unlist(strsplit(MRTab$exposure," - "))
MRTab[,exposursSetting := dummy1[seq(2,length(dummy1),4)]]
MRTab[,exposureInstruments := dummy1[seq(3,length(dummy1),4)]]
MRTab[,exposureApproach:= dummy1[seq(4,length(dummy1),4)]]
MRTab[,exposure := dummy1[seq(1,length(dummy1),4)]]

dummy2 = unlist(strsplit(MRTab$outcome," - "))
MRTab[,outcomeSex := dummy2[seq(2,length(dummy2),4)]]
MRTab[,outcomeSource := dummy2[seq(3,length(dummy2),4)]]
MRTab[,outcomeModel:= dummy2[seq(4,length(dummy2),4)]]
MRTab[,outcome := dummy2[seq(1,length(dummy2),4)]]

MRTab[,pval_adj := p.adjust(p=pval,method = "fdr")]
MRTab[outcome == "CAD",outcome := "Coronary Artery Disease"]
MRTab[outcome == "BC",outcome := "Breast Cancer"]
MRTab[outcome == "BCS",outcome := "Breast Cancer Survival"]
MRTab[outcome == "PLD",outcome := "Parental Longevity (combined parental age at death)"]

names(MRTab)
MRTab = MRTab[, c(1,11:13,2,14:16,3:7,17,8:10)]

#' S4: MR-ratio table
stab4 = copy(MRTab)
stab4 = stab4[exposureApproach == "rs562556",]

#' remove columns that are always the same (or NA)
#' 
stab4[,exposureInstruments := NULL]
stab4[,exposureApproach := NULL]
stab4[,nSNPs := NULL]
stab4[,method := NULL]
stab4[,Q := NULL]
stab4[,pval_Q := NULL]
stab4[,outcomeSource := gsub("Mei et al.","meta-analysis",outcomeSource)]
stab4

#' S5: MR-IVW table
stab5 = copy(MRTab)
stab5 = stab5[exposureApproach != "rs562556",]

#' remove columns that are always the same (or NA)
#' 
stab5[,exposureApproach := NULL]
stab5[,outcomeModel := NULL]
stab5

#' # Sup Tab 6 ####
#' ***
#' MVMR-IVW results
#' 
load("../results/05_MVMR.RData")
MVMRTab
MVMRTab[,pval_adj_exp1 := p.adjust(p=pval_exp1,method = "fdr")]
MVMRTab[,pval_adj_exp2 := p.adjust(p=pval_exp2,method = "fdr")]

dummy2 = unlist(strsplit(MVMRTab$outcome," - "))
MVMRTab[,outcomeSex := dummy2[seq(2,length(dummy2),4)]]
MVMRTab[,outcomeSource := dummy2[seq(3,length(dummy2),4)]]
MVMRTab[,outcome := dummy2[seq(1,length(dummy2),4)]]
MVMRTab[outcome == "CAD",outcome := "Coronary Artery Disease"]
MVMRTab[outcome == "BC",outcome := "Breast Cancer"]
MVMRTab[outcome == "BCS",outcome := "Breast Cancer Survival"]
MVMRTab[outcome == "PLD",outcome := "Parental Longevity (combined parental age at death)"]

names(MVMRTab)

stab6 = copy(MVMRTab)
stab6 = stab6[,c(1,19,20,2:8,17,9:13,18,14:16)]
setnames(stab6,"sex","exposureSex")
setnames(stab6,"SNPset","instrumentLocus")
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
myList = list(stab0,stab1,stab2,stab3,stab4,stab5,stab6)
write_xlsx(x = myList, 
           path = paste0("../results/SupTables.xlsx"),col_names = T,format_headers= F)
save(stab1,stab2,stab3,stab4,stab5,stab6,file = "../results/SupTables.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
