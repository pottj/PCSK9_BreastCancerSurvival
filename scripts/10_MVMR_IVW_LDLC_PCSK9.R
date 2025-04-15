#' ---
#' title: "MVMR of LDL-C and PCSK9 on outcomes"
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
#' Used data sources: 
#' 
#' - exposure 1: LDL-C in women (GLGC data, Kanoni et al. 2022)
#' - exposure 2: PCSK9 in women (Pott et al. 2024)
#' - exposure 3: PCSK9 in statin-free individuals (Pott et al. 2024)
#' - outcome 1: breast cancer survival (BCAC data, Morra et al. 2021)
#' - outcome 2: breast cancer (FinnGen & UKB)
#' - outcome 3: longevity, parents' age at death (UKB, Pilling et al. 2017)
#' - outcome 4: CAD (FinnGen & UKB)
#'
#' Instrument selection: 
#' 
#' - HMGCR and PCSK9 SNPs --> 5 SNPs
#' - PCSK9 SNPs only --> 3 SNPs
#' - HMGCR and PCSK9 SNPs + rs562556 --> 6 SNPs
#' - PCSK9 SNPs only + rs562556 --> 4 SNPs
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
load("../temp/08_MR_LDLC_input_beforePruning.RData")

mySNPs = c("rs2495491","rs11591147","rs693668","rs11583680","rs562556","rs12916","rs1525764")

LDLC = LDLC[rsID %in% mySNPs,]
BCAC = BCAC[rsID %in% mySNPs,]
BC_FinnGen_UKB = BC_FinnGen_UKB[rsID %in% mySNPs,]
CAD = CAD[rsID %in% mySNPs,]
ParentalLongevity_Death = ParentalLongevity_Death[rsID %in% mySNPs,]

#' # Get PCSK9 
myFiles = list.files(path = GWAS_data_PCSK9_pQTLs)
myFiles = myFiles[!grepl("gz",myFiles)]
myFiles = myFiles[c(1,4)]
myFiles

PCSK9_females = fread(paste0(GWAS_data_PCSK9_pQTLs,myFiles[1]))
PCSK9_females[, rsID := gsub(":.*","",markername)]
PCSK9_females = PCSK9_females[rsID %in% LDLC$rsID]
stopifnot(PCSK9_females$rsID == LDLC$rsID)

PCSK9_free = fread(paste0(GWAS_data_PCSK9_pQTLs,myFiles[2]))
PCSK9_free[, rsID := gsub(":.*","",markername)]
PCSK9_free = PCSK9_free[rsID %in% LDLC$rsID]
stopifnot(PCSK9_free$rsID == LDLC$rsID)

#' Harmonize alleles (get same coding as in LDLC data set)
#' 
filt = LDLC$EA != PCSK9_females$EA
PCSK9_free[filt,EAF := 1-EAF]
PCSK9_free[filt,beta := beta * (-1)]
PCSK9_females[filt,EAF := 1-EAF]
PCSK9_females[filt,beta := beta * (-1)]
plot(LDLC$EAF,PCSK9_females$EAF)
abline(0,1)
diff = abs(LDLC$EAF - PCSK9_females$EAF)
hist(diff)

#' Get necessary parameters
modes = c("HMGCR_PCSK9","PCSK9","HMGCR_PCSK9_rs562556","PCSK9_rs562556")
modes_SNPlist = list(mySNPs[-5],mySNPs[1:4],mySNPs,mySNPs[1:5])
OutcomeList = list(BCAC,BC_FinnGen_UKB,CAD,ParentalLongevity_Death)

#' # Loop for PCSK9 females ####
#' ***
dumTab1 = foreach(i = 1:length(modes))%do%{
  #i=1
  message("Working on mode ",i," (",modes[i],")")
  
  data_x1 = copy(LDLC)
  data_x2 = copy(PCSK9_females)
  data_x3 = copy(PCSK9_free)
  
  instruments = modes_SNPlist[[i]]
  data_x1 = data_x1[rsID %in% instruments,]
  data_x2 = data_x2[rsID %in% instruments,]
  data_x3 = data_x3[rsID %in% instruments,]
  
  dumTab2 = foreach(j = 1:4)%do%{
    #j=1
    message("   Working on outcome ",j)
    
    data_y = OutcomeList[[j]]
    data_y = data_y[rsID %in% instruments,]
    
    myOutcome = unique(data_y$phenotype)
    myOutcome2 = gsub(" ","",myOutcome)
    
    # running MR
    data_beta1 = cbind(data_x2$beta,data_x1$beta)
    data_SE1 = cbind(data_x2$SE,data_x1$se)
    data_beta2 = cbind(data_x3$beta,data_x1$beta)
    data_SE2 = cbind(data_x3$SE,data_x1$se)
    
    mvmr_obj1 = mr_mvinput(bx = as.matrix(data_beta1),
                          bxse = as.matrix(data_SE1),
                          by = data_y$beta, 
                          byse = data_y$se,
                          snps = data_y$rsID,
                          exposure = c("PCSK9 (females)","LDLC"),
                          outcome = myOutcome)
    mvmr_obj2 = mr_mvinput(bx = as.matrix(data_beta2),
                           bxse = as.matrix(data_SE2),
                           by = data_y$beta, 
                           byse = data_y$se,
                           snps = data_y$rsID,
                           exposure = c("PCSK9 (statin-free)","LDLC"),
                           outcome = myOutcome)
    n1 = data_x1[,max(nSamples)]
    n2 = data_x2[,max(nSamples)]
    n3 = data_x3[,max(nSamples)]
    
    res1 = mr_mvivw(mvmr_obj1,nx = c(n2,n1))
    res2 = mr_mvivw(mvmr_obj2,nx = c(n3,n1))
    
    res3 = data.table(outcome = res1@Outcome, 
                      NR_SNPs_total = res1@SNPs,
                      exposure1 = res1@Exposure[1],
                      beta_exp1 = res1@Estimate[1],
                      SE_exp1 = res1@StdError[1],
                      pval_exp1 = res1@Pvalue[1],
                      condFstat_exp1 = res1@CondFstat[1],
                      exposure2 = res1@Exposure[2],
                      beta_exp2 = res1@Estimate[2],
                      SE_exp2 = res1@StdError[2],
                      pval_exp2 = res1@Pvalue[2],
                      condFstat_exp2 = res1@CondFstat[2],
                      HeteroStat = res1@Heter.Stat[1],
                      HeteroStat_pval = res1@Heter.Stat[2])

    res4 = data.table(outcome = res2@Outcome, 
                      NR_SNPs_total = res2@SNPs,
                      exposure1 = res2@Exposure[1],
                      beta_exp1 = res2@Estimate[1],
                      SE_exp1 = res2@StdError[1],
                      pval_exp1 = res2@Pvalue[1],
                      condFstat_exp1 = res2@CondFstat[1],
                      exposure2 = res2@Exposure[2],
                      beta_exp2 = res2@Estimate[2],
                      SE_exp2 = res2@StdError[2],
                      pval_exp2 = res2@Pvalue[2],
                      condFstat_exp2 = res2@CondFstat[2],
                      HeteroStat = res2@Heter.Stat[1],
                      HeteroStat_pval = res2@Heter.Stat[2])
    
    res = rbind(res3,res4)
    res

  }
  
  MVMRTab = rbindlist(dumTab2)
  MVMRTab[,comment := modes[i]]
  MVMRTab
}
MVMR_IVW = rbindlist(dumTab1)
MVMR_IVW[pval_exp1 <0.05]

#' # Save data ####
#' ***
save(MVMR_IVW, file = "../results/10_MVMR_IVW_LDLC_PCSK9.RData")

LDLC = LDLC[rsID %in% mySNPs,]
outcomeData = rbind(BCAC, BC_FinnGen_UKB, ParentalLongevity_Death, CAD,fill=T)
outcomeData = outcomeData[rsID %in% mySNPs,]

save(LDLC,PCSK9_females,PCSK9_free,outcomeData,file = "../temp/10_MVMR_IVW_input.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
