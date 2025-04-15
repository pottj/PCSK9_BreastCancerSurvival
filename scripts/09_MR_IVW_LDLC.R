#' ---
#' title: "MR of LDL-C on outcomes"
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
#' - exposure: LDL-C in women (GLGC data, Kanoni et al. 2022)
#' - outcome 1: breast cancer survival (BCAC data, Morra et al. 2021)
#' - outcome 2: breast cancer (FinnGen & UKB)
#' - outcome 3: longevity, parents' age at death (UKB, Pilling et al. 2017)
#' - outcome 4: CAD (FinnGen & UKB)
#'
#' Instrument selection: 
#' 
#' - genome-wide SNPs, priority pruning by position (exclude SNPs within +- 1Mb of variant with lowest p-value)
#' - same SNPs as in PCSK9 pQTL MR approach (3 SNPs, as one variant is missing in LDLC data set)
#' - HMGCR SNPs
#' - HMGCR and PCSK9 SNPs
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
load("../temp/08_MR_LDLC_input_Pruned.RData")
mySNPs_pruned = LDLC$rsID

load("../temp/08_MR_LDLC_input_beforePruning.RData")
mySNPs_PCSK9 = c("rs2495491","rs11591147","rs693668","rs11583680")
mySNPs_HMGCR = c("rs12916","rs1525764")
mySNPs_PCSK9_HMGCR = c(mySNPs_PCSK9,mySNPs_HMGCR)

modes = c("pruned","PCSK9","HMGCR","HMGCR_PCSK9")
modes_SNPlist = list(mySNPs_pruned,mySNPs_PCSK9,mySNPs_HMGCR,mySNPs_PCSK9_HMGCR)
OutcomeList = list(BCAC,BC_FinnGen_UKB,CAD,ParentalLongevity_Death)

#' # Loop ####
#' ***
dumTab1 = foreach(i = 1:length(modes))%do%{
  #i=1
  message("Working on mode ",i," (",modes[i],")")
  
  data_x = copy(LDLC)
  instruments = modes_SNPlist[[i]]
  data_x = data_x[rsID %in% instruments,]
  
  dumTab2 = foreach(j = 1:4)%do%{
    #j=1
    message("   Working on outcome ",j)
    
    data_y = OutcomeList[[j]]
    data_y = data_y[rsID %in% instruments,]
    
    myOutcome = unique(data_y$phenotype)
    myOutcome2 = gsub(" ","",myOutcome)
    
    # running MR
    mrob = mr_input(bx = data_x$beta, 
                    bxse = data_x$se,
                    by = data_y$beta, 
                    byse = data_y$se,
                    snps = paste0(data_x$rsID,"\n"), 
                    exposure = "LDLC",
                    outcome = myOutcome)
    mod1= mr_ivw(mrob)
    
    mr_plot(mrob, interactive = F,labels = T)
    ggsave(paste0("../results/09_ScatterPlot_IVW/LDLC_",myOutcome2,"_",modes[i],".png"),
           height = 7, width = 7)
    
    res = data.table(exposure = mod1@Exposure,
                     outcome = mod1@Outcome,
                     nSNPs = dim(data_x)[1],
                     beta_IVW = mod1@Estimate,
                     se_IVW = mod1@StdError,
                     pval_IVW = mod1@Pvalue,
                     FStat = mod1@Fstat, 
                     Q = mod1@Heter.Stat[1],
                     pval_Q = mod1@Heter.Stat[2])
    res
    
  }
  
  MRTab = rbindlist(dumTab2)
  MRTab[,comment := modes[i]]
  MRTab
}
MR_IVW = rbindlist(dumTab1)
MR_IVW[pval_IVW <0.05]

#' # Save data ####
#' ***
save(MR_IVW, file = "../results/09_MR_IVW_LDLC.RData")

mySNPs = unique(c(mySNPs_pruned,mySNPs_PCSK9,mySNPs_HMGCR))
LDLC = LDLC[rsID %in% mySNPs,]
outcomeData = rbind(BCAC, BC_FinnGen_UKB, ParentalLongevity_Death, CAD,fill=T)
outcomeData = outcomeData[rsID %in% mySNPs,]

save(LDLC,outcomeData,file = "../temp/09_MR_IVW_input.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
