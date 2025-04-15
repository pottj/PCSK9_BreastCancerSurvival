#' ---
#' title: "Test candidate SNP only"
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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = T

source("../SourceFile.R")
.libPaths()

#' # Get data from Mei et al. ####
#' ***
#' I extracted manually the hazard ratios from the publication. All numbers are taken from the Supplemental Figure S1, panels as indicated per entry. 
#' 
myTab = read_xlsx("../temp/Mei_2024_rs561556_BCS.xlsx")
setDT(myTab)
myTab

#' # Meta-analyses ####
#' ***
#' There are 3 types of HR: one for all, one in EUR only, and one in data adjusted for population structure. In the main MR, I only want the EUR only, but I will meta-analyse all of them.
#' 
#' ## All samples
#' 
filt1 = myTab$comment == "all" #& myTab$Study != "Nik-Zainal et al., 2016"

meta1 = metagen(TE = myTab[filt1,log(HR)], lower = myTab[filt1,log(CI_low)],
                upper = myTab[filt1,log(CI_up)], sm = "HR", studlab = myTab[filt1,Study])
summary(meta1)
forest(meta1)

#' ## EUR only
#' 
filt2 = myTab$comment == "EUR only" #& myTab$Study != "Nik-Zainal et al., 2016"

meta2 = metagen(TE = myTab[filt2,log(HR)], lower = myTab[filt2,log(CI_low)],
                upper = myTab[filt2,log(CI_up)], sm = "HR", studlab = myTab[filt2,Study])
summary(meta2)
forest(meta2)

#' ## Adjusted for population structure
#' 
filt3 = myTab$comment == "PopStruc adjusted" #& myTab$Study != "Nik-Zainal et al., 2016"

meta3 = metagen(TE = myTab[filt3,log(HR)], lower = myTab[filt3,log(CI_low)],
                upper = myTab[filt3,log(CI_up)], sm = "HR", studlab = myTab[filt3,Study])
summary(meta3)
forest(meta3)

#' ## Combine results
outcomeData_recessive = data.table(setting = c("all","EUR only","PopStruc adjusted"), 
                 logHR = c(meta1$TE.fixed,meta2$TE.fixed,meta3$TE.fixed), 
                 SE_logHR = c(meta1$seTE.fixed,meta2$seTE.fixed,meta3$seTE.fixed),
                 pval = c(meta1$pval.fixed,meta2$pval.fixed,meta3$pval.fixed))
outcomeData_recessive

#' # Load instruments ####
#' ***
#' I reuse the data set I generated in the previous script. 
load("../results/06_SumStats_rs562556.RData")
IVData

#' # Wald ratio ####
#' ***
myTraits = unique(IVData$phenotype)
myOutcomes = unique(outcomeData_recessive$setting)

dumTab2 = foreach(i = 1:length(myOutcomes))%do%{
  #i=1
  myRow = outcomeData_recessive[i,]
  
  dumTab3 = foreach(j = 1:length(myTraits))%do%{
    #j=1
    myX = IVData[j,]
    
    # get ratio and SE
    test = MRfunction_jp(betaX = myX$beta, seX = myX$se, betaY = myRow$logHR, seY = myRow$SE_logHR)
    
    # prep output (SNP info, statistics exposure, statistics outcome, Wald ratio)
    res = cbind(myX, myRow,test)
    res
  }
  tab3 = rbindlist(dumTab3)
  tab3
}
MR_ratio = rbindlist(dumTab2)
MR_ratio[setting=="EUR only" & p_IV2<0.05]

#' # Save data ####
#' ***
save(MR_ratio, file = "../results/07_MR_ratio_rs562556_Mei.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
