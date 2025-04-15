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
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = list.files(path = GWAS_data_PCSK9_pQTLs)
myFiles
myFiles = myFiles[!grepl("gz",myFiles)]
myFiles = myFiles[c(1,4)]
myFiles

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data0 = fread(paste0(GWAS_data_PCSK9_pQTLs,myFiles[i]))
  data0 = data0[chr == 1,]
  data0 = data0[bp_hg19 < 55505647 + 500000]
  data0 = data0[bp_hg19 > 55505647 - 500000]
  data0
}
pQTLData = rbindlist(dumTab1)
pQTLData[,table(phenotype,pval<5e-8)]
pQTLData[,table(phenotype,pval<1e-6)]
pQTLData[,min(pval,na.rm = T),by=phenotype]
pQTLData = pQTLData[!is.na(pval),]
pQTLData[, rsID := gsub(":.*","",markername)]
save(pQTLData, file = "../temp/PCSK9_pQTLs.RData")

pQTLData[ rsID == "rs562556",]
pQTLData = pQTLData[rsID == "rs562556",]
pQTLData

load("../temp/GTExV8_selectedTissues_PCSK9.RData")
GTExData = GTExData[pos_b37 %in% pQTLData$bp_hg19,]
GTExData = GTExData[pval<0.05,]
GTExData

#' # Format ####
#' ***
#' I want the same format as before
matched = match(pQTLData$bp_hg19,GTExData$pos_b37)
table(is.na(matched))
pQTLData[,pos_b38 := GTExData[matched,pos_b38]]
names(pQTLData)
pQTLData = pQTLData[, c(17,2,3,18,5,4, 16,6,8,10,11,12)]
head(pQTLData)
names(pQTLData) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                    "phenotype","EAF","nSamples","beta","se","pval")

matched = match(GTExData$pos_b37,pQTLData$pos_b37)
GTExData[,rsID := pQTLData[matched,rsID]]
GTExData[,EA_pQTLs := pQTLData[matched,EA]]
GTExData[,EAF_pQTLs := pQTLData[matched,EAF]]
filt = GTExData$effect_allele == GTExData$EA_pQTLs
table(filt)
GTExData[,eaf := maf]
GTExData[EAF_pQTLs>0.5,eaf := 1-maf]

names(GTExData)
GTExData = GTExData[, c(24,10,17,11,12,13, 22,27,14,7,8,6)]
head(GTExData)
names(GTExData) = names(pQTLData)

IVData = rbind(GTExData,pQTLData)
IVData[pval<1e-6,]
IVData[pval<0.05,]

#' # Load outcome data
#' ***
load("../results/03_outcome_harmonized.RData")
outcomeData = outcomeData[rsID == "rs562556",]

#' Save data 
save(IVData, outcomeData,file = "../results/06_SumStats_rs562556.RData")

#' # Wald ratio ####
#' ***
myTraits = unique(IVData$phenotype)
myOutcomes = unique(outcomeData$phenotype)

dumTab2 = foreach(i = 1:length(myOutcomes))%do%{
  #i=1
  myRow = outcomeData[i,]
  
  dumTab3 = foreach(j = 1:length(myTraits))%do%{
    #j=1
    myX = IVData[j,]
    
    # get ratio and SE
    test = MRfunction_jp(betaX = myX$beta, seX = myX$se, betaY = myRow$beta, seY = myRow$se)
    
    # prep output (SNP info, statistics exposure, statistics outcome, Wald ratio)
    res = cbind(myX, myRow[,7:13],test)
    res
  }
  tab3 = rbindlist(dumTab3)
  tab3
}
MR_ratio = rbindlist(dumTab2)
MR_ratio[pval<0.05,table(p_IV2<0.05)]
names(MR_ratio)[c(13:19)] = paste("outcome",names(MR_ratio)[c(13:19)],sep="_")

#' # Save data ####
#' ***
save(MR_ratio, file = "../results/06_MR_ratio_rs562556.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

