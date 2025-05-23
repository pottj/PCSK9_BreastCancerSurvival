#' ---
#' title: "MVMR: PCSK9 and LDLC on outcomes"
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
load_meta = T

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/03_Exposure_for_MR_pruned.RData")
load("../results/02_Outcome_for_MR.RData")
OutcomeData = OutcomeData[rsID %in% ExposureData$rsID,]
OutcomeData[,chr := as.numeric(chr)]

ExposureData = ExposureData[!grepl("gene expression",phenotype)]
ExposureData = ExposureData[!grepl("ratio",setting)]
ExposureData[,.N,by=c("phenotype","setting")]

SNPList = copy(ExposureData)
SNPList = SNPList[,c(1:6)]
SNPList = distinct(SNPList)
setDT(SNPList)
SNPList[,chrPos := paste(chr,pos_b37,sep = "_")]

#' # Load PCSK9 data ####
#' ***
myFiles = c(Pott_PCSK9_females,Pott_PCSK9_males)

dumTab2 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data0 = fread(myFiles[i])
  
  # filter data for relevant SNPs
  data0[,rsID := gsub(":.*","",markername)]
  data0 = data0[rsID %in% ExposureData$rsID,]
  data0[,setting := gsub("PCSK9_","",phenotype)]
  data0[,phenotype := "PCSK9 protein levels"]
  
  # select relevant columns
  names(data0)
  data0 = data0[,c(17,2,3,5,4,16,18,6,8,10,11,12)]
  names(data0) = c("rsID","chr","pos_b37","OA","EA",
                   "phenotype","setting","EAF","nSamples","beta","se","pval")
  
  # return
  data0
}
pQTLData = rbindlist(dumTab2)

#' Meta-analyse males and females
mySNPs = unique(pQTLData$rsID)
dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  filt = pQTLData$rsID == mySNPs[i]
  meta = metagen(TE = pQTLData[filt,beta],
                 seTE = pQTLData[filt,se],
                 studlab = pQTLData[filt,setting])
  # summary(meta)
  # forest(meta)
  
  mySamples = pQTLData[filt,nSamples]
  myEAFs = pQTLData[filt,EAF]
  mySamples2 = mySamples/sum(mySamples)
  myEAFs2 = round(sum(myEAFs * mySamples/sum(mySamples)),4)
  
  data = copy(pQTLData)
  data = data[!grepl("females",setting) & rsID == mySNPs[i],]
  data[, beta := meta$TE.fixed]
  data[, se := meta$seTE.fixed]
  data[, pval := meta$pval.fixed]
  data[, nSamples := sum(mySamples)]
  data[, EAF := myEAFs2]
  data[, setting := "all"]
  
  data
}
pQTLData2 = rbindlist(dumTab2)
pQTLData = rbind(pQTLData[grepl("females",setting),],pQTLData2)
setorder(pQTLData,setting,chr,pos_b37)

matched = match(pQTLData$rsID,SNPList$rsID)
table(is.na(matched))
table(SNPList[matched,rsID] == pQTLData$rsID)
table(SNPList[matched,EA] == pQTLData$EA)
pQTLData[,EA2 := SNPList[matched,EA]]
pQTLData[,OA2 := SNPList[matched,OA]]
pQTLData[EA != EA2, EAF := 1-EAF]
pQTLData[EA != EA2, beta := beta*(-1)]
pQTLData[, EA := EA2]
pQTLData[, OA := OA2]
pQTLData[, EA2 := NULL]
pQTLData[, OA2 := NULL]

#' # Load LDLC data ####
#' ***
myFiles = c(GLGC_LDLC_females,GLGC_LDLC_all)
mySettings = c("females","all")

dumTab3 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data0 = fread(myFiles[i])
  
  # filter data for relevant SNPs
  data0 = data0[!is.na(rsID),]
  data0 = data0[POOLED_ALT_AF>0.01,]
  data0 = data0[rsID %in% pQTLData$rsID,]
  data0[,setting := mySettings[i]]
  data0[,phenotype := "LDL-C levels"]
  
  # select relevant columns
  names(data0)
  data0 = data0[,c(1,2:5,16,15,8,6,9,10,12)]
  names(data0) = c("rsID","chr","pos_b37","OA","EA",
                   "phenotype","setting","EAF","nSamples","beta","se","pval")
  
  # return
  data0
}
LDLCData = rbindlist(dumTab3)
LDLCData[,.N,by=setting]
LDLCData[,chr := as.numeric(chr)]

matched = match(LDLCData$rsID,SNPList$rsID)
table(is.na(matched))
table(SNPList[matched,rsID] == LDLCData$rsID)
table(SNPList[matched,EA] == LDLCData$EA)

#' # Run MVMR ####
#' ***
pQTLData[setting == "free",setting := "all"]
setorder(pQTLData,chr,pos_b37,setting)
setorder(LDLCData,chr,pos_b37,setting)
plot(pQTLData$EAF,LDLCData$EAF)

ExposureData[setting == "free",setting := "all"]
ExposureData[grepl("gw",setting),setting2 := "genome-wide"]
ExposureData[grepl("PCSK9",setting),setting2 := "PCSK9_v2"]
ExposureData[grepl("HMGCR",setting),setting2 := "PCSK9_HMGCR"]
ExposureData[is.na(setting2),setting2 := "PCSK9_v1"]
ExposureData = ExposureData[rsID %in% pQTLData$rsID,]

ToDoList = data.table(sex = rep(c("females","all"),each=4),
                      SNPset = rep(c("PCSK9_v1","PCSK9_v2","PCSK9_HMGCR","genome-wide"),2))
OutcomeData[,trait := paste(phenotype,setting,sep = " - ")]

dumTab1 = foreach(i = 1:dim(ToDoList)[1])%do%{
  #i=1
  myRow = ToDoList[i,]
  
  # get sample set (females or all)
  exposure1 = copy(pQTLData)
  exposure1 = exposure1[setting == myRow$sex]
  exposure2 = copy(LDLCData)
  exposure2 = exposure2[setting == myRow$sex]
  
  # get SNP set
  exposure3 = copy(ExposureData)
  if(myRow$SNPset=="PCSK9_HMGCR"){
    exposure3 = exposure3[(setting2 == myRow$SNPset | setting2 == "PCSK9_v2") & grepl(myRow$sex,setting),]
    
  }else{
    exposure3 = exposure3[setting2 == myRow$SNPset & grepl(myRow$sex,setting),]
  }
  exposure1 = exposure1[rsID %in% exposure3$rsID,]
  exposure2 = exposure2[rsID %in% exposure3$rsID,]
  
  # check SNP order and alleles
  stopifnot(exposure1$rsID == exposure3$rsID)
  stopifnot(exposure2$rsID == exposure3$rsID)
  stopifnot(exposure1$EA == exposure2$EA)
  
  # check outcome
  outcome2 = copy(OutcomeData)
  outcome2 = outcome2[rsID %in% exposure3$rsID]
  myOutcomes = unique(outcome2$trait)
  
  dumTab4 = foreach(j = 1:length(myOutcomes))%do%{
    #j=1
    outcome = copy(outcome2)
    outcome = outcome[trait == myOutcomes[j],]
    setorder(outcome,chr,pos_b37)
    stopifnot(outcome$rsID == exposure3$rsID)
    stopifnot(outcome$EA == exposure1$EA)
    
    # running MR
    data_beta = cbind(exposure1$beta,exposure2$beta)
    data_SE = cbind(exposure1$se,exposure2$se)
    
    mvmr_obj = mr_mvinput(bx = as.matrix(data_beta),
                           bxse = as.matrix(data_SE),
                           by = outcome$beta, 
                           byse = outcome$se,
                           snps = outcome$rsID,
                           exposure = c("PCSK9","LDLC"),
                           outcome = myOutcomes[j])
    
    n1 = exposure1[,max(nSamples)]
    n2 = exposure2[,max(nSamples)]
    
    res1 = mr_mvivw(mvmr_obj,nx = c(n1,n2))
    
    res3 = data.table(outcome = res1@Outcome, 
                      sex = myRow$sex,
                      SNPset = myRow$SNPset,
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
    
    res3
    
  }
  
  MVMRTab1 = rbindlist(dumTab4)
  MVMRTab1
  
}
MVMRTab = rbindlist(dumTab1)
MVMRTab[pval_exp1 <0.05]

#' # Save data ####
#' ***
save(MVMRTab, file = "../results/05_MVMR.RData")

stopifnot(pQTLData$rsID == LDLCData$rsID)
stopifnot(pQTLData$setting == LDLCData$setting)

ExposureData_MVMR = cbind(pQTLData[,c(1:5,7,6,8:12)],LDLCData[,c(6,8:12)])
names(ExposureData_MVMR)[7] = "exposure1"
names(ExposureData_MVMR)[8:12] = paste0("exp1_",names(ExposureData_MVMR)[8:12])
names(ExposureData_MVMR)[13] = "exposure2"
names(ExposureData_MVMR)[14:18] = paste0("exp2_",names(ExposureData_MVMR)[14:18])

save(ExposureData_MVMR,file = "../results/05_MVMR_input.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

