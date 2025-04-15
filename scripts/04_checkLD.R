#' ---
#' title: "Check IVs LD"
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
#' Here I will extract the rsID per tissue or subgroup and check their pairwise LD in [LDlink](https://ldlink.nih.gov/?tab=ldmatrix). 
#' 
#' But first I will exclude rs562556, as the missense variant will only be used in the MR-ratio approach. I will also exclude testis tissue, as this is not relevant in the female-specific breast cancer situation.  
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
load("../results/03_IVs_updated.RData")
IVData = IVData[rsID != "rs562556",]
IVData = IVData[phenotype != "Testis",]
myTraits = IVData[,unique(phenotype)]
myTraits

#' # Check NR of SNPs ####
#' ***
IVData[,.N,by=phenotype]

#' **Summary**: 
#' 
#' - For liver, pancreas, and the two skin tissues I will use only 1 instrument. 
#' - For protein levels, I will use the four SNPs as described in my previous paper
#' - for the remaining tissues (adipose, brain, lung, nerve, and whole blood), I will extract the SNP effects and upload to FUMA to get all independent SNPs
#' 
IVData[phenotype %in% myTraits[c(4,7:9,11,12)],flag := T]
IVData[rsID == "rs17111503" & phenotype==myTraits[8],flag := F]

myTraits = IVData[is.na(flag),unique(phenotype)]
myTraits

#' # LD check ####
#' ***
#' Create files for FUMA upload. In addition, get SNP list for LD link plot of all SNPs

for(i in 1:length(myTraits)){
  #i=1
  data = copy(IVData)
  data = data[phenotype == myTraits[i],]
  setorder(data,pval)
  message("Working on trait ",i,": ",myTraits[i])
  print(data$rsID,quote = F)
  
  filename = paste0("../temp/04_SNPList_",myTraits[i],".txt")
  write.table(data,file = filename,col.names = T,row.names = F,quote = F)
  
}

#' Download FUMA results and check lead SNPs (LD r2<0.1 in EUR)
#' 
myFiles = list.files("../temp/FUMA_jobs/",pattern = "lead")

for(i in 1:length(myFiles)){
  #i=1
  leadSNPs = fread(paste0("../temp/FUMA_jobs/",myFiles[i]))
  myTrait = gsub("leadSNPs_","",myFiles[i])
  myTrait = gsub(".txt","",myTrait)
  IVData[phenotype==myTrait,flag:=F]
  IVData[phenotype==myTrait & rsID %in% leadSNPs$rsID,flag:=T]
  
  message("Working on trait ",i,": ",myTraits[i]," - pruned:")
  print(IVData[flag==T & phenotype==myTrait, rsID],quote = F)
  
}

#' # Save data ####
#' ***
IVData = IVData[flag==T,]
save(IVData,file = "../results/04_IVs_pruned.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

