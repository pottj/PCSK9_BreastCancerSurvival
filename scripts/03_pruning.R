#' ---
#' title: "Pruning"
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
#' For the pQTLs, I do not need to prune, because I will reuse published instruments. 
#' 
#' For the other exposures, I will check here: 
#' 
#' - Gene expression: LD pruning, and use independent SNPs only (LD r^2<0.1)
#' - LDL-C: position pruning, and eliminate SNPs within 1Mb around lead variant
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

#' # Load data ####
#' ***
load("../results/02_Exposure_for_MR.RData")
ExposureData

GTEx = copy(ExposureData)
GTEx = GTEx[phenotype == "PCSK9 gene expression" & !grepl("ratio",setting),]

LDLC = copy(ExposureData)
LDLC = LDLC[phenotype == "LDL-C levels"  & !grepl("ratio",setting), ]

#' # LD pruning for GTEx ####
#' ***
#' I will use FUMA for this. Here, I will first create the FUMA input per tissue, upload the data to FUMA, and then load the results (lead SNPs) and filter for them. 
#' 
GTEx[,.N,setting]

myTissues = GTEx[,unique(setting)]
myTissues2 = gsub("[()]","",myTissues)
myTissues2 = gsub(" - ","_",myTissues2)
myTissues2 = gsub(" ","_",myTissues2)
myTissues2

for(i in 1:length(myTissues2)){
  #i=1
  data = copy(GTEx)
  data = data[setting == myTissues[i],]

  filename = paste0("../temp/FUMA_input/SNPList_",myTissues2[i],".txt")
  write.table(data[,c(1,2,3,5,6,9:13)],file = filename,
              col.names = T,row.names = F,quote = F)
}

#' Download FUMA results and check lead SNPs (LD r2<0.1 in EUR)
#' 
myFiles = list.files("../temp/FUMA_results/",pattern = "lead")

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  leadSNPs = fread(paste0("../temp/FUMA_results/",myFiles[i]))
  
  myTissue = gsub("leadSNPs_job6218.._","",myFiles[i])
  myTissue = gsub(".txt","",myTissue)
  filt = grep(myTissue,myTissues2)
  
  data = copy(GTEx)
  data = data[setting == myTissues[filt],]
  data = data[rsID %in% leadSNPs$rsID,]
  
  data
}
GTEx_pruned = rbindlist(dumTab1)
GTEx[,.N,setting]
GTEx_pruned[,.N,setting]

#' # Position pruning for LDLC ####
#' ***
mySettings = unique(LDLC$setting)
mySNPs_PCSK9 = fread("../temp/SNPList_PCSK9_Yang2024.txt",header = F)
mySNPs_HMGCR = fread("../temp/SNPList_HMGCR_Yang2024.txt",header = F)

dumTab2 = foreach(i = 1:length(mySettings))%do%{
  #i=1
  data = copy(LDLC)
  data = data[setting == mySettings[i],]
  
  # get PCSK9 and HMGCR specific sets
  data2 = copy(data)
  data2 = data2[rsID %in% mySNPs_PCSK9$V1,]
  data2[,setting := paste(setting, "PCSK9 SNPs")]
  data3 = copy(data)
  data3 = data3[rsID %in% mySNPs_HMGCR$V1,]
  data3[,setting := paste(setting, "HMGCR SNPs")]
  
  # do pruning
  myChrs = unique(as.numeric(data$chr))
  
  dumTab4 = foreach(s2 = myChrs) %do% {
    # s2 = myChrs[1]
    subdata2 = copy(data)
    subdata2 = subdata2[chr == s2, ]
    
    setkey(subdata2, pos_b37)
    
    if(dim(subdata2)[1]<=1){
      subdata2[, keep := T]
      subdata2[, NR_SNPs := 0]
    }else{
      subdata2[, keep := NA]
      subdata2[, NR_SNPs := as.numeric(NA)]
      
      smallestDist = getSmallestDist(subdata2[, pos_b37])
      while(smallestDist < 500000) {
        maxZ = max(subdata2[is.na(keep), absZScore])
        mybp = subdata2[absZScore == maxZ & is.na(keep), pos_b37]
        if(length(mybp)>1){
          mybp = mybp[1]
        }
        subdata2[pos_b37 == mybp, keep := T]
        
        #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
        myFilt = (subdata2[, pos_b37] < (mybp - 500000)) | 
          (subdata2[, pos_b37] > (mybp + 500000)) | 
          subdata2[, keep] 
        myFilt[is.na(myFilt)] = FALSE
        subdata2 = subdata2[myFilt == TRUE, ]
        
        subdata2[pos_b37 == mybp, NR_SNPs := sum(myFilt==F)]
        smallestDist = getSmallestDist(subdata2[, pos_b37])
      }
      
      #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
      subdata2[is.na(keep), NR_SNPs := 0]
      subdata2[is.na(keep), keep := TRUE]
    }
    
    subdata2
  }
  data1 = rbindlist(dumTab4)
  data1 = data1[NR_SNPs>9,]
  data1[,setting := paste(setting, "gw SNPs")]
  
  # return
  data1 = data1[,c(1:14)]
  data5 = rbind(data1,data2,data3)
  data5
}
LDLC_pruned = rbindlist(dumTab2)
LDLC[,.N,setting]
LDLC_pruned[,.N,setting]

#' # Save data ####
#' ***
ExposureData = ExposureData[!(phenotype == "PCSK9 gene expression" & !grepl("ratio",setting)),]
ExposureData = ExposureData[!(phenotype == "LDL-C levels" & !grepl("ratio",setting)),]

ExposureData = rbind(ExposureData,GTEx_pruned,LDLC_pruned)
ExposureData[,.N,by = c("phenotype","setting")]

save(ExposureData,file = "../results/03_Exposure_for_MR_pruned.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

