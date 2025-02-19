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
load("../results/03_IVs_updated.RData")
myTraits = IVData[,unique(phenotype)]
myTraits

#' # LD check ####
#' ***
#' 1) get rsIDs per phenotype
#' 2) check pairwise LD in [LDlink](https://ldlink.nih.gov/?tab=ldmatrix)
#' 3) pick best associated & independent SNPs stepwise
#' 
IVData[,flag := F]

#' ## Tissue 1: Brain cerebellar hemisphere
#' 
data = copy(IVData)
data = data[phenotype == myTraits[1],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are two cluster, and the two signals are rs11583680 and rs77947184 (LD r2 = 0.494)
#' 
IVData[phenotype == myTraits[1] & rsID %in% c("rs11583680","rs77947184"),flag := T] 
IVData[flag == T,]

dumTab1 = data.table(SNP1 = c("rs77947184"),
                    OA1 = c("C"),
                    EA1 = c("T"),
                    SNP2 = c("rs11583680"), 
                    OA2 = c("C"),
                    EA2 = c("T"),
                    p00 = c(0.789),
                    p0x = c(0.866),
                    p1x = c(0.134),
                    px0 = c(0.801), 
                    px1 = c(0.199))
dumTab1

#' ## Tissue 2: Brain cerebellum
#' 
data = copy(IVData)
data = data[phenotype == myTraits[2],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are two cluster, and the two signals are rs11583680 and rs7546522 (LD r2 = 0.566)
#' 
IVData[phenotype == myTraits[2] & rsID %in% c("rs11583680","rs7546522"),flag := T] 
IVData[flag == T,]

dumTab2 = data.table(SNP1 = c("rs11583680"),
                     OA1 = c("C"),
                     EA1 = c("T"),
                     SNP2 = c("rs7546522"), 
                     OA2 = c("C"),
                     EA2 = c("T"),
                     p00 = c(0.840),
                     p0x = c(0.866),
                     p1x = c(0.134),
                     px0 = c(0.871), 
                     px1 = c(0.129))
dumTab2

#' ## Tissue 3: Liver
#' 
data = copy(IVData)
data = data[phenotype == myTraits[3],]
print(data$rsID,quote = F)

#' There is only one SNP: rs553741
#' 
IVData[phenotype == myTraits[3] & rsID %in% c("rs553741"),flag := T] 
IVData[flag == T,]

#' ## Tissue 4: Lung
#' 
data = copy(IVData)
data = data[phenotype == myTraits[4],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are 6 cluster, and the 6 signals are rs4500361, rs2479395, rs11206510, rs3976734, rs2495502, and rs2479396 (pairwise LD r2 < 0.6)
#' 
IVData[phenotype == myTraits[4] & rsID %in% c("rs4500361","rs2479395","rs11206510",
                                              "rs3976734", "rs2495502", "rs2479396"),flag := T] 
IVData[flag == T,]

dumTab4 = data.table(SNP1 = c(rep("rs2479396",5),rep("rs2479395",4),
                             rep("rs2495502",3),rep("rs3976734",2),"rs4500361"),
                    OA1 = c(rep("A",5), rep("C",4),rep("C",3), rep("A",2),"C"),
                    EA1 = c(rep("G",5), rep("T",4),rep("G",3), rep("G",2),"T"),
                    SNP2 = c("rs2479395","rs2495502","rs3976734","rs4500361","rs11206510",
                             "rs2495502","rs3976734","rs4500361","rs11206510",
                             "rs3976734","rs4500361","rs11206510",
                             "rs4500361","rs11206510",
                             "rs11206510"), 
                    OA2 = c("C","C","A","C","C",    "C","A","C","C",   "A","C","C",  "C","C","C"),
                    EA2 = c("T","G","G","T","T",    "G","G","T","T",   "G","T","T",  "T","T","T"),
                    p00 = c(0.597, 0.232, 0.480, 0.098, 0.048,
                            0.192, 0.500, 0.037, 0.038, 
                            0.291, 0.238, 0.163,
                            0.038, 0.007,
                            0.166),
                    p0x = c(rep(0.672,5),rep(0.668,4),rep(0.518,3),rep(0.645,2),0.243),
                    p1x = c(rep(0.328,5),rep(0.332,4),rep(0.482,3),rep(0.355,2),0.757),
                    px0 = c(0.668, 0.518, 0.645, 0.243, 0.172,
                            0.518, 0.645, 0.243, 0.172,
                            0.645, 0.243, 0.172,
                            0.243, 0.172, 
                            0.172), 
                    px1 = c(0.332, 0.482, 0.355, 0.757, 0.828,
                            0.482, 0.355, 0.757, 0.828,
                            0.355, 0.757, 0.828,
                            0.757, 0.828, 
                            0.828))
dumTab4

#' ## Tissue 5: Nerve
#' 
data = copy(IVData)
data = data[phenotype == myTraits[5],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are 7 cluster, and the 7 signals are rs4500361, rs11206510, rs2479395, rs3976734, rs12046679, rs11488228, and rs2495502 (pairwise LD r2 < 0.6)
#' 
IVData[phenotype == myTraits[5] & rsID %in% c("rs4500361","rs11206510","rs2479395","rs3976734",
                                              "rs12046679","rs11488228","rs2495502"),flag := T] 
IVData[flag == T,]

dumTab5 = data.table(SNP1 = c(rep("rs2479395",6),rep("rs2495502",5),
                              rep("rs3976734",4),rep("rs12046679",3),
                              rep("rs4500361",2),"rs11488228"),
                     OA1 = c(rep("C",6),rep("C",5), rep("A",4),
                             rep("A",3), rep("C",2),"A"),
                     EA1 = c(rep("T",6),rep("G",5), rep("G",4),
                             rep("G",3), rep("T",2),"C"),
                     SNP2 = c("rs2495502","rs3976734","rs12046679","rs4500361","rs11488228","rs11206510",
                              "rs3976734","rs12046679","rs4500361","rs11488228","rs11206510",      
                              "rs12046679","rs4500361","rs11488228","rs11206510",      
                              "rs4500361","rs11488228","rs11206510",
                              "rs11488228","rs11206510",
                              "rs11206510"), 
                     OA2 = c("C","A","A","C","A","C", 
                             "A","A","C","A","C",
                             "A","C","A","C",
                             "C","A","C",
                             "A","C",
                             "C"),
                     EA2 = c("G","G","G","T","C","T",
                             "G","G","T","C","T",
                             "G","T","C","T",
                             "T","C","T",
                             "C","T",
                             "T"),
                     p00 = c(0.192, 0.500, 0.501, 0.037, 0.120, 0.038,
                             0.291, 0.241, 0.238, 0.119, 0.163,
                             0.386, 0.038, 0.119, 0.007,
                             0.242, 0.001, 0.172,
                             0.001, 0.165, 
                             0),
                     p0x = c(rep(0.668,6),rep(0.518,5),rep(0.645,4),rep(0.716,3),rep(0.243,2),0.122),
                     p1x = c(rep(0.332,6),rep(0.482,5),rep(0.355,4),rep(0.284,3),rep(0.757,2),0.878),
                     px0 = c(0.518, 0.645, 0.716, 0.243, 0.122, 0.172,
                             0.645, 0.716, 0.243, 0.122, 0.172,
                             0.716, 0.243, 0.122, 0.172,
                             0.243, 0.122, 0.172,
                             0.122, 0.172,
                             0.172), 
                     px1 = c(0.482, 0.355, 0.284, 0.757, 0.878, 0.828,
                             0.355, 0.284, 0.757, 0.878, 0.828,
                             0.284, 0.757, 0.878, 0.828,
                             0.757, 0.878, 0.828,
                             0.878, 0.828, 
                             0.828))
dumTab5

#' ## Tissue 6: Pancreas
#' 
data = copy(IVData)
data = data[phenotype == myTraits[6],]
print(data$rsID,quote = F)

#' There is only one SNP: rs535471
#' 
IVData[phenotype == myTraits[6] & rsID %in% c("rs535471"),flag := T] 
IVData[flag == T,]

#' ## Setting 1: females
#' 
data = copy(IVData)
data = data[phenotype == myTraits[7],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are 4 cluster, and the 4 best signals are rs11591147, rs28385704, rs472495, and rs2495491 (pairwise LD r2 < 0.6). But I want to use the same SNPs throughout for my pQTLs. Hence I stick to rs11591147, rs693668, rs11583680, and rs2495491.
#' 
IVData[phenotype == myTraits[7] & rsID %in% c("rs11591147","rs693668","rs11583680","rs2495491"),flag := T] 
IVData[flag == T,]

dumTab7 = data.table(SNP1 = c(rep("rs2495491",3),rep("rs11591147",2),"rs11583680"),
                     OA1 = c(rep("G",3), rep("G",2),"C"),
                     EA1 = c(rep("T",3), rep("T",2),"T"),
                     SNP2 = c("rs11591147","rs11583680","rs693668",
                              "rs11583680","rs693668",
                              "rs693668"), 
                     OA2 = c("G","C","A",  "C","A", "A"),
                     EA2 = c("T","T","G",  "T","G", "G"),
                     p00 = c(0.247, 0.243, 0.189,
                             0.845, 0.614,
                             0.572),
                     p0x = c(rep(0.247,3),rep(0.979,2),0.866),
                     p1x = c(rep(0.753,3),rep(0.021,2),0.134),
                     px0 = c(0.979, 0.866, 0.625,
                             0.866, 0.625, 
                             0.625), 
                     px1 = c(0.021, 0.134, 0.375,
                             0.134, 0.375, 
                             0.375))
dumTab7

#' ## Setting 2: females statin-free
#' 
data = copy(IVData)
data = data[phenotype == myTraits[8],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are 4 cluster, and the 4 best signals are rs11591147, rs28385704, rs693668, and rs2495491 (pairwise LD r2 < 0.6). But I want to use the same SNPs throughout for my pQTLs. Hence I stick to rs11591147, rs693668, rs11583680, and rs2495491.
#' 
IVData[phenotype == myTraits[8] & rsID %in% c("rs11591147","rs693668","rs11583680","rs2495491"),flag := T] 
IVData[flag == T,]

dumTab8 = data.table(SNP1 = c(rep("rs2495491",3),rep("rs11591147",2),"rs11583680"),
                     OA1 = c(rep("G",3), rep("G",2),"C"),
                     EA1 = c(rep("T",3), rep("T",2),"T"),
                     SNP2 = c("rs11591147","rs11583680","rs693668",
                              "rs11583680","rs693668",
                              "rs693668"), 
                     OA2 = c("G","C","A",  "C","A", "A"),
                     EA2 = c("T","T","G",  "T","G", "G"),
                     p00 = c(0.247, 0.243, 0.189,
                             0.845, 0.614,
                             0.572),
                     p0x = c(rep(0.247,3),rep(0.979,2),0.866),
                     p1x = c(rep(0.753,3),rep(0.021,2),0.134),
                     px0 = c(0.979, 0.866, 0.625,
                             0.866, 0.625, 
                             0.625), 
                     px1 = c(0.021, 0.134, 0.375,
                             0.134, 0.375, 
                             0.375))
dumTab8

#' ## Setting 3: statin-free (sex-combined)
#' 
data = copy(IVData)
data = data[phenotype == myTraits[9],]
print(data$rsID,quote = F)
setorder(data,pval)
print(data$rsID,quote = F)

#' There are 4 cluster, and the 4 best signals are rs11591147, rs693668, rs11583680, and rs2495491 (pairwise LD r2 < 0.6)
#' 
IVData[phenotype == myTraits[9] & rsID %in% c("rs11591147","rs693668","rs11583680","rs2495491"),flag := T] 
IVData[flag == T,]

dumTab9 = data.table(SNP1 = c(rep("rs2495491",3),rep("rs11591147",2),"rs11583680"),
                     OA1 = c(rep("G",3), rep("G",2),"C"),
                     EA1 = c(rep("T",3), rep("T",2),"T"),
                     SNP2 = c("rs11591147","rs11583680","rs693668",
                              "rs11583680","rs693668",
                              "rs693668"), 
                     OA2 = c("G","C","A",  "C","A", "A"),
                     EA2 = c("T","T","G",  "T","G", "G"),
                     p00 = c(0.247, 0.243, 0.189,
                             0.845, 0.614,
                             0.572),
                     p0x = c(rep(0.247,3),rep(0.979,2),0.866),
                     p1x = c(rep(0.753,3),rep(0.021,2),0.134),
                     px0 = c(0.979, 0.866, 0.625,
                             0.866, 0.625, 
                             0.625), 
                     px1 = c(0.021, 0.134, 0.375,
                             0.134, 0.375, 
                             0.375))
dumTab9

#' # Save data ####
#' ***
IVData = IVData[flag==T,]
save(IVData,file = "../results/04_IVs_pruned.RData")

dumTab1[,phenotype := myTraits[1]]
dumTab2[,phenotype := myTraits[2]]
dumTab4[,phenotype := myTraits[4]]
dumTab5[,phenotype := myTraits[5]]
dumTab7[,phenotype := myTraits[7]]
dumTab8[,phenotype := myTraits[8]]
dumTab9[,phenotype := myTraits[9]]

dumTab = rbind(dumTab1,dumTab2,dumTab4,dumTab5,dumTab7, dumTab8, dumTab9)
dumTab[,D := p00 - p0x*px0]                    
dumTab[,r := D / sqrt(p0x*p1x*px0*px1)]
dumTab[,r2 := r^2]

save(dumTab,file = "../temp/04_IV_correlation.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

