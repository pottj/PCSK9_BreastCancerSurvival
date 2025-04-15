#' ---
#' title: "Check PCSK9 meta-GWAS data"
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
#' In this script, I load my PCSK9 summary statistics as generated in the [sex- and statin-stratified meta-GWAS project](https://github.com/GenStatLeipzig/GWAMA_PCSK9_strat) (script *01_createDataFiles.R*) 
#' 
#' I will focus on females, statin-free females (to check if statins having an effect), and the sex-combined statin-free set (to increase power as the sex-effect might change the estimate, but not the effect direction).
#'  
#' I will first load each data set, filter for PCSK9 gene region, and get the number of genome-wide significant SNPs (potential IVs). Finally, I save the significant pQTLs for later use. 
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

#' Okay, for both data sets I have some genome-wide significant hits. I will use the independent signals as described in my previous paper. 
#' 
#' # Harmonize data ####
#' ***
#' I want to generate a harmonized format for all significant SNPs (either p- or eQTL). For each SNP I want: 
#' 
#' - rsID, bp_hg19, bp_hg38, EA, OA
#' 
#' For each tissue (eQTLs) or setting (pQTLs) I want
#' 
#' - EAF, sample size, beta, SE, pval 
#' 
#' ## Harmonize eQTLs data set 
#'  
load("../temp/GTExV8_potentialIVs.RData")
table(is.element(GTExData$pos_b37, pQTLData$bp_hg19))

#' Okay, only 3 SNP cannot be found. I will look-up its rsID manually in the [GTEx data portal](https://gtexportal.org/home/) (using *variant_id_b38*). I will also do this for any variant without rsID in the pQTL data (will have rsID == 1)
#' 
matched = match(GTExData$pos_b37,pQTLData$bp_hg19)
GTExData[,rsID := pQTLData[matched,rsID]]
GTExData[rsID == 1,]
GTExData[rsID == 1, rsID := "rs536607133"]
GTExData[is.na(rsID),]
GTExData[is.na(rsID), rsID := c("rs41294823","rs35595930","rs35595930")]

#' Now I add the effect allele as reported in my pQTL data. I only do this so I can generate the EAF instead of the MAF. 
#' 
GTExData[,EA_pQTLs := pQTLData[matched,EA]]
GTExData[,EAF_pQTLs := pQTLData[matched,EAF]]
filt = GTExData$effect_allele == GTExData$EA_pQTLs
table(filt)
plot(GTExData$EAF_pQTLs,GTExData$maf)
abline(0,1)
abline(1,-1)
GTExData[,eaf := maf]
GTExData[EAF_pQTLs>0.5,eaf := 1-maf]
plot(GTExData$EAF_pQTLs,GTExData$eaf)
abline(0,1)
abline(1,-1)

GTExData[rsID == "rs41294823"]

#' According to dbSNP the variant rs41294823 has alleles G/T and the T allele is the minor allele. No need to add changes here. 

#' Now I can filter for the relevant columns
#' 
names(GTExData)
eQTLData = GTExData[, c(24,10,17,11,12,13, 22,27,14,7,8,6)]
head(eQTLData)
names(eQTLData)

#' ## Harmonize pQTL data set
#' 
load("../temp/GTExV8_selectedTissues_PCSK9.RData")
pQTLData = pQTLData[pval<1e-6,]

#' From my meta-GWAS: these were the independent variants according to GCTA COJO slct (*03_GCTA_COJO.R*) or the missense mutation reported by Mei et al.
#' 
pQTLData = pQTLData[rsID %in% c("rs2495491","rs11591147","rs11583680","rs693668","rs562556")]
matched = match(pQTLData$bp_hg19,GTExData$pos_b37)
table(is.na(matched))

pQTLData[,pos_b38 := GTExData[matched,pos_b38]]
pQTLData[,table(invalidAssocs)]
pQTLData[invalidAssocs==T,]

#' There are 2 invalid associations in the sex-combined setting. This heterogeneity is caused by the sex-biased SNP effect (effect in men stronger than in women). I only use this phenotype to increase power, so the heterogeneity might not be a real problem, since it is still a true positive. 

names(pQTLData)
pQTLData = pQTLData[, c(17,2,3,18,5,4, 16,6,8,10,11,12)]
head(pQTLData)
pQTLData[,max(pval,na.rm = T),by=phenotype]
pQTLData[,min(pval,na.rm = T),by=phenotype]

#' # Save data ####
#' ***
names(eQTLData) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                    "phenotype","EAF","nSamples","beta","se","pval")
names(pQTLData) = names(eQTLData)

IVData = rbind(eQTLData,pQTLData)
IVData[,length(unique(rsID))]
save(IVData,file = "../results/02_potentialIVs.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

