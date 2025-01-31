#' ---
#' title: "Check outcome data"
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
#' In this script, I load the outcome data sets and check for overlapping SNPs (maybe not all potential IVs available). Then I will harmonize the matching SNPs (same effect allele, check EAF).
#'  
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
BCAC = fread(GWAS_data_BCAC)

BC_FinnGen_UKB = fread(GWAS_data_finngen)

load("../results/02_potentialIVs.RData")

#' # Filter data ####
#' ***
#' BCAC is given in hg19 = b37
BCAC
BCAC = BCAC[chromosome == 1,]
BCAC = BCAC[position %in% IVData$pos_b37,]

#' FinnGen & UKB meta is given in b38
BC_FinnGen_UKB
setnames(BC_FinnGen_UKB,"#CHR","CHR")
BC_FinnGen_UKB = BC_FinnGen_UKB[CHR == 1,]
BC_FinnGen_UKB = BC_FinnGen_UKB[POS %in% IVData$pos_b38,]
BC_FinnGen_UKB[duplicated(POS)]
table(is.na(BC_FinnGen_UKB$rsid))
BC_FinnGen_UKB = BC_FinnGen_UKB[!is.na(rsid),]

#' check rsID in FinnGen & UKB meta
table(BC_FinnGen_UKB$rsid %in% IVData$rsID)
matched = match(BC_FinnGen_UKB$POS,IVData$pos_b38)
table(BC_FinnGen_UKB$rsid == IVData[matched,rsID])
BC_FinnGen_UKB[rsid != IVData[matched,rsID]]
IVData[pos_b38 %in% BC_FinnGen_UKB[rsid != IVData[matched,rsID],POS],]

#' Look-up in dbSNP: 
#' 
#' - rs115749208 (IVData) has merged into rs111835689 (FinnGen & UKB)
#' - rs149790879 has merged into rs3118722
#' - rs113675231 has merged into rs3118721 (rs145663020 not found?!)
#' - rs201732208 has merged into rs4500361
#' - rs370423508 has merged into rs59267263
#' 
#' Okay, I want to use these update rsIDs. 
#' 
matched = match(IVData$pos_b38,BC_FinnGen_UKB$POS)
IVData[,rsID := BC_FinnGen_UKB[matched,rsid]]

#' Add updated rsID to BCAC
#' 
matched = match(BCAC$position,IVData$pos_b37)
BCAC[,rsID := IVData[matched,rsID]]

#' # Check alleles ####
#' ***
table(BCAC$rsID == BC_FinnGen_UKB$rsid)
table(BCAC$a1 == BC_FinnGen_UKB$ALT)
plot(BCAC$eaf_a1_onco,BC_FinnGen_UKB$FINNGEN_af_alt)
plot(BCAC$eaf_a1_onco,BC_FinnGen_UKB$UKBB_af_alt)

SNPList = copy(IVData)
setorder(SNPList,pos_b38)
SNPList = SNPList[!duplicated(pos_b38)]
table(BCAC$a1 == SNPList$EA)
plot(BCAC$eaf_a1_onco,SNPList$EAF)
plot(BC_FinnGen_UKB$UKBB_af_alt,SNPList$EAF)

#' Okay, looks good! 
#' 
#' # Get outcome in long format ####
#' ***
#' Now I want my outcome data set in a long format. Same columns as for the exposure (plus info for number of cases if possible). I want to split the data into 6 outcomes: 
#' 
#' - BC survival (BCS): onco, cogs, meta
#' - BC prevalence (BCP): FinnGen, UKB, meta
#' 
names(IVData)

#' ## BCAC data
#' I do not know the sample size per array, but I found the sample size for the main analysis in the Supplemental Table S3 of [Morra et al., 2021](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-021-01450-7#availability-of-data-and-materials)
#' 
BCAC[,pos_b38 := SNPList$pos_b38]
names(BCAC)

BCAC1 = copy(BCAC)
BCAC1[,phenotype := "BCS_cogs"]
BCAC1 = BCAC1[,c(18,2,3,19,4,5, 20,6,8,9,10)]
head(BCAC1)
names(BCAC1) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                 "phenotype","EAF","beta","se","pval")

BCAC2 = copy(BCAC)
BCAC2[,phenotype := "BCS_onco"]
BCAC2 = BCAC2[,c(18,2,3,19,4,5, 20,7,11,12,13)]
head(BCAC2)
names(BCAC2) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                 "phenotype","EAF","beta","se","pval")

BCAC3 = copy(BCAC)
BCAC3[,phenotype := "BCS_meta"]
BCAC3[,nSamples := 99217]
BCAC3[,EAF := (eaf_a1_cogs+eaf_a1_onco)/2]
BCAC3[eaf_a1_cogs==0,EAF := eaf_a1_onco]
BCAC3[eaf_a1_onco==0,EAF := eaf_a1_cogs]
BCAC3 = BCAC3[,c(18,2,3,19,4,5, 20,22,21,14,15,16)]
head(BCAC3)
names(BCAC3) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                 "phenotype","EAF","nSamples","beta","se","pval")

names(IVData)
names(BC_FinnGen_UKB)

#' ## BC data 
#' 
#' Sample sizes are taken from the [FinnGen + UKB webpage](https://public-metaresults-fg-ukbb.finngen.fi/) (filtering for "Malignant neoplasm of breast (controls excluding all cancers)")
#' 
BC_FinnGen_UKB[,pos_b37 := SNPList$pos_b37]
names(BC_FinnGen_UKB)

BC1 = copy(BC_FinnGen_UKB)
BC1[,phenotype := "BCP_FinnGen"]
BC1[,nSamples := 18786 + 182927]
BC1[,nCases := 18786]
BC1 = BC1[,c(22,1,23,2,3,4,  24,9,25,26,6,7,8)]
head(BC1)
names(BC1) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
               "phenotype","EAF","nSamples","nCases","beta","se","pval")

BC2 = copy(BC_FinnGen_UKB)
BC2[,phenotype := "BCP_UKB"]
BC2[,nSamples := 11807 + 205913]
BC2[,nCases := 11807]
BC2 = BC2[,c(22,1,23,2,3,4,  24,15,25,26,12,13,14)]
head(BC2)
names(BC2) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
               "phenotype","EAF","nSamples","nCases","beta","se","pval")

BC3 = copy(BC_FinnGen_UKB)
BC3[,phenotype := "BCP_meta"]
BC3[,nSamples := 11807 + 205913 + 18786 + 182927]
BC3[,nCases := 11807 + 18786]
BC3[,EAF := (FINNGEN_af_alt+UKBB_af_alt)/2]
BC3[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]
BC3 = BC3[,c(22,1,23,2,3,4,  24,27,25,26,17,18,19)]
head(BC3)
names(BC3) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
               "phenotype","EAF","nSamples","nCases","beta","se","pval")

#' Merge data sets info one data table
outcomeData = rbind(BC3, BC1, BC2, BCAC3,BCAC1,BCAC2,fill=T)
outcomeData[is.na(beta)]
outcomeData[beta==0,]

outcomeData = outcomeData[beta != 0,]
outcomeData = outcomeData[!is.na(beta),]

#' # Save data ####
#' ***
save(IVData,file = "../results/03_IVs_updated.RData")
save(outcomeData,file = "../results/03_outcome_harmonized.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

