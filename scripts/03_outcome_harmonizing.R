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
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
BCAC = fread(GWAS_data_BCAC)

BC_FinnGen_UKB = fread(GWAS_data_finngen_BC)

ParentalLongevity_Death = fread(GWAS_data_GWASCatalog)

CAD = fread(GWAS_data_finngen_CAD)

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

#' Check Longevity (in hg38)
ParentalLongevity_Death
ParentalLongevity_Death = ParentalLongevity_Death[hm_chrom == 1,]
ParentalLongevity_Death = ParentalLongevity_Death[hm_pos %in% IVData$pos_b38,]

#' Aragam CAD sex-stratified meta is given in b37
CAD
setnames(CAD,"#CHR","CHR")
CAD = CAD[CHR == 1,]
CAD = CAD[POS %in% IVData$pos_b38,]
CAD[duplicated(POS)]
table(is.na(CAD$rsid))
CAD = CAD[!is.na(rsid),]

#' # Check alleles ####
#' ***
BC_FinnGen_UKB = BC_FinnGen_UKB[rsid %in% BCAC$rsID,]
CAD = CAD[rsid %in% BCAC$rsID,]

table(BCAC$rsID == BC_FinnGen_UKB$rsid)
table(BCAC$a1 == BC_FinnGen_UKB$ALT)
plot(BCAC$eaf_a1_onco,BC_FinnGen_UKB$FINNGEN_af_alt)
plot(BCAC$eaf_a1_onco,BC_FinnGen_UKB$UKBB_af_alt)

table(BCAC$rsID == CAD$rsid)
table(BCAC$a1 == CAD$ALT)
plot(BCAC$eaf_a1_onco,CAD$FINNGEN_af_alt)
plot(BCAC$eaf_a1_onco,CAD$UKBB_af_alt)

SNPList = copy(IVData)
setorder(SNPList,pos_b38)
SNPList = SNPList[!duplicated(pos_b38)]
SNPList = SNPList[rsID %in% BCAC$rsID,]
table(BCAC$a1 == SNPList$EA)
plot(BCAC$eaf_a1_onco,SNPList$EAF)
plot(BC_FinnGen_UKB$UKBB_af_alt,SNPList$EAF)
plot(CAD$UKBB_af_alt,SNPList$EAF)

matched = match(SNPList$rsID,ParentalLongevity_Death$hm_rsid)
table(is.na(matched))
plot(ParentalLongevity_Death$hm_effect_allele_frequency[matched],SNPList$EAF)
abline(0,1)
abline(1,-1)

#' Okay, looks good! 
#' 
#' # Get outcome in long format ####
#' ***
#' Now I want my outcome data set in a long format. Same columns as for the exposure (plus info for number of cases if possible). I want to split the data into 6 outcomes: 
#' 
#' - BC survival (BCS): meta only
#' - Breast cancer (BC): meta only
#' - Coronary atherosclerosis (CAD): meta only
#' - Parental age at death (PAAD): UKB only
#' 
names(IVData)

#' ## BCAC data
#' I extracted the sample size and breast cancer specific deaths from the abstract of [Morra et al., 2021](https://breast-cancer-research.biomedcentral.com/articles/10.1186/s13058-021-01450-7)
#' 
BCAC[,pos_b38 := SNPList$pos_b38]
names(BCAC)

BCAC[,phenotype := "BCS"]
BCAC[,nSamples := 91686 + 7531]
BCAC[,nCases := 7531]
BCAC[,EAF := (eaf_a1_cogs+eaf_a1_onco)/2]
BCAC[eaf_a1_cogs==0,EAF := eaf_a1_onco]
BCAC[eaf_a1_onco==0,EAF := eaf_a1_cogs]
BCAC = BCAC[,c(18,2,3,19,4,5, 20,23,21,22,14,15,16)]
head(BCAC)
names(BCAC) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                 "phenotype","EAF","nSamples","nCases","beta","se","pval")

#' ## BC data 
#' 
names(IVData)
names(BC_FinnGen_UKB)

#' Sample sizes are taken from the [FinnGen + UKB webpage](https://public-metaresults-fg-ukbb.finngen.fi/) (filtering for "Malignant neoplasm of breast (controls excluding all cancers)")
#' 
BC_FinnGen_UKB[,pos_b37 := SNPList$pos_b37]
names(BC_FinnGen_UKB)

BC_FinnGen_UKB[,phenotype := "BC"]
BC_FinnGen_UKB[,nSamples := 11807 + 205913 + 18786 + 182927]
BC_FinnGen_UKB[,nCases := 11807 + 18786]
BC_FinnGen_UKB[,EAF := (FINNGEN_af_alt+UKBB_af_alt)/2]
BC_FinnGen_UKB[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]
BC_FinnGen_UKB = BC_FinnGen_UKB[,c(22,1,23,2,3,4,  24,27,25,26,17,18,19)]
head(BC_FinnGen_UKB)
names(BC_FinnGen_UKB) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
               "phenotype","EAF","nSamples","nCases","beta","se","pval")

#' ## Parental longevity of parents death data 
#' 
names(IVData)
names(ParentalLongevity_Death)

matched = match(ParentalLongevity_Death$hm_rsid,SNPList$rsID)
ParentalLongevity_Death[, pos_b37 := SNPList$pos_b37[matched] ]
ParentalLongevity_Death[, phenotype := "PAAD"]
ParentalLongevity_Death[, nSamples := 208118 ]

ParentalLongevity_Death = ParentalLongevity_Death[,c(2,3,25,4,5,6,26,11,27,7,20,21)]
names(ParentalLongevity_Death) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                             "phenotype","EAF","nSamples","beta","se","pval")

#' ## CAD data 
#' 
names(IVData)
names(CAD)

#' Sample sizes are taken from the [FinnGen + UKB webpage](https://public-metaresults-fg-ukbb.finngen.fi/) (filtering for "Coronary atherosclerosis")
#' 
CAD[,pos_b37 := SNPList$pos_b37]
names(CAD)

CAD[,phenotype := "CAD"]
CAD[,nSamples := 51589 + 31198 + 343079 + 382052]
CAD[,nCases := 51589 + 31198]
CAD[,EAF := (FINNGEN_af_alt+UKBB_af_alt)/2]
CAD[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]
CAD = CAD[,c(22,1,23,2,3,4,  24,27,25,26,17,18,19)]
head(CAD)
names(CAD) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
               "phenotype","EAF","nSamples","nCases","beta","se","pval")


#' Merge data sets info one data table
outcomeData = rbind(BCAC, BC_FinnGen_UKB, ParentalLongevity_Death, CAD,fill=T)
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

