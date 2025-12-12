#' ---
#' title: "MR: PCSK9 on outcomes"
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
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/03_Exposure_for_MR_pruned.RData")
load("../results/02_Outcome_for_MR.RData")
OutcomeData = OutcomeData[rsID %in% ExposureData$rsID,]
OutcomeData[,chr := as.numeric(chr)]

#' # Load FinnGen data ####
#' ***
FinnGen = fread(FinnGen_BCS_add)
head(FinnGen)
FinnGen_rec = fread(FinnGen_BCS_rec)
head(FinnGen_rec)

BCAC = copy(OutcomeData)
BCAC = BCAC[source == "Morra et al.",]
matched = match(FinnGen$rsID,BCAC$rsID)
stopifnot(BCAC$rsID[matched] == FinnGen$rsID)
BCAC = BCAC[matched,]

stopifnot(BCAC$EA == FinnGen$EA)
plot(BCAC$EAF,FinnGen$EAF)
plot(BCAC$beta,FinnGen$logHR)
cor.test(BCAC$beta,FinnGen$logHR)
plot(FinnGen$logHR)

filt = FinnGen$pval<0.05 | BCAC$pval<0.05
plot(BCAC$beta[filt],FinnGen$logHR[filt])
abline(0,1)
cor.test(BCAC$beta[filt],FinnGen$logHR[filt])

#' Okay, BCAC and FinnGen are not correlated - maybe due to sample selection (that is what I wanted to test here).
#' Other thing: there are some really big effects for 5 SNPs (logHR>0.5) with low MAF - I will consider to remove them, because I think they are not plausible. 
#' 
#' Otherwise, I just try to make the data set as similar to BCAC as possible. 
names(BCAC)
names(FinnGen)
names(FinnGen_rec)

FinnGen[,phenotype := "BCS"]
FinnGen[,setting := "females"]
FinnGen[,source := "FinnGen"]
FinnGen[,geneticModel := "additive"]
FinnGen[,nSamples := 4648]
FinnGen[,nCases := 288]
setnames(FinnGen,"logHR","beta")
FinnGen[,outcomeID := paste(phenotype,setting,source,geneticModel,sep=" - ")]
table(names(FinnGen) %in% names(BCAC))
setcolorder(FinnGen,names(BCAC))

FinnGen_rec[,phenotype := "BCS"]
FinnGen_rec[,setting := "females"]
FinnGen_rec[,source := "FinnGen"]
FinnGen_rec[,geneticModel := "recessive"]
FinnGen_rec[,nSamples := 4648]
FinnGen_rec[,nCases := 288]
setnames(FinnGen_rec,"logHR","beta")
FinnGen_rec[,outcomeID := paste(phenotype,setting,source,geneticModel,sep=" - ")]
table(names(FinnGen_rec) %in% names(BCAC))
setcolorder(FinnGen_rec,names(BCAC))

filt = FinnGen$beta > 0.5
table(filt)
FinnGen = FinnGen[!filt,]

#' # Merga and save ####
#' ***
#' Merge with outcome data 
OutcomeData = rbind(OutcomeData,FinnGen,FinnGen_rec)

#' Reduce to overlapping SNPs (I want the same SNP set in all my analyses!)
ExposureData[!is.element(rsID,FinnGen$rsID), ]
ExposureData = ExposureData[is.element(rsID,FinnGen$rsID),]
OutcomeData = OutcomeData[is.element(rsID,ExposureData$rsID),]

save(ExposureData,file = "../results/03b_Exposure_for_MR_filtered.RData")
save(OutcomeData,file = "../results/03b_Outcome_for_MR_filtered.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
