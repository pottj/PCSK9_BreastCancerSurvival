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

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../temp/PCSK9_pQTLs.RData")

pQTLData[ rsID == "rs11583680",]
pQTLData = pQTLData[rsID == "rs11583680",]
pQTLData

load("../temp/GTExV8_selectedTissues_PCSK9.RData")
GTExData = GTExData[pos_b37 %in% pQTLData$bp_hg19,]
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
save(IVData,file = "../results/08_SumStats_rs11583680.RData")

#' # Load outcome data
#' ***
BCAC = fread(GWAS_data_BCAC)
BCAC = BCAC[chromosome == 1,]
BCAC = BCAC[position %in% IVData$pos_b37,]

BC_FinnGen_UKB = fread(GWAS_data_finngen)
setnames(BC_FinnGen_UKB,"#CHR","CHR")
BC_FinnGen_UKB = BC_FinnGen_UKB[CHR == 1,]
BC_FinnGen_UKB = BC_FinnGen_UKB[POS %in% IVData$pos_b38,]

#' reformate BCAC
BCAC[, rsID := "rs11583680"]
BCAC[,pos_b38 := BC_FinnGen_UKB$POS]
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

#' ## BC data 
#' 
#' Sample sizes are taken from the [FinnGen + UKB webpage](https://public-metaresults-fg-ukbb.finngen.fi/) (filtering for "Malignant neoplasm of breast (controls excluding all cancers)")
#' 
BC_FinnGen_UKB[,pos_b37 := BCAC$position]
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

#' Save data 
save(IVData, outcomeData,file = "../results/08_SumStats_rs11583680.RData")

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
MR_ratio[p_IV1<0.05]

#' # Save data ####
#' ***
save(MR_ratio, file = "../results/08_MR_ratio_rs11583680.RData")

#' # Forest plot ####
#' ***
#' ## Breast cancer survival 

plotData = copy(MR_ratio)
plotData = plotData[pval<0.05,]
names(plotData)[13:19] = paste(names(plotData)[13:19], "outcome", sep="_")
plotData = plotData[phenotype_outcome == "BCS_meta",]

plotData[,lowerCI95 := beta_IV-1.96*se_IV2]
plotData[,upperCI95 := beta_IV+1.96*se_IV2]

plotData$` ` <- paste(rep(" ", 50), collapse = " ")
plotData$`Estimate \n[95% CI]` <- ifelse(is.na(plotData$se_IV2), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta_IV, plotData$lowerCI95, plotData$upperCI95))

setorder(plotData,beta_IV)
plotData$phenotype
plotData$Exposure = c("GE (brain, cerebellum)",
                      "GE (brain, cerebellar hemisphere)",
                      "GE (pancreas)",
                      "GE (nerve, tibial)",
                      "GE (lung)",
                      "PE (men, statin-treated)",
                      "PE (statin-treated)",
                      "PE (men)",
                      "PE (women, statin-free)",
                      "PE (statin-free)",
                      "PE (women)",
                      "PE (men, statin-free)",
                      "PE (women, statin-treated)")

setnames(plotData,"Exposure","Exposure subgroups")
p2<- forest(plotData[,c(29,27,28)],
            est = plotData$beta_IV,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            xlab = "Causal effect estimate",
            title = "MR of PCSK9 levels on breast cancer survival (rs11583680)")

plot(p2)

filename = paste0("../results/08_ForestPlots_rs11583680_BCS.png")
png(filename = filename,width = 1800, height = 900, res=200)
plot(p2)
dev.off()

#' ## Breast cancer

plotData = copy(MR_ratio)
plotData = plotData[pval<0.05,]
names(plotData)[13:19] = paste(names(plotData)[13:19], "outcome", sep="_")
plotData = plotData[phenotype_outcome == "BCP_meta",]

plotData[,lowerCI95 := beta_IV-1.96*se_IV2]
plotData[,upperCI95 := beta_IV+1.96*se_IV2]

plotData$` ` <- paste(rep(" ", 50), collapse = " ")
plotData$`Estimate \n[95% CI]` <- ifelse(is.na(plotData$se_IV2), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 plotData$beta_IV, plotData$lowerCI95, plotData$upperCI95))

setorder(plotData,-beta_IV)
plotData$phenotype
plotData$Exposure = c("GE (brain, cerebellum)",
                      "GE (brain, cerebellar hemisphere)",
                      "GE (pancreas)",
                      "GE (nerve, tibial)",
                      "GE (lung)",
                      "PE (men, statin-treated)",
                      "PE (statin-treated)",
                      "PE (men)",
                      "PE (women, statin-free)",
                      "PE (statin-free)",
                      "PE (women)",
                      "PE (men, statin-free)",
                      "PE (women, statin-treated)")

setnames(plotData,"Exposure","Exposure subgroups")
p2<- forest(plotData[,c(29,27,28)],
            est = plotData$beta_IV,
            lower = plotData$lowerCI95, 
            upper = plotData$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            xlab = "Causal effect estimate",
            title = "MR of PCSK9 levels on breast cancer (rs11583680)")

plot(p2)

filename = paste0("../results/08_ForestPlots_rs11583680_BC.png")
png(filename = filename,width = 1800, height = 900, res=200)
plot(p2)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

