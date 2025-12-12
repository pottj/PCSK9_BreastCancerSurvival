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
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = T

source("../SourceFile.R")

#' # Load exposure data ####
#' ***
load("../results/01_Exposure_for_MR.RData")
ExposureData[,dumID1 := paste(chr,pos_b37,EA,OA,sep="_")]
ExposureData[,dumID2 := paste(chr,pos_b37,OA,EA,sep="_")]
ExposureData[,dumID3 := paste0(chr,":",pos_b37,"_",EA,"_",OA)]
ExposureData[,dumID4 := paste0(chr,":",pos_b37,"_",OA,"_",EA)]

SNPList = copy(ExposureData)
SNPList = SNPList[,c(1,3:6,19:22)]
SNPList = distinct(SNPList)
SNPList[,chrPos := paste(chr,pos_b37,sep = "_")]

#' # Load breast cancer survival data from BCAC ####
#' ***
#' Load data
BCAC = fread(BCAC_BCS)
names(BCAC)

#' Remove duplicated chr_pos (potential triallelic SNPs or other artifacts)
BCAC[,chrPos := paste(Chromosome,Position,sep="_")]
badSNPs = BCAC[duplicated(chrPos),chrPos]
BCAC = BCAC[!is.element(chrPos,badSNPs),]

#' Filter for overlapping SNPs
BCAC = BCAC[Variant %in% SNPList$dumID1 | Variant %in% SNPList$dumID2,]
BCAC[,EAF := (exp_freq_a1_iCOGS + exp_freq_a1_OncoArray)/2]
BCAC[exp_freq_a1_iCOGS==0,EAF := exp_freq_a1_OncoArray]
BCAC[is.na(exp_freq_a1_iCOGS),EAF := exp_freq_a1_OncoArray]
BCAC[exp_freq_a1_OncoArray==0,EAF := exp_freq_a1_iCOGS]
BCAC[is.na(exp_freq_a1_OncoArray),EAF := exp_freq_a1_iCOGS]

#' Check alleles
SNPList[,flag_BCS_BCAC := F]
SNPList[,flag_BCS_BCAC := F]
SNPList[dumID1 %in% BCAC$Variant | dumID2 %in% BCAC$Variant,flag_BCS_BCAC := T]
SNPList[,table(flag_BCS_BCAC)]

matched = match(BCAC$chrPos,SNPList$chrPos)
table(is.na(matched))
table(SNPList[matched,chrPos] == BCAC$chrPos)
BCAC[,rsID := SNPList[matched,rsID]]
BCAC[,EA := SNPList[matched,EA]]
BCAC[,OA := SNPList[matched,OA]]

BCAC[EA != a1, EAF := 1-EAF]
BCAC[EA != a1,Beta := Beta*(-1)]

#' Get all relevant columns
BCAC[,phenotype := "BCS"]
BCAC[,setting := "females"]
BCAC[,source := "Morra et al."]
BCAC[,nSamples := 91686]
BCAC[,nCases := 7531]
names(BCAC)[c(20,2,3,22,21,23,24,25,19,26,27,14,15,16)]
BCAC = BCAC[,c(20,2,3,22,21,23,24,25,19,26,27,14,15,16)]
head(BCAC)
names(BCAC) = c("rsID","chr","pos_b37","OA","EA",
                "phenotype","setting","source","EAF","nSamples","nCases","beta","se","pval")

#' Sanity check
matched = match(ExposureData$rsID,BCAC$rsID)
table(is.na(matched))
plot(ExposureData$EAF,BCAC[matched,EAF])

#' there are still some SNPs with high difference in allele frequency. I will remove them too
diff = abs(ExposureData$EAF - BCAC[matched,EAF])
hist(diff)
table(diff>0.1)
ExposureData[diff>0.1,table(phenotype,setting)]
badSNPs = ExposureData[diff>0.1,rsID]
SNPList[rsID %in% badSNPs,flag_BCS_BCAC := F]
BCAC = BCAC[!is.element(rsID, badSNPs),]
matched = match(ExposureData$rsID,BCAC$rsID)
table(is.na(matched))
plot(ExposureData$EAF,BCAC[matched,EAF])

save(SNPList,BCAC,file = "../temp/Outcomes_unpruned.RData")

#' # Load breast cancer survival data from Mei et al. ####
#' ***
#' Load data
BCS_Mei = read_xlsx("../temp/Mei_2024_rs561556_BCS.xlsx")
setDT(BCS_Mei)
BCS_Mei

#' Meta- analysis of EUR only
filt = BCS_Mei$comment == "EUR only"
meta = metagen(TE = BCS_Mei[filt,log(HR)], 
               lower = BCS_Mei[filt,log(CI_low)], upper = BCS_Mei[filt,log(CI_up)],
               sm = "HR", studlab = BCS_Mei[filt,Study])
summary(meta)

forest(meta)

filename = paste0("../results/02_ForestPlot_Mei_data.png")
png(filename = filename,width = 1900, height = 600, res=200)
forest(meta)
dev.off()

#' Summarize in the same way as BCAC. Note: Mei et al. used A as effect allele, so I have to switch the direction of the beta estimate
BCS = copy(BCAC)
BCS = BCS[rsID == "rs562556",]
BCS[,source := "Mei et al."]
BCS[,nSamples := NA]
BCS[,nCases := NA]
BCS[,beta := -meta$TE.fixed]
BCS[,se := meta$seTE.fixed]
BCS[,pval := meta$pval.fixed]

BCS = rbind(BCS,BCS,BCS,BCS,BCS)
BCS[2:5, source := meta$studlab]
BCS[2:5, beta := -meta$TE]
BCS[2:5, se := meta$seTE]
BCS[2:5, pval := meta$pval]
BCS[,source := gsub(", .*","",source)]
BCS

save(SNPList,BCAC,BCS,file = "../temp/Outcomes_unpruned.RData")

#' # Load breast cancer data ####
#' ***
#' Load data
BC = fread(FinnGenUKB_BC)
setnames(BC,"#CHR","CHR")

#' Remove duplicated chr_pos (potential triallelic SNPs or other artifacts)
BC[,chrPos := paste(CHR,POS,sep="_")]
badSNPs = BC[duplicated(chrPos),chrPos]
BC = BC[!is.element(chrPos,badSNPs),]

#' Filter for overlapping SNPs
BC = BC[rsid %in% SNPList$rsID ,]
BC[,table(duplicated(rsid))]
BC[,EAF := (FINNGEN_af_alt + UKBB_af_alt)/2]
BC[is.na(UKBB_af_alt),EAF := FINNGEN_af_alt]
BC[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]

#' Check alleles
SNPList[,flag_BC := F]
SNPList[rsID %in% BC$rsid,flag_BC := T]
SNPList[,table(flag_BC)]

matched = match(BC$rsid,SNPList$rsID)
table(is.na(matched))
table(SNPList[matched,rsID] == BC$rsid)
BC[,pos_b37 := SNPList[matched,pos_b37]]
BC[,EA := SNPList[matched,EA]]
BC[,OA := SNPList[matched,OA]]

BC[EA != ALT, EAF := 1-EAF]
BC[EA != ALT, all_inv_var_meta_beta := all_inv_var_meta_beta*(-1)]

#' Get all relevant columns
BC[,phenotype := "BC"]
BC[,setting := "females"]
BC[,source := "FinnGen + UKB"]
BC[,nSamples := 18786 + 11807 + 182927 + 205913]
BC[,nCases := 18786 + 11807]
names(BC)[c(22,1,25,2,27,26,28,29,30,24,31,32,17,18,19)]
BC = BC[,c(22,1,25,2,27,26,28,29,30,24,31,32,17,18,19)]
head(BC)
names(BC) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
              "phenotype","setting","source","EAF","nSamples","nCases","beta","se","pval")

#' Sanity check
matched = match(ExposureData$rsID,BC$rsID)
table(is.na(matched))
plot(ExposureData$EAF,BC[matched,EAF])

#' there are still some SNPs with high difference in allele frequency. I will remove them too
diff = abs(ExposureData$EAF - BC[matched,EAF])
hist(diff)
table(diff>0.1)
ExposureData[diff>0.1,table(phenotype,setting)]
badSNPs = ExposureData[diff>0.1,rsID]
SNPList[rsID %in% badSNPs,flag_BC := F]
BC = BC[!is.element(rsID, badSNPs),]
matched = match(ExposureData$rsID,BC$rsID)
table(is.na(matched))
plot(ExposureData$EAF,BC[matched,EAF])

save(SNPList,BCAC,BCS, BC, file = "../temp/Outcomes_unpruned.RData")

#' # Load parental longevity data ####
#' ***
#' Load data
PLD = fread(Pilling_ParentalLongevity)

#' Remove duplicated chr_pos (potential triallelic SNPs or other artefacts)
PLD[,chrPos := paste(hm_chrom,hm_pos,sep="_")]
badSNPs = PLD[duplicated(chrPos),chrPos]
PLD = PLD[!is.element(chrPos,badSNPs),]

#' Filter for overlapping SNPs
PLD = PLD[hm_rsid %in% SNPList$rsID ,]
PLD[,table(duplicated(hm_rsid))]

#' Check alleles
SNPList[,flag_PLD := F]
SNPList[rsID %in% PLD$hm_rsid,flag_PLD := T]
SNPList[,table(flag_PLD)]

matched = match(PLD$hm_rsid,SNPList$rsID)
table(is.na(matched))
table(SNPList[matched,rsID] == PLD$hm_rsid)
PLD[,pos_b37 := SNPList[matched,pos_b37]]
PLD[,EA := SNPList[matched,EA]]
PLD[,OA := SNPList[matched,OA]]

PLD[,EAF := hm_effect_allele_frequency]
PLD[,beta2 := hm_beta]
PLD[EA != hm_effect_allele, EAF := 1-EAF]
PLD[EA != hm_effect_allele, beta2 := beta2*(-1)]

#' Get all relevant columns
PLD[,phenotype := "PLD"]
PLD[,setting := "all"]
PLD[,source := "Pilling et al."]
PLD[,nSamples := 208118]
names(PLD)[c(2,3,26,4,28,27,31,32,33,29,34,30,20,21)]
PLD = PLD[,c(2,3,26,4,28,27,31,32,33,29,34,30,20,21)]
head(PLD)
names(PLD) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
              "phenotype","setting","source","EAF","nSamples","beta","se","pval")

#' Sanity check
matched = match(ExposureData$rsID,PLD$rsID)
table(is.na(matched))
plot(ExposureData$EAF,PLD[matched,EAF])

#' there are still some SNPs with high difference in allele frequency. I will remove them too
diff = abs(ExposureData$EAF - PLD[matched,EAF])
hist(diff)
table(diff>0.1)
ExposureData[diff>0.1,table(phenotype,setting)]
badSNPs = ExposureData[diff>0.1,rsID]
SNPList[rsID %in% badSNPs,flag_PLD := F]
PLD = PLD[!is.element(rsID, badSNPs),]
matched = match(ExposureData$rsID,PLD$rsID)
table(is.na(matched))
plot(ExposureData$EAF,PLD[matched,EAF])

save(SNPList,BCAC,BCS, BC, PLD, file = "../temp/Outcomes_unpruned.RData")

#' # Load coronary artery disease data (females) ####
#' ***
#' Load data
CAD_f = fread(Aragam_CAD_females)

#' Remove duplicated chr_pos (potential triallelic SNPs or other artefacts)
CAD_f[,chrPos := paste(CHR,BP,sep="_")]
badSNPs = CAD_f[duplicated(chrPos),chrPos]
CAD_f = CAD_f[!is.element(chrPos,badSNPs),]

#' Filter for overlapping SNPs
CAD_f = CAD_f[rsid_ukb %in% SNPList$rsID ,]
CAD_f[,table(duplicated(rsid_ukb))]

#' Check alleles
SNPList[,flag_CAD_fem := F]
SNPList[rsID %in% CAD_f$rsid_ukb,flag_CAD_fem := T]
SNPList[,table(flag_CAD_fem)]

matched = match(CAD_f$rsid_ukb,SNPList$rsID)
table(is.na(matched))
table(SNPList[matched,rsID] == CAD_f$rsid_ukb)
table(SNPList[matched,pos_b37] == CAD_f$BP)
CAD_f[,pos_b37 := SNPList[matched,pos_b37]]
CAD_f[,EA := SNPList[matched,EA]]
CAD_f[,OA := SNPList[matched,OA]]

CAD_f[,EAF := female_eaf]
CAD_f[,beta2 := female_beta]
CAD_f[EA != reference_allele, EAF := 1-EAF]
CAD_f[EA != reference_allele, beta2 := beta2*(-1)]

#' Get all relevant columns
CAD_f[,phenotype := "CAD"]
CAD_f[,setting := "females"]
CAD_f[,source := "Aragam et al."]
names(CAD_f)[c(38,41,44,46,45,49,50,51,47,35,48,29,33)]
CAD_f = CAD_f[,c(38,41,44,46,45,49,50,51,47,35,48,29,33)]
head(CAD_f)
names(CAD_f) = c("rsID","chr","pos_b37","OA","EA",
               "phenotype","setting","source","EAF","nSamples","beta","se","pval")

#' Sanity check
matched = match(ExposureData$rsID,CAD_f$rsID)
table(is.na(matched))
plot(ExposureData$EAF,CAD_f[matched,EAF])

#' there are still some SNPs with high difference in allele frequency. I will remove them too
diff = abs(ExposureData$EAF - CAD_f[matched,EAF])
hist(diff)
table(diff>0.1)
ExposureData[diff>0.1,table(phenotype,setting)]
badSNPs = ExposureData[diff>0.1,rsID]
SNPList[rsID %in% badSNPs,flag_CAD_fem := F]
CAD_f = CAD_f[!is.element(rsID, badSNPs),]
matched = match(ExposureData$rsID,CAD_f$rsID)
table(is.na(matched))
plot(ExposureData$EAF,CAD_f[matched,EAF])

save(SNPList,BCAC,BCS, BC, PLD, CAD_f, file = "../temp/Outcomes_unpruned.RData")

#' # Load coronary artery disease data (sex-combined) ####
#' ***
#' Load data
CAD_a = fread(Aragam_CAD_all)

#' Remove duplicated chr_pos (potential triallelic SNPs or other artefacts)
CAD_a[,chrPos := paste(CHR,BP,sep="_")]
badSNPs = CAD_a[duplicated(chrPos),chrPos]
CAD_a = CAD_a[!is.element(chrPos,badSNPs),]

#' Filter for overlapping SNPs
CAD_a = CAD_a[MarkerName %in% SNPList$dumID3 | MarkerName %in% SNPList$dumID4 ,]
CAD_a[,table(duplicated(MarkerName))]

#' Check alleles
SNPList[,flag_CAD_all := F]
SNPList[dumID3 %in% CAD_a$MarkerName | dumID4 %in% CAD_a$MarkerName,flag_CAD_all := T]
SNPList[,table(flag_CAD_all)]

matched = match(CAD_a$chrPos,SNPList$chrPos)
table(is.na(matched))
table(SNPList[matched,chrPos] == CAD_a$chrPos)
table(SNPList[matched,dumID3] == CAD_a$MarkerName)
table(SNPList[matched,dumID4] == CAD_a$MarkerName)
table(SNPList[matched,pos_b37] == CAD_a$BP)
CAD_a[,rsID := SNPList[matched,rsID]]
CAD_a[,EA := SNPList[matched,EA]]
CAD_a[,OA := SNPList[matched,OA]]

CAD_a[,EAF := Freq1]
CAD_a[,beta2 := Effect]
CAD_a[EA != toupper(Allele1), EAF := 1-EAF]
CAD_a[EA != toupper(Allele1), beta2 := beta2*(-1)]

#' Get all relevant columns
CAD_a[,phenotype := "CAD"]
CAD_a[,setting := "all"]
CAD_a[,source := "Aragam et al."]
names(CAD_a)[c(23,2,3,25,24,28,29,30,26,20,18,27,11,12)]
CAD_a = CAD_a[,c(23,2,3,25,24,28,29,30,26,20,18,27,11,12)]
head(CAD_a)
names(CAD_a) = c("rsID","chr","pos_b37","OA","EA",
                 "phenotype","setting","source","EAF","nSamples","nCases","beta","se","pval")

#' Sanity check
matched = match(ExposureData$rsID,CAD_a$rsID)
table(is.na(matched))
plot(ExposureData$EAF,CAD_a[matched,EAF])

#' there are still some SNPs with high difference in allele frequency. I will remove them too
diff = abs(ExposureData$EAF - CAD_a[matched,EAF])
hist(diff)
table(diff>0.1)
ExposureData[diff>0.1,table(phenotype,setting)]
badSNPs = ExposureData[diff>0.1,rsID]
SNPList[rsID %in% badSNPs,flag_CAD_all := F]
CAD_a = CAD_a[!is.element(rsID, badSNPs),]
matched = match(ExposureData$rsID,CAD_a$rsID)
table(is.na(matched))
plot(ExposureData$EAF,CAD_a[matched,EAF])

save(SNPList,BCAC,BCS, BC, PLD, CAD_f, CAD_a, 
     file = "../temp/Outcomes_unpruned.RData")

#' # Final checks ####
#' ***
#' Ideally, I would use the same SNP set per exposure on all outcomes. Then I could at least make a fair comparison between the outcomes. 
#' 
#' I will check if I will lose an outcome because of this. 
#' 
SNPList[,flag := F]
SNPList[flag_BCS_BCAC==T & flag_BC==T & flag_PLD==T & flag_CAD_fem==T & flag_CAD_all==T,flag := T]
SNPList[,table(flag)]

ExposureData[,.N,by=c("phenotype","setting")]

ExposureData[rsID %in% SNPList[flag==T,rsID],.N,by=c("phenotype","setting")]
ExposureData[!is.element(rsID,SNPList[flag==T,rsID]),.N,by=c("phenotype","setting")]

ExposureData[rsID %in% SNPList[flag==T,rsID],max(absZScore),by=c("phenotype","setting")]
ExposureData[!is.element(rsID,SNPList[flag==T,rsID]),max(absZScore),by=c("phenotype","setting")]

#' Summary: 
#' 
#' - all SNPs for the pQTLs still available
#' - still a lot SNPs for the LDL-C analyses available
#' - some less eQTLs, but I only lose 1 tissue (Adipose - Subcutaneous), and only in Nerve - Tibial I lose the best eQTL. 
#' 
#' Decision: I will restrict to those SNPs available in all outcomes!
#' 
ExposureData = ExposureData[rsID %in% SNPList[flag==T,rsID],]

BCAC[,geneticModel := "additive"]
BC[,geneticModel := "additive"]
PLD[,geneticModel := "additive"]
CAD_a[,geneticModel := "additive"]
CAD_f[,geneticModel := "additive"]
BCS[,geneticModel := "recessive"]

OutcomeData = rbind(BCAC,BCS, BC, PLD, CAD_f, CAD_a,fill=T)
OutcomeData = OutcomeData[rsID %in% SNPList[flag==T,rsID],]

ExposureData[rsID=="rs11591147",]
OutcomeData[rsID=="rs11591147",]
ExposureData[rsID=="rs562556",]
OutcomeData[rsID=="rs562556",]

ExposureData[,table(rsID %in% OutcomeData$rsID)]
matched = match(ExposureData$rsID,OutcomeData[phenotype=="BC",rsID])
table(is.na(matched))
table(ExposureData$rsID == OutcomeData[phenotype=="BC",rsID][matched])
ExposureData[,pos_b38 := OutcomeData[phenotype=="BC",pos_b38][matched]]
ExposureData = ExposureData[,c(1,3,4,16,5:8,14,15,9:13,17,18)]

OutcomeData[,table(rsID %in% ExposureData$rsID)]
matched = match(OutcomeData$rsID,ExposureData[,rsID])
table(is.na(matched))
table(OutcomeData$rsID == ExposureData[,rsID][matched])
OutcomeData[,pos_b38 := ExposureData[,pos_b38][matched]]
OutcomeData = OutcomeData[,c(1:3,16,4:7,14,15,8:13)]
OutcomeData[,outcomeID := paste(phenotype,setting,source,geneticModel,sep=" - ")]

#' # Save data ####
#' ***
save(ExposureData,file = "../results/02_Exposure_for_MR.RData")
save(OutcomeData, file = "../results/02_Outcome_for_MR.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

