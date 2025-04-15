#' ---
#' title: "MR of LDL-C on outcomes"
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
#' Used data sources: 
#' 
#' - exposure: LDL-C in women (GLGC data, Kanoni et al. 2022)
#' - outcome 1: breast cancer survival (BCAC data, Morra et al. 2021)
#' - outcome 2: breast cancer (FinnGen & UKB)
#' - outcome 3: longevity, parents' age at death (UKB, Pilling et al. 2017)
#'
#' Instrument selection: 
#' 
#' - genome-wide SNPs, priority pruning by position (exclude SNPs within +- 1Mb of variant with lowest p-value)
#' - same SNPs as in PCSK9 pQTL MR approach (3 SNPs, as one variant is missing in LDLC data set)
#' - HMGCR SNPs
#' - HMGCR and PCSK9 SNPs
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
LDLC = fread(GWAS_data_LDLC)
LDLC

BCAC = fread(GWAS_data_BCAC)
BC_FinnGen_UKB = fread(GWAS_data_finngen_BC)
CAD = fread(GWAS_data_finngen_CAD)
ParentalLongevity_Death = fread(GWAS_data_GWASCatalog)

#' # Filter data ####
#' ***
#' ## Step 1: just genome-wide significant SNPs
LDLC = LDLC[pvalue_neg_log10_GC> -log10(5e-8)]
LDLC[,dumID1 := paste(CHROM, POS_b37,REF,ALT,sep="_")]
LDLC[,dumID2 := paste(CHROM, POS_b37,ALT,REF,sep="_")]
setorder(LDLC,CHROM,POS_b37)

table(is.element(LDLC$dumID1,BCAC$snp))
table(is.element(LDLC$dumID2,BCAC$snp))
BCAC = BCAC[snp %in% LDLC$dumID1 | snp %in% LDLC$dumID2,]
setorder(BCAC,chromosome,position)

BC_FinnGen_UKB = BC_FinnGen_UKB[rsid %in% LDLC$rsID]
BC_FinnGen_UKB = BC_FinnGen_UKB[!is.na(rsid)]
setnames(BC_FinnGen_UKB,"#CHR","chr")
setorder(BC_FinnGen_UKB,chr,POS)

CAD = CAD[rsid %in% LDLC$rsID]
CAD = CAD[!is.na(rsid)]
setnames(CAD,"#CHR","chr")
setorder(CAD,chr,POS)

ParentalLongevity_Death = ParentalLongevity_Death[hm_rsid %in% LDLC$rsID,]
ParentalLongevity_Death = ParentalLongevity_Death[!is.na(hm_rsid)]
ParentalLongevity_Death[,hm_chrom := as.numeric(hm_chrom)]
setorder(ParentalLongevity_Death,hm_chrom,hm_pos)

#' ## Step 2: reduce to SNPs in overlap
LDLC = LDLC[dumID1 %in% BCAC$snp | dumID2 %in% BCAC$snp,]
LDLC = LDLC[rsID %in% BC_FinnGen_UKB$rsid,]
LDLC = LDLC[rsID %in% CAD$rsid,]
LDLC = LDLC[rsID %in% ParentalLongevity_Death$hm_rsid,]
duplicatedSNPs = LDLC[duplicated(rsID),rsID]
LDLC[rsID %in% duplicatedSNPs]

#' Okay, those are tri-allelic SNPs, and I do not want them - remove all of them!
LDLC = LDLC[!is.element(rsID,duplicatedSNPs),]

BCAC = BCAC[snp %in% LDLC$dumID1 | snp %in% LDLC$dumID2,]
BC_FinnGen_UKB = BC_FinnGen_UKB[rsid %in% LDLC$rsID]
stopifnot(LDLC$rsID == BC_FinnGen_UKB$rsid)
CAD = CAD[rsid %in% LDLC$rsID]
stopifnot(LDLC$rsID == CAD$rsid)

LDLC[,POS_b38 := BC_FinnGen_UKB[,POS]]
LDLC[,dumID3 := paste(CHROM, POS_b38,REF,ALT,sep=":")]
LDLC[,dumID4 := paste(CHROM, POS_b38,ALT,REF,sep=":")]
LDLC[,dumID5 := paste(CHROM, POS_b38,REF,ALT,sep="_")]
LDLC[,dumID6 := paste(CHROM, POS_b38,ALT,REF,sep="_")]

BC_FinnGen_UKB = BC_FinnGen_UKB[SNP %in% LDLC$dumID3 | SNP %in% LDLC$dumID4,]
CAD = CAD[SNP %in% LDLC$dumID3 | SNP %in% LDLC$dumID4,]
LDLC = LDLC[rsID %in% BC_FinnGen_UKB$rsid,]
BCAC = BCAC[snp %in% LDLC$dumID1 | snp %in% LDLC$dumID2,]
ParentalLongevity_Death = ParentalLongevity_Death[hm_variant_id %in% LDLC$dumID5 | hm_variant_id %in% LDLC$dumID6,]

stopifnot(LDLC$POS_b37 == BCAC$position)
stopifnot(LDLC$CHROM == BCAC$chromosome)
BCAC[,rsID := LDLC$rsID]

#' ## Step 3: harmonize format
#' I want rsID, chr, pos b37, pos b38, OA, EA, phenotype, EAF, nSamples, beta, se, and pval
names(LDLC)
LDLC[,phenotype := "LDLC"]
LDLC[,pval2 := 10^-pvalue_neg_log10]
LDLC = LDLC[,c(1,2,3,17,4,5,22,8,6,9,10,23)]
names(LDLC) = c("rsID","chr","pos_b37","pos_b38","OA","EA","phenotype","EAF","nSamples","beta","se","pval")
LDLC 

names(BCAC)
BCAC[,phenotype := "BCS"]
BCAC[,pos_b38 := LDLC$pos_b38]
BCAC[,nSamples := 91686 + 7531]
BCAC[,nCases := 7531]
BCAC = BCAC[,c(18,2,3,20,4,5,19,7,21,22,14,15,16)]
names(BCAC) = c("rsID","chr","pos_b37","pos_b38","OA","EA","phenotype","EAF","nSamples","nCases","beta","se","pval")
BCAC 

stopifnot(BC_FinnGen_UKB$rsid == LDLC$rsID)
names(BC_FinnGen_UKB) 
BC_FinnGen_UKB[,phenotype := "BC"]
BC_FinnGen_UKB[,nSamples := 11807 + 205913 + 18786 + 182927]
BC_FinnGen_UKB[,nCases := 11807 + 18786]
BC_FinnGen_UKB[,EAF := (FINNGEN_af_alt+UKBB_af_alt)/2]
BC_FinnGen_UKB[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]
BC_FinnGen_UKB[,pos_b37 := LDLC$pos_b37]
BC_FinnGen_UKB = BC_FinnGen_UKB[,c(22,1,26,2,3,4,23,26,24,25,17,18,19)]
names(BC_FinnGen_UKB) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                          "phenotype","EAF","nSamples","nCases","beta","se","pval")

stopifnot(CAD$rsid == LDLC$rsID)
names(CAD) 
CAD[,phenotype := "CAD"]
CAD[,nSamples := 51589 + 31198 + 343079 + 382052]
CAD[,nCases := 51589 + 31198]
CAD[,EAF := (FINNGEN_af_alt+UKBB_af_alt)/2]
CAD[is.na(FINNGEN_af_alt),EAF := UKBB_af_alt]
CAD[,pos_b37 := LDLC$pos_b37]
CAD = CAD[,c(22,1,26,2,3,4,23,26,24,25,17,18,19)]
names(CAD) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                          "phenotype","EAF","nSamples","nCases","beta","se","pval")

stopifnot(ParentalLongevity_Death$hm_rsid == LDLC$rsID)
names(ParentalLongevity_Death)
ParentalLongevity_Death[,phenotype := "PAAD"]
ParentalLongevity_Death[,nSamples := 208118]
ParentalLongevity_Death[,pos_b37 := LDLC$pos_b37]
ParentalLongevity_Death = ParentalLongevity_Death[,c(2,3,27,4,5,6,25,11,26,7,20,21)]
names(ParentalLongevity_Death) = c("rsID","chr","pos_b37","pos_b38","OA","EA",
                                   "phenotype","EAF","nSamples","beta","se","pval")

#' ## Step 4: harmonize alleles
#' 
#' LDLC data stays at it is, outcome data will be changed
#' 
#' SNPs will be excluded if their allele frequency differ more than 0.1
#' 
filt = LDLC$OA == BCAC$OA
table(filt)
BCAC[!filt, OA := LDLC[!filt,OA]]
BCAC[!filt, EA := LDLC[!filt,EA]]
BCAC[!filt, EAF := 1-EAF]
BCAC[!filt, beta := beta * (-1)]
plot(LDLC$EAF,BCAC$EAF)
abline(0,1)
diff = abs(LDLC$EAF - BCAC$EAF)
hist(diff)
filt1 = diff>0.1
table(filt1)

filt = LDLC$OA == BC_FinnGen_UKB$OA
table(filt)
BC_FinnGen_UKB[!filt, OA := LDLC[!filt,OA]]
BC_FinnGen_UKB[!filt, EA := LDLC[!filt,EA]]
BC_FinnGen_UKB[!filt, EAF := 1-EAF]
BC_FinnGen_UKB[!filt, beta := beta * (-1)]
plot(LDLC$EAF,BC_FinnGen_UKB$EAF)
abline(0,1)
diff = abs(LDLC$EAF - BC_FinnGen_UKB$EAF)
hist(diff)
filt2 = diff>0.1
table(filt2)

filt = LDLC$OA == CAD$OA
table(filt)
CAD[!filt, OA := LDLC[!filt,OA]]
CAD[!filt, EA := LDLC[!filt,EA]]
CAD[!filt, EAF := 1-EAF]
CAD[!filt, beta := beta * (-1)]
plot(LDLC$EAF,CAD$EAF)
abline(0,1)
diff = abs(LDLC$EAF - CAD$EAF)
hist(diff)
filt3 = diff>0.1
table(filt3)

filt = LDLC$OA == ParentalLongevity_Death$OA
table(filt)
ParentalLongevity_Death[!filt, OA := LDLC[!filt,OA]]
ParentalLongevity_Death[!filt, EA := LDLC[!filt,EA]]
ParentalLongevity_Death[!filt, EAF := 1-EAF]
ParentalLongevity_Death[!filt, beta := beta * (-1)]
plot(LDLC$EAF,ParentalLongevity_Death$EAF)
abline(0,1)
diff = abs(LDLC$EAF - ParentalLongevity_Death$EAF)
hist(diff)
filt4 = diff>0.1
table(filt4)

filt = filt1 | filt2 | filt3 | filt4
table(filt)

LDLC = LDLC[!filt,]
BCAC = BCAC[!filt,]
BC_FinnGen_UKB = BC_FinnGen_UKB[!filt,]
CAD = CAD[!filt,]
ParentalLongevity_Death = ParentalLongevity_Death[!filt,]

save(LDLC,BCAC,BC_FinnGen_UKB,CAD,ParentalLongevity_Death,
     file = "../temp/08_MR_LDLC_input_beforePruning.RData")

#' # Pruning ####
#' ***
myChrs = unique(as.numeric(LDLC$chr))

result.22 = foreach(s2 = myChrs) %do% {
  # s2 = myChrs[2]
  subdata2 = copy(LDLC)
  subdata2 = subdata2[chr == s2, ]
  
  setkey(subdata2, pos_b38)
  
  if(dim(subdata2)[1]<=1){
    subdata2[, keep := T]
    subdata2[, NR_SNPs := 0]
  }else{
    subdata2[, keep := NA]
    subdata2[, NR_SNPs := as.numeric(NA)]
    
    smallestDist = getSmallestDist(subdata2[, pos_b38])
    while(smallestDist < 500000) {
      minP = min(subdata2[is.na(keep), pval])
      mybp_hg38 = subdata2[minP == pval & is.na(keep), pos_b38]
      if(length(mybp_hg38)>1){
        mybp_hg38 = mybp_hg38[1]
      }
      subdata2[pos_b38 == mybp_hg38, keep := T]
      
      #filter for SNPs that can stay within the set (outside the +- 500 kb range or keep==T)
      myFilt = (subdata2[, pos_b38] < (mybp_hg38 - 500000)) | 
        (subdata2[, pos_b38] > (mybp_hg38 + 500000)) | 
        subdata2[, keep] 
      myFilt[is.na(myFilt)] = FALSE
      subdata2 = subdata2[myFilt == TRUE, ]
      
      subdata2[pos_b38 == mybp_hg38, NR_SNPs := sum(myFilt==F)]
      smallestDist = getSmallestDist(subdata2[, pos_b38])
    }
    
    #stopifnot(sum(is.na(subdata2[,keep])) <= 1)
    subdata2[is.na(keep), NR_SNPs := 0]
    subdata2[is.na(keep), keep := TRUE]
  }
  
  subdata2
}
result.22 = rbindlist(result.22)
table(result.22$chr)
result.22
result.22 = result.22[NR_SNPs>9,]

LDLC = LDLC[rsID %in% result.22$rsID,]
BCAC = BCAC[rsID %in% result.22$rsID,]
BC_FinnGen_UKB = BC_FinnGen_UKB[rsID %in% result.22$rsID,]
CAD = CAD[rsID %in% result.22$rsID,]
ParentalLongevity_Death = ParentalLongevity_Death[rsID %in% result.22$rsID,]

save(LDLC,BCAC,BC_FinnGen_UKB,CAD,ParentalLongevity_Death,
     file = "../temp/08_MR_LDLC_input_Pruned.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
