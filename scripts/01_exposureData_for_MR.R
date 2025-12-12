#' ---
#' title: "Check GTEx v10 data"
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
#' In this script, I extract all necessary SNP data for the exposures: PCSK9 gene expression ([GTEx v10](https://gtexportal.org/home/downloads/adult-gtex/qtl)), PCSK9 protein levels ([Pott et al.](https://zenodo.org/records/10600167)), and LDL-C levels ([GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/))
#' 
#' The GTEx data will be LD-based pruned using [FUMA](https://fuma.ctglab.nl/). For the PCSK9 protein levels, I will use the same SNPs as reported independent in the respective publication. For LDL-C, I will perform position-based pruning (I cannot upload the data to FUMA, as FUMA does not tolerate p-value of 0). However, pruning will be performed after I filter for outcome availability.   
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = T

source("../SourceFile.R")

#' # Load GTEx data ####
#' ***
#' GTEx data is given in GRCh38, but I want GRCh37 (BCAC data is given with chr pos37 only!). I will load the GLGC data, and map the rsIDs. Missing SNPs will be discarded. 
#' 
GTEx = read.csv(GTEx_PCSK9)
setDT(GTEx)
GTEx = GTEx[P.Value<5e-8 | SNP.Id == "rs562556",]

GTEx[, table(Tissue)]
GTEx = GTEx[Tissue != "Testis"]

data0 = fread(GLGC_LDLC_all)
data0 = data0[rsID %in% GTEx$SNP.Id,]
GTEx = GTEx[SNP.Id %in% data0$rsID]
matched = match(GTEx$SNP.Id,data0$rsID)
GTEx[,posb37 := data0[matched,POS_b37]]

GTEx[,tissue2 := gsub(" - ","_",Tissue)]
GTEx[,tissue2 := gsub(" ","_",tissue2)]
GTEx[,tissue2 := gsub("[())]","",tissue2)]
GTEx[, .N, by = tissue2]

myTissues = unique(GTEx$tissue2)
myTissues_GTExNotation = unique(GTEx$Tissue)
myGeneID = unique(GTEx$Gencode.Id)

dumTab1 = foreach(i = 1:length(myTissues))%do%{
  #i=1
  data0 = read_parquet(paste0(GTEx_tissues,myTissues[i],".v10.eQTLs.signif_pairs.parquet"))
  setDT(data0)
  
  # filter data
  data0 = data0[gene_id == myGeneID,]
  data0 = data0[variant_id %in% GTEx[tissue2 == myTissues[i],Variant.Id]]
  
  # get information on sample size
  data0[,maf := af]
  data0[af>0.5,maf := 1-af]
  data0[,sampleSize := round(ma_count/(2*maf),0)]
  
  # get allele coding
  dummy = unlist(strsplit(data0$variant_id,"_"))
  data0[,chr := dummy[seq(1,length(dummy),5)]]
  data0[,chr := gsub("chr","",chr)]
  data0[,chr := as.numeric(chr)]
  data0[,pos_hg38 := dummy[seq(2,length(dummy),5)]]
  data0[,other_allele := dummy[seq(3,length(dummy),5)]]
  data0[,effect_allele := dummy[seq(4,length(dummy),5)]]
  setnames(data0,"af","eaf")
  
  # get relevant columns
  matched = match(data0$variant_id,GTEx$Variant.Id)
  stopifnot(data0$variant_id == GTEx[matched,Variant.Id])
  data0[,rsID := GTEx[matched,SNP.Id]]
  data0[,pos_b37 := GTEx[matched,posb37]]
  data0[,tissue := myTissues_GTExNotation[i]]
  data0[,phenotype := "PCSK9 GE levels"]
  
  # select relevant columns
  names(data0)
  data0 = data0[,c(19,2,15,20,16:18,22,21,4,14,8,9,7)]
  names(data0) = c("rsID","study_ID","chr","pos_b37","pos_b38","OA","EA",
                   "phenotype","setting","EAF","nSamples","beta","se","pval")
  
  # return
  data0[rsID == "rs562556",setting := paste(setting, "ratio")]
  data0
}
GTExData = rbindlist(dumTab1)
GTExData[pval<5e-8,.N,by=setting]
GTExData[,.N,by=setting]
GTExData[,min(pval),by=setting]
GTExData[,max(pval),by=setting]
GTExData[,setting := gsub(" - "," ",setting)]
GTExData[,setting := gsub("[()]","",setting)]
GTExData[,source := "GTExv10"]
GTExData[,MRapproach := "QTL"]
GTExData[grepl("ratio",setting),MRapproach := "rs562556"]
GTExData[,setting := gsub(" ratio","",setting)]

save(GTExData,file = "../temp/GTExV10_potentialIVs.RData")

#' # Load Pott et al. data ####
#' ***
#' Here, I want to use the same instruments as in the publication: rs11591147, rs693668, rs11583680, and rs2495491. In addition, I extract rs562556 for the MR-ratio approach. 
#' 
mySNPs = c("rs11591147", "rs693668", "rs11583680", "rs2495491", "rs562556")
myFiles = c(Pott_PCSK9_females,Pott_PCSK9_males)

dumTab2 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data0 = fread(myFiles[i])
  
  # filter data for relevant SNPs
  data0 = data0[chr == 1,]
  data0[,rsID := gsub(":.*","",markername)]
  data0 = data0[rsID %in% mySNPs,]
  data0[,setting := gsub("PCSK9_","",phenotype)]
  data0[,phenotype := "PCSK9 PE levels"]
  
  # select relevant columns
  names(data0)
  data0 = data0[,c(17,1,2,3,5,4,16,18,6,8,10,11,12)]
  names(data0) = c("rsID","study_ID","chr","pos_b37","OA","EA",
                   "phenotype","setting","EAF","nSamples","beta","se","pval")
  
  # return
  data0[rsID == "rs562556",setting := paste(setting, "ratio")]
  data0
}
pQTLData = rbindlist(dumTab2)
pQTLData[pval<5e-8,.N,by=setting]
pQTLData[,.N,by=setting]
pQTLData[,min(pval),by=setting]
pQTLData[,max(pval),by=setting]

#' Meta-analyse males and females
dumTab2 = foreach(i = 1:length(mySNPs))%do%{
  #i=1
  filt = pQTLData$rsID == mySNPs[i]
  meta = metagen(TE = pQTLData[filt,beta],
                 seTE = pQTLData[filt,se],
                 studlab = pQTLData[filt,setting])
  # summary(meta)
  # forest(meta)
  
  mySamples = pQTLData[filt,nSamples]
  myEAFs = pQTLData[filt,EAF]
  mySamples2 = mySamples/sum(mySamples)
  myEAFs2 = round(sum(myEAFs * mySamples/sum(mySamples)),4)
  
  data = copy(pQTLData)
  data = data[!grepl("females",setting) & rsID == mySNPs[i],]
  data[, beta := meta$TE.fixed]
  data[, se := meta$seTE.fixed]
  data[, pval := meta$pval.fixed]
  data[, nSamples := sum(mySamples)]
  data[, EAF := myEAFs2]
  data[, setting := "all"]
  
  data
}
pQTLData2 = rbindlist(dumTab2)
pQTLData2[rsID=="rs562556",setting := paste(setting,"ratio")]

pQTLData = rbind(pQTLData[grepl("females",setting),],pQTLData2)
setorder(pQTLData,setting,chr,pos_b37)

pQTLData[,source := "Pott et al."]
pQTLData[,MRapproach := "QTL"]
pQTLData[grepl("ratio",setting),MRapproach := "rs562556"]
pQTLData[,setting := gsub(" ratio","",setting)]
pQTLData
save(pQTLData,file = "../temp/Pott_potentialIVs_MR.RData")

#' # Load GLGC data ####
#' ***
myFiles = c(GLGC_LDLC_females,GLGC_LDLC_all)
mySettings = c("females","all")
mySNPs_PCSK9 = fread("../temp/SNPList_PCSK9_Yang2024.txt",header = F)
mySNPs_HMGCR = fread("../temp/SNPList_HMGCR_Yang2024.txt",header = F)

dumTab3 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  data0 = fread(myFiles[i])
  
  # filter data for relevant SNPs
  stopifnot(data0[rsID == "rs562556",pvalue_neg_log10_GC]> -log10(5e-8))
  data0 = data0[pvalue_neg_log10_GC> -log10(5e-8)]
  data0 = data0[POOLED_ALT_AF>0.01,]
  data0 = data0[!is.na(rsID),]
  
  data0[,setting := mySettings[i]]
  data0[,phenotype := "LDL-C levels"]
  
  # select relevant columns
  names(data0)
  data0 = data0[,c(1,1:5,16,15,8,6,9,10,12)]
  names(data0) = c("rsID","study_ID","chr","pos_b37","OA","EA",
                   "phenotype","setting","EAF","nSamples","beta","se","pval")
  save(data0,file = paste0("../temp/GLGC_potentialIVs_unpruned_",mySettings[i],".RData"))
  
  # return
  data0
}
LDLCData = rbindlist(dumTab3)
LDLCData[,.N,by=setting]
LDLCData[grepl("female",setting),source := "Kanoni et al."]
LDLCData[grepl("all",setting),source := "Graham et al."]
LDLCData[,MRapproach := "QTL"]
LDLCData[rsID == "rs562556",MRapproach := "rs562556"]
LDLCData

#' # Merge data ####
#' ***
ExposureData = rbind(pQTLData,GTExData,LDLCData,fill=T,use.names=T)
ExposureData[, absZScore := abs(beta/se)]

ExposureData[,min(absZScore),by=c("phenotype", "setting")]
ExposureData[,max(absZScore),by=c("phenotype", "setting")]

#' I want to harmonize the coding, so that the same SNP for different exposures is coded with the same effect allele throughout. I will use the minor allele for all SNPs. 
#' 
#' Step 1: turn alleles & effects in case effect allele is major allele
ExposureData[,EA_ori := EA]
ExposureData[,OA_ori := OA]
ExposureData[EAF>0.5, beta := beta * (-1)]
ExposureData[EAF>0.5, EA := OA_ori]
ExposureData[EAF>0.5, OA := EA_ori]
ExposureData[EAF>0.5, EAF := 1-EAF]
ExposureData[,EA_ori := NULL]
ExposureData[,OA_ori := NULL]

#' Step 2: get unique SNP list
SNPList = copy(ExposureData)
SNPList = SNPList[,c(1,3:6)]
SNPList = distinct(SNPList)

#' Step 3: get chr_pos_EA_OA and chr_pos_OA_EA and check for overlaps
SNPList[,dumID1 := paste0(chr,"_",pos_b37,"_",EA,"_",OA)]
SNPList[,dumID2 := paste0(chr,"_",pos_b37,"_",OA,"_",EA)]
table(is.element(SNPList$dumID1,SNPList$dumID2))
SNPList[dumID1 %in% dumID2,.N,by=chr]

#' OK, this affects mainly the LDLC SNPs - I have plenty of those, so I will remove them all
problemIDs = SNPList[dumID1 %in% dumID2,rsID]
SNPList = SNPList[!is.element(rsID, problemIDs),]

#' Let's check if there are any duplicated rsIDs left
SNPList[duplicated(rsID),]
problemIDs = SNPList[duplicated(rsID),rsID]
SNPList[rsID %in% problemIDs,]

#' OK, these are triallelic SNPs - just remove them, they always cause trouble
SNPList = SNPList[!is.element(rsID, problemIDs),]

#' Let's check if there are any duplicated positions left
SNPList[,chrPos := paste(chr,pos_b37,sep=":")]
SNPList[duplicated(chrPos),]
problemIDs = SNPList[duplicated(chrPos),chrPos]
SNPList[chrPos %in% problemIDs,]

#' OK, these are triallelic SNPs with different rsIDs - just remove them, they always cause trouble
SNPList = SNPList[!is.element(chrPos, problemIDs),]

#' Now I filter for those SNPs in the exposure data
ExposureData[,.N,by=c("phenotype","setting")]
ExposureData[rsID %in% SNPList$rsID,.N,by=c("phenotype","setting")]
ExposureData = ExposureData[rsID %in% SNPList$rsID,]

#' # Save data ####
#' ***
ExposureData[,exposureID := paste(phenotype,setting,source,MRapproach,sep=" - ")]
save(ExposureData, file = "../results/01_Exposure_for_MR.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

