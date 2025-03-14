#' ---
#' title: "Sup Fig: Forest plot per outcome"
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
load("../results/SupTables.RData")

MR_IVW = copy(stab2)
MR_ratio = copy(stab3)

#' # Filter ####
#' ***
MR_IVW[,unique(outcome)]
MR_IVW[outcome=="BreastCancer [FinnGen+UKB]", outcome := "Breast Cancer [FinnGen+UKB]"]
MR_IVW[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]
MR_IVW[outcome=="AgeAtDeath [UKB]", outcome := "Parents' age at death [UKB]"]

MR_ratio[,unique(outcome)]
MR_ratio[outcome=="BreastCancer [FinnGen+UKB]", outcome := "Breast Cancer [FinnGen+UKB]"]
MR_ratio[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]
MR_ratio[outcome=="AgeAtDeath [UKB]", outcome := "Parents' age at death [UKB]"]
MR_ratio[outcome=="BCSurvival [Mei et al.]", outcome := "BC Survival [Mei et al.]"]

#' Merge data sets - what columns do I really need?
MR_ratio[,nSNPs := 1]
setnames(MR_ratio,"beta_ratio","beta_IVW")
setnames(MR_ratio,"se_ratio","se_IVW")
setnames(MR_ratio,"pval_ratio","pval_IVW")
setnames(MR_ratio,"Fstat","FStat")
MR_merged = rbind(MR_IVW[,1:7],MR_ratio[,c(7,13,23,19,20,21,22)])
MR_merged[,method := c(rep("MR-IVW",24),rep("MR-ratio",20))]

#' add nice name for phenotype
unique(MR_merged$exposure)
MR_merged[exposure=="Brain_Cerebellum", exposure:="Brain Cb"]
MR_merged[exposure=="Brain_Cerebellar_Hemisphere", exposure:="Brain CH"]
MR_merged[exposure=="Nerve_Tibial", exposure:="Nerve tib"]
MR_merged[exposure=="Esophagus_Muscularis", exposure:="Esophagus musc"]

MR_merged[!grepl("PCSK9",exposure), exposure := paste0("GE - ",exposure)]
MR_merged[grepl("PCSK9",exposure), exposure := gsub("PCSK9_","PE - ",exposure)]
MR_merged[exposure=="PE - free", exposure := "PE - statin-free individuals"]

myOutcomes = unique(MR_merged$outcome)
myOutcomes
myOutcomes2 = c("BC_FinnGenUKB","BCS_BCAC","Parents_AgeAtDeath_UKB","BCS_Mei")
myEstimate = c("OR","HR","beta","HR")
myColor_GE = c("#FBE3D6","#C2F1C8","#F2CFEE","#D9F2D0")
myColor_PE = c("#F2AA84","#47D45A","#D86ECC","#8ED973")

MR_merged[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_merged[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_merged$` ` <- paste(rep(" ", 50), collapse = " ")
MR_merged$`Estimate [95% CI]` <- ifelse(is.na(MR_merged$se_IVW), "",
                                          sprintf("%.2f [%.2f, %.2f]",
                                                  MR_merged$beta_IVW, MR_merged$lowerCI95, MR_merged$upperCI95))

#' # Plots per outcome
#' ***
for(i in 1:length(myOutcomes)){
  #i=1
  plotData1 = copy(MR_merged)
  plotData1 = plotData1[outcome == myOutcomes[i]]
  setorder(plotData1,method,exposure)
  plotData1[,Fstat2 := round(FStat,1)]
  plotData1[,Fstat2 := as.character(Fstat2)]
  
  dummy = data.table(exposure = c("Gene expression - MR-IVW",
                                  "Protein expression - MR-IVW",
                                  "Gene expression - MR-ratio",
                                  "Protein expression - MR-ratio"))
  if(i==4)  dummy = data.table(exposure = c("Gene expression - MR-ratio",
                                            "Protein expression - MR-ratio"))
  
  plotData1 = rbind(plotData1,dummy, fill=T)
  if(i==4){
    plotData1 = plotData1[c(6,1:3,7,4,5)]
  }else{
    plotData1 = plotData1[c(14,1:6,15,7,8,16,9:11,17,12,13)]
  }
  plotData1[is.na(outcome),` ` := ""]
  plotData1[is.na(outcome),`Estimate [95% CI]`:= ""]
  plotData1[is.na(outcome),Fstat2 := ""]
  plotData1[!is.na(outcome),exposure := paste0("   ",exposure)]
  
  dummy2 = plotData1$exposure
  if(i==4){
    dummy2[c(1,5)] = "white"
  }else{
    dummy2[c(1,8,11,15)] = "white"
  }
  dummy2[grepl("GE",dummy2)] = myColor_GE[i]
  dummy2[grepl("PE",dummy2)] = myColor_PE[i]
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))
  
  plotData1[,exposure := gsub("GE - ","",exposure)]
  plotData1[,exposure := gsub("PE - ","",exposure)]
  setnames(plotData1,"exposure","Tissue / subgroup")
  setnames(plotData1,"Fstat2","F-stat")
  
  myTitle = paste0("PCSK9 levels on ",myOutcomes[i])
  if(myEstimate[i]=="OR") myXlab = "log(OR) per PCSK9 unit increase" 
  if(myEstimate[i]=="HR") myXlab = "log(HR) per PCSK9 unit increase"
  if(myEstimate[i]=="beta") myXlab = "beta per PCSK9 unit increase" 
  
  p1 <- forest(plotData1[,c(1,11,12,13)],
               est = plotData1$beta_IVW,
               lower = plotData1$lowerCI95, 
               upper = plotData1$upperCI95,
               sizes = 0.5,
               ci_column = 2,
               ref_line = 0,
               theme = tm1,
               xlab = myXlab,
               #xlim = c(-0.7,0.3),
               title = myTitle)
  
  plot(p1)
  
  filename = paste0("../results/ForestPlots/SupFig1_",myOutcomes2[i],".png")
  myHeight = 1200
  if(i==4) myHeight = 700
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p1)
  dev.off()
  
}

#' # Plots per outcome - MR-ratio
#' ***
for(i in 1:length(myOutcomes)){
  #i=1
  plotData1 = copy(MR_merged)
  plotData1 = plotData1[outcome == myOutcomes[i]]
  plotData1 = plotData1[method == "MR-ratio",]
  setorder(plotData1,method,exposure)
  plotData1[,Fstat2 := round(FStat,1)]
  plotData1[,Fstat2 := as.character(Fstat2)]
  
  dummy = data.table(exposure = c("Gene expression",
                                  "Protein expression"))

  plotData1 = rbind(plotData1,dummy, fill=T)
  plotData1 = plotData1[c(6,1:3,7,4,5)]
  plotData1[is.na(outcome),` ` := ""]
  plotData1[is.na(outcome),`Estimate [95% CI]`:= ""]
  plotData1[is.na(outcome),Fstat2 := ""]
  plotData1[!is.na(outcome),exposure := paste0("   ",exposure)]
  
  dummy2 = plotData1$exposure
  dummy2[c(1,5)] = "white"
  dummy2[grepl("GE",dummy2)] = myColor_GE[i]
  dummy2[grepl("PE",dummy2)] = myColor_PE[i]
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))
  
  plotData1[,exposure := gsub("GE - ","",exposure)]
  plotData1[,exposure := gsub("PE - ","",exposure)]
  setnames(plotData1,"exposure","Tissue / subgroup")
  setnames(plotData1,"Fstat2","F-stat")
  
  myTitle = paste0("PCSK9 levels on ",myOutcomes[i],", MR-ratio")
  if(myEstimate[i]=="OR") myXlab = "log(OR) per PCSK9 unit increase" 
  if(myEstimate[i]=="HR") myXlab = "log(HR) per PCSK9 unit increase"
  if(myEstimate[i]=="beta") myXlab = "beta per PCSK9 unit increase" 
  
  p1 <- forest(plotData1[,c(1,11,12,13)],
               est = plotData1$beta_IVW,
               lower = plotData1$lowerCI95, 
               upper = plotData1$upperCI95,
               sizes = 0.5,
               ci_column = 2,
               ref_line = 0,
               theme = tm1,
               xlab = myXlab,
               #xlim = c(-0.7,0.3),
               title = myTitle)
  
  plot(p1)
  
  filename = paste0("../results/ForestPlots/SupFig3_MRratio_",myOutcomes2[i],".png")
  png(filename = filename,width = 1900, height = 700, res=200)
  plot(p1)
  dev.off()
  
}

#' # Plots per outcome - MR-IVW
#' ***
for(i in 1:(length(myOutcomes)-1)){
  #i=1
  plotData1 = copy(MR_merged)
  plotData1 = plotData1[outcome == myOutcomes[i]]
  plotData1 = plotData1[method == "MR-IVW",]
  setorder(plotData1,method,exposure)
  plotData1[,Fstat2 := round(FStat,1)]
  plotData1[,Fstat2 := as.character(Fstat2)]
  
  dummy = data.table(exposure = c("Gene expression",
                                  "Protein expression"))

  plotData1 = rbind(plotData1,dummy, fill=T)
  plotData1 = plotData1[c(9,1:6,10,7,8)]
  plotData1[is.na(outcome),` ` := ""]
  plotData1[is.na(outcome),`Estimate [95% CI]`:= ""]
  plotData1[is.na(outcome),Fstat2 := ""]
  plotData1[!is.na(outcome),exposure := paste0("   ",exposure)]
  
  dummy2 = plotData1$exposure
  dummy2[c(1,8)] = "white"
  dummy2[grepl("GE",dummy2)] = myColor_GE[i]
  dummy2[grepl("PE",dummy2)] = myColor_PE[i]
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))
  
  plotData1[,exposure := gsub("GE - ","",exposure)]
  plotData1[,exposure := gsub("PE - ","",exposure)]
  setnames(plotData1,"exposure","Tissue / subgroup")
  setnames(plotData1,"Fstat2","F-stat")
  
  myTitle = paste0("PCSK9 levels on ",myOutcomes[i],", MR-IVW")
  if(myEstimate[i]=="OR") myXlab = "log(OR) per PCSK9 unit increase" 
  if(myEstimate[i]=="HR") myXlab = "log(HR) per PCSK9 unit increase"
  if(myEstimate[i]=="beta") myXlab = "beta per PCSK9 unit increase" 
  
  p1 <- forest(plotData1[,c(1,11,12,13)],
               est = plotData1$beta_IVW,
               lower = plotData1$lowerCI95, 
               upper = plotData1$upperCI95,
               sizes = 0.5,
               ci_column = 2,
               ref_line = 0,
               theme = tm1,
               xlab = myXlab,
               #xlim = c(-0.7,0.3),
               title = myTitle)
  
  plot(p1)
  
  filename = paste0("../results/ForestPlots/SupFig2_MRIVW_",myOutcomes2[i],".png")
  png(filename = filename,width = 1900, height = 800, res=200)
  plot(p1)
  dev.off()
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
