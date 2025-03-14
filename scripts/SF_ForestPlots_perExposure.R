#' ---
#' title: "Sup Fig: Forest plot per exposure"
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

MR_merged[,exposure2 := "gene expression"]
MR_merged[grepl("PE",exposure),exposure2 := "protein expression"]
myExposures = unique(MR_merged$exposure2)
myExposures2 = c("eQTLs","pQTLs")
myColor_GE = c("#FBE3D6","#C2F1C8","#F2CFEE","#D9F2D0")
myColor_PE = c("#F2AA84","#47D45A","#D86ECC","#8ED973")

#' # Plots per outcome - MR-ratio
#' ***
for(i in 1:length(myExposures)){
  #i=2
  plotData1 = copy(MR_merged)
  plotData1 = plotData1[exposure2 == myExposures[i]]
  plotData1 = plotData1[method == "MR-ratio",]
  setorder(plotData1,outcome,exposure)
  plotData1[,Fstat2 := round(FStat,1)]
  plotData1[,Fstat2 := as.character(Fstat2)]
  
  dummy = data.table(exposure = unique(plotData1$outcome))
  
  plotData1 = rbind(plotData1,dummy, fill=T)
  if(i==1){
    plotData1 = plotData1[c(13,1:3,14,4:6,15,7:9,16,10:12)]
  }else{
    plotData1 = plotData1[c(9,1:2,10,3:4,11,5:6,12,7:8)]
  }
  plotData1[is.na(outcome),` ` := ""]
  plotData1[is.na(outcome),`Estimate [95% CI]`:= ""]
  plotData1[is.na(outcome),Fstat2 := ""]
  plotData1[!is.na(outcome),exposure := paste0("   ",exposure)]
  
  dummy2 = plotData1$exposure
  if(i==1){
    dummy2[c(1,5,9,13)] = "white"
    dummy2[2:4] = myColor_GE[2]
    dummy2[6:8] = myColor_GE[4]
    dummy2[10:12] = myColor_GE[1]
    dummy2[14:16] = myColor_GE[3]
  }else{
    dummy2[c(1,4,7,10)] = "white"
    dummy2[2:3] = myColor_PE[2]
    dummy2[5:6] = myColor_PE[4]
    dummy2[8:9] = myColor_PE[1]
    dummy2[11:12] = myColor_PE[3]
  }
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))
  
  plotData1[,exposure := gsub("GE - ","",exposure)]
  plotData1[,exposure := gsub("PE - ","",exposure)]
  
  if(i==1) setnames(plotData1,"exposure","Outcome [source] / tissue")
  if(i==2) setnames(plotData1,"exposure","Outcome [source] / subgroup")
  
  setnames(plotData1,"Fstat2","F-stat")
  
  myTitle = paste0("PCSK9 ",myExposures[i]," on outcomes, MR-ratio (rs562556)")

  if(i==1){
    p1 <- forest(plotData1[,c(1,11,12,14)],
                 est = plotData1$beta_IVW,
                 lower = plotData1$lowerCI95, 
                 upper = plotData1$upperCI95,
                 sizes = 0.5,
                 ci_column = 2,
                 ref_line = 0,
                 theme = tm1,
                 xlab = "Estimate per PCSK9 unit increase",
                 #xlim = c(-0.7,0.3),
                 title = myTitle)
  }else{
    p1 <- forest(plotData1[,c(1,11,12,14)],
                 est = plotData1$beta_IVW,
                 lower = plotData1$lowerCI95, 
                 upper = plotData1$upperCI95,
                 sizes = 0.5,
                 ci_column = 2,
                 ref_line = 0,
                 theme = tm1,
                 xlab = "Estimate per PCSK9 unit increase",
                 xlim = c(-2,23),
                 title = myTitle)
    
  }
  
  plot(p1)
  
  filename = paste0("../results/ForestPlots/SupFig5_MRratio_",myExposures2[i],".png")
  myHeight = 1100
  if(i==2) myHeight = 800
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p1)
  dev.off()
  
}

#' # Plots per outcome - MR-IVW
#' ***
for(i in 1:length(myExposures)){
  #i=2
  plotData1 = copy(MR_merged)
  plotData1 = plotData1[exposure2 == myExposures[i]]
  plotData1 = plotData1[method == "MR-IVW",]
  setorder(plotData1,method,outcome)
  plotData1[,Fstat2 := round(FStat,1)]
  plotData1[,Fstat2 := as.character(Fstat2)]
  
  dummy = data.table(exposure = unique(plotData1$outcome))
  
  plotData1 = rbind(plotData1,dummy, fill=T)
  if(i==1){
    plotData1 = plotData1[c(19,1:6,20,7:12,21,13:18)]
  }else{
    plotData1 = plotData1[c(7,1:2,8,3:4,9,5:6)]
  }
  
  plotData1[is.na(outcome),` ` := ""]
  plotData1[is.na(outcome),`Estimate [95% CI]`:= ""]
  plotData1[is.na(outcome),Fstat2 := ""]
  plotData1[!is.na(outcome),exposure := paste0("   ",exposure)]
  
  dummy2 = plotData1$exposure
  if(i==1){
    dummy2[c(1,8,15)] = "white"
    dummy2[2:7] = myColor_GE[2]
    dummy2[9:14] = myColor_GE[1]
    dummy2[16:21] = myColor_GE[3]
  }else{
    dummy2[c(1,4,7)] = "white"
    dummy2[2:3] = myColor_PE[2]
    dummy2[5:6] = myColor_PE[1]
    dummy2[8:9] = myColor_PE[3]
  }
  tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))
  
  plotData1[,exposure := gsub("GE - ","",exposure)]
  plotData1[,exposure := gsub("PE - ","",exposure)]
  if(i==1) setnames(plotData1,"exposure","Outcome [source] / tissue")
  if(i==2) setnames(plotData1,"exposure","Outcome [source] / subgroup")
  setnames(plotData1,"Fstat2","F-stat")
  
  myTitle = paste0("PCSK9 ",myExposures[i]," on outcomes, MR-IVW")
  
  p1 <- forest(plotData1[,c(1,11,12,14)],
               est = plotData1$beta_IVW,
               lower = plotData1$lowerCI95, 
               upper = plotData1$upperCI95,
               sizes = 0.5,
               ci_column = 2,
               ref_line = 0,
               theme = tm1,
               xlab = "Estimate per PCSK9 unit increase",
               #xlim = c(-0.7,0.3),
               title = myTitle)
  
  plot(p1)
  
  filename = paste0("../results/ForestPlots/SupFig4_MRIVW_",myExposures2[i],".png")
  myHeight = 1300
  if(i==2) myHeight = 800
  png(filename = filename,width = 1900, height = myHeight, res=200)
  plot(p1)
  dev.off()
  
}

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
