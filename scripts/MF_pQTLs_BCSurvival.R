#' ---
#' title: "Get main figure"
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
#' Main Figure: Forest plot of MR-IVW and MR-RE results of PCSK9 pQTLs on BCS
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
MR_IVW = MR_IVW[grepl("Survival",outcome),]
MR_IVW = MR_IVW[grepl("PCSK9",exposure),]
MR_IVW[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]

MR_ratio[,unique(outcome)]
MR_ratio = MR_ratio[grepl("Survival",outcome),]
MR_ratio = MR_ratio[grepl("PCSK9",exposure),]
MR_ratio[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]
MR_ratio[outcome=="BCSurvival [Mei et al.]", outcome := "BC Survival [Mei et al.]"]

#' Merge data sets - what columns do I really need?
MR_ratio[,nSNPs := 1]
setnames(MR_ratio,"beta_ratio","beta_IVW")
setnames(MR_ratio,"se_ratio","se_IVW")
setnames(MR_ratio,"pval_ratio","pval_IVW")
setnames(MR_ratio,"Fstat","FStat")
MR_merged = rbind(MR_IVW[,1:7],MR_ratio[,c(7,13,23,19,20,21,22)])
MR_merged[,method := c(rep("MR-IVW",2),rep("MR-ratio",4))]

#' add nice name for phenotype
unique(MR_merged$exposure)
MR_merged[, exposure := gsub("PCSK9_","",exposure)]
MR_merged[exposure=="free", exposure := "statin-free individuals"]

MR_merged[,outcome2 := gsub("BC Survival","MR-IVW",outcome)]
MR_merged[method=="MR-ratio",outcome2 := gsub("MR-IVW","MR-ratio",outcome2)]

myColor_PE = c("#84E291","#47D45A","#8ED973")

MR_merged[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_merged[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_merged$` ` <- paste(rep(" ", 50), collapse = " ")
MR_merged$`Estimate [95% CI]` <- ifelse(is.na(MR_merged$se_IVW), "",
                                        sprintf("%.2f [%.2f, %.2f]",
                                                MR_merged$beta_IVW, MR_merged$lowerCI95, MR_merged$upperCI95))

#' # Plot 3: MR results for BC survival only
#' ***
plotData3 = copy(MR_merged)
plotData3[,Fstat2 := round(FStat,1)]
plotData3[,Fstat2 := as.character(Fstat2)]

dummy = data.table(exposure = unique(plotData3$outcome2))
plotData3 = rbind(plotData3,dummy, fill=T)
plotData3 = plotData3[c(7,1,2,8,3,4,9,5,6)]
plotData3[is.na(outcome),` ` := ""]
plotData3[is.na(outcome),`Estimate [95% CI]`:= ""]
plotData3[is.na(outcome),Fstat2 := ""]
plotData3[!is.na(outcome),exposure := paste0("   ",exposure)]

dummy2 = plotData3$exposure
dummy2[c(1,4,7)] = "white"
dummy2[c(2,3)] = myColor_PE[1]
dummy2[c(5,6)] = myColor_PE[2]
dummy2[c(8,9)] = myColor_PE[3]
tm3<- forest_theme(core=list(bg_params=list(fill = dummy2)))

setnames(plotData3,"exposure","Method [source] / subgroup")
setnames(plotData3,"Fstat2","F-stat")

p3 <- forest(plotData3[,c(1,12,13,14)],
             est = plotData3$beta_IVW,
             lower = plotData3$lowerCI95, 
             upper = plotData3$upperCI95,
             sizes = 0.5,
             ci_column = 2,
             ref_line = 0,
             theme = tm3,
             xlab = "log(HR) per PCSK9 unit increase",
             xlim = c(-1,23),
             title = "PCSK9 protein expression on breast cancer survival")

plot(p3)

filename = paste0("../results/ForestPlots/MainFig_pQTLs_BCSurvival.png")
png(filename = filename,width = 1900, height = 700, res=200)
plot(p3)
dev.off()


#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
