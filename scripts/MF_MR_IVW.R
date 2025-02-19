#' ---
#' title: "Get main figures and main tables for MR-IVW"
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
#' Main Figure 1: Forest plot of MR-IVW results of PCSK9 on BCS
#' 
#' Main Figure 2: Forest plot of MR-IVW results of PCSK9 on BC
#' 
#' Main Table 1: Table with MR-IVW results per PCSK9 trait on BCS and BC
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
load("../results/05_MR_IVW.RData")

#' # Filter ####
#' ***
#' Only meta-results per outcome
MR_IVW = MR_IVW[grepl("meta",outcome),]

#' remove statin-free women
MR_IVW = MR_IVW[exposure!="PCSK9_females_free",]

#' add nice name for phenotype
MR_IVW$exposure
MR_IVW[grepl("PCSK9_",exposure),group := "Protein expression"]
MR_IVW[!grepl("PCSK9_",exposure),group := "Gene expression"]
MR_IVW[,exposure2 := c(rep("   Brain (CH)",2),
                       rep("   Brain (Cb)",2),
                       rep("   Liver",2),
                       rep("   Lung",2),
                       rep("   Nerve (tibial)",2),
                       rep("   Pancreas",2),
                       rep("   Women",2),
                       rep("   Statin-free",2))]

MR_IVW[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_IVW[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_IVW$` ` <- paste(rep(" ", 50), collapse = " ")
MR_IVW$`Estimate \n[95% CI]` <- ifelse(is.na(MR_IVW$se_IVW), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 MR_IVW$beta_IVW, MR_IVW$lowerCI95, MR_IVW$upperCI95))

#' # Plot 1: Breast Cancer Survival
#' ***
plotData1 = copy(MR_IVW)
plotData1 = plotData1[outcome == "BCS_meta",]
plotData1[,Fstat2 := round(FStat,1)]
plotData1[,Fstat2 := as.character(Fstat2)]

setorder(plotData1,group,-beta_IVW)

dummy = data.table(exposure2 = unique(plotData1$group))
plotData1 = rbind(plotData1,dummy, fill=T)
plotData1 = plotData1[c(9,1:6,10,7,8)]
plotData1[is.na(exposure),` ` := ""]
plotData1[is.na(exposure),`Estimate \n[95% CI]`:= ""]
plotData1[is.na(exposure),Fstat2 := ""]

dummy2 = plotData1$exposure
dummy2[is.na(dummy2)] = "white"
dummy2[!grepl("PCSK9",dummy2) & dummy2!="white"] = "#C2F1C8"
dummy2[grepl("PCSK9",dummy2)] = "#CAEEFB"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))

setnames(plotData1,"exposure2","Exposure subgroups")
setnames(plotData1,"Fstat2","F-stat")

p1 <- forest(plotData1[,c(9,12,13,14)],
            est = plotData1$beta_IVW,
            lower = plotData1$lowerCI95, 
            upper = plotData1$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            theme = tm1,
            xlab = "Causal effect estimate",
            title = "MR-IVW of PCSK9 levels on breast cancer survival")

plot(p1)

filename = paste0("../results/ForestPlots/MR_IVW_BreastCancerSurvival.png")
png(filename = filename,width = 1800, height = 800, res=200)
plot(p1)
dev.off()

#' # Plot 2: Breast Cancer
#' ***
plotData2 = copy(MR_IVW)
plotData2 = plotData2[outcome == "BCP_meta",]
plotData2[,Fstat2 := round(FStat,1)]
plotData2[,Fstat2 := as.character(Fstat2)]

setorder(plotData2,group,-beta_IVW)

dummy = data.table(exposure2 = unique(plotData2$group))
plotData2 = rbind(plotData2,dummy, fill=T)
plotData2 = plotData2[c(9,1:6,10,7,8)]
plotData2[is.na(exposure),` ` := ""]
plotData2[is.na(exposure),`Estimate \n[95% CI]`:= ""]
plotData2[is.na(exposure),Fstat2 := ""]

dummy2 = plotData2$exposure
dummy2[is.na(dummy2)] = "white"
dummy2[!grepl("PCSK9",dummy2) & dummy2!="white"] = "#C2F1C8"
dummy2[grepl("PCSK9",dummy2)] = "#CAEEFB"
tm2<- forest_theme(core=list(bg_params=list(fill = dummy2)))

setnames(plotData2,"exposure2","Exposure subgroups")
setnames(plotData2,"Fstat2","F-stat")

p2 <- forest(plotData2[,c(9,12,13,14)],
             est = plotData2$beta_IVW,
             lower = plotData2$lowerCI95, 
             upper = plotData2$upperCI95,
             sizes = 0.5,
             ci_column = 2,
             ref_line = 0,
             theme = tm2,
             xlab = "Causal effect estimate",
             title = "MR-IVW of PCSK9 levels on breast cancer")

plot(p2)

filename = paste0("../results/ForestPlots/MR_IVW_BreastCancer.png")
png(filename = filename,width = 1800, height = 800, res=200)
plot(p2)
dev.off()

#' # Table 1: Overview of MR-IVW results
#' ***
#' I want Exposure,	#IVs, F-stat, beta, SE, pval per outcome
setorder(plotData1,group,exposure)
setorder(plotData2,group,exposure)

tab1 = cbind(plotData1[3:10,c(9,3,14,4:6)],plotData2[3:10,4:6])
names(tab1)[4:6] = paste0("BCS_",names(tab1)[4:6]) 
names(tab1)[7:9] = paste0("BC_",names(tab1)[7:9]) 
names(tab1) = gsub("_IVW","",names(tab1))

tab1[,BCS_beta := round(BCS_beta,3)]
tab1[,BCS_se := round(BCS_se,3)]
tab1[,BC_beta := round(BC_beta,3)]
tab1[,BC_se := round(BC_se,3)]
tab1[,BCS_pval := signif(BCS_pval,3)]
tab1[,BC_pval := signif(BC_pval,3)]

tab1[1:6,`Exposure subgroups` := gsub("   ","GE - ",`Exposure subgroups`)]
tab1[7:8,`Exposure subgroups` := gsub("   ","PE - ",`Exposure subgroups`)]

tab1

WriteXLS(x = "tab1", 
         ExcelFileName=paste0("../results/MainTable1.xlsx"), 
         SheetNames="Table1", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
