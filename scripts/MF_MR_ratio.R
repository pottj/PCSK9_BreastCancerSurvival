#' ---
#' title: "Get main figures and main tables for ratio estimates"
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
#' Main Figure 3: Forest plot of ratio estimates of PCSK9 on BCS in Morra et al.
#' 
#' Main Figure 4: Forest plot of ratio estimates of PCSK9 on BCS in Mei et al.
#' 
#' Main Table 2: Table with ratio estimates per PCSK9 trait on BCS in Morra et al. and Mei et al.
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
load("../results/06_MR_ratio_rs562556.RData")
MR_ratio_BCAC = copy(MR_ratio)

load("../results/07_MR_ratio_rs562556_Mei.RData")
MR_ratio_Mei = copy(MR_ratio)

MR_ratio_BCAC[,outcome_source := "Morra et al."]
MR_ratio_Mei[,outcome_source := "Mei et al."]

names(MR_ratio_BCAC)
names(MR_ratio_Mei)
names(MR_ratio_Mei)[13:15] = c("phenotype","beta","se" )

MR_ratio = rbind(MR_ratio_BCAC,MR_ratio_Mei,use.names=T,fill=T)
names(MR_ratio)[13:19] = paste0("outcome_",names(MR_ratio)[13:19])

#' # Filter ####
#' ***
#' Only meta-results per outcome
MR_ratio = MR_ratio[outcome_phenotype %in% c("BCS_meta","EUR only"),]

#' remove statin-free women
MR_ratio = MR_ratio[phenotype %in% c("Brain_Cerebellar_Hemisphere","Brain_Cerebellum","Colon_Transverse","Esophagus_Mucosa","Esophagus_Muscularis",
                                    "Liver","Lung","Nerve_Tibial","Pancreas","Small_Intestine_Terminal_Ileum","Spleen","PCSK9_females","PCSK9_free"),]

#' add nice name for phenotype
MR_ratio$phenotype
MR_ratio[grepl("PCSK9_",phenotype),group := "Protein expression"]
MR_ratio[!grepl("PCSK9_",phenotype),group := "Gene expression"]
MR_ratio[,exposure2 := rep(c("   Brain (CH)","   Brain (Cb)","   Colon (Transv)","   Esophagus (Muc)","   Esophagus (Musc)",
                             "   Liver","   Lung","   Nerve (tibial)","   Pancreas","   Small int. (TI)","   Spleen",
                             "   Women","   Statin-free"),2)]

MR_ratio[,lowerCI95 := beta_IV-1.96*se_IV2]
MR_ratio[,upperCI95 := beta_IV+1.96*se_IV2]

MR_ratio$` ` <- paste(rep(" ", 50), collapse = " ")
MR_ratio$`Estimate \n[95% CI]` <- ifelse(is.na(MR_ratio$se_IV2), "",
                                         sprintf("%.2f [%.2f, %.2f]",
                                                 MR_ratio$beta_IV, MR_ratio$lowerCI95, MR_ratio$upperCI95))

MR_ratio = MR_ratio[pval<0.05,]

#' # Plot 1: Protein expression 
#' ***
plotData1 = copy(MR_ratio)
plotData1[,Fstat2 := round(Fstat,1)]
plotData1[,Fstat2 := as.character(Fstat2)]
plotData1 = plotData1[group == "Protein expression"]

setorder(plotData1,outcome_source,-beta_IV)

dummy = data.table(exposure2 = unique(plotData1$outcome_source))
plotData1 = rbind(plotData1,dummy, fill=T)
plotData1 = plotData1[c(5,1,2,6,3,4)]
plotData1[is.na(phenotype),` ` := ""]
plotData1[is.na(phenotype),`Estimate \n[95% CI]`:= ""]
plotData1[is.na(phenotype),Fstat2 := ""]

dummy2 = plotData1$outcome_source
dummy2[is.na(dummy2)] = "white"
dummy2[grepl("Mei ",dummy2)] = "#F5C0A5"
dummy2[grepl("Morra ",dummy2)] = "#F4F6A4"
tm1<- forest_theme(core=list(bg_params=list(fill = dummy2)))

setnames(plotData1,"exposure2","Protein expression \n per sample stratification")
setnames(plotData1,"Fstat2","F-stat")

p1 <- forest(plotData1[,c(28,31,32,33)],
            est = plotData1$beta_IV,
            lower = plotData1$lowerCI95, 
            upper = plotData1$upperCI95,
            sizes = 0.5,
            ci_column = 2,
            ref_line = 0,
            theme = tm1,
            xlab = "Causal effect estimate",
            title = "MR ratio estimates of PCSK9 rs562556 on breast cancer survival")

plot(p1)

filename = paste0("../results/ForestPlots/MR_ratio_ProteinExpression.png")
png(filename = filename,width = 1800, height = 600, res=200)
plot(p1)
dev.off()

#' # Plot 2: Breast Cancer
#' ***
plotData2 = copy(MR_ratio)
plotData2[,Fstat2 := round(Fstat,1)]
plotData2[,Fstat2 := as.character(Fstat2)]
plotData2 = plotData2[group != "Protein expression"]

setorder(plotData2,outcome_source,-beta_IV)

dummy = data.table(exposure2 = unique(plotData2$outcome_source))
plotData2 = rbind(plotData2,dummy, fill=T)
plotData2 = plotData2[c(7,1:3,8,4:6)]
plotData2[is.na(phenotype),` ` := ""]
plotData2[is.na(phenotype),`Estimate \n[95% CI]`:= ""]
plotData2[is.na(phenotype),Fstat2 := ""]

dummy2 = plotData2$outcome_source
dummy2[is.na(dummy2)] = "white"
dummy2[grepl("Mei ",dummy2)] = "#F5C0A5"
dummy2[grepl("Morra ",dummy2)] = "#F4F6A4"
tm2<- forest_theme(core=list(bg_params=list(fill = dummy2)))

setnames(plotData2,"exposure2","Gene expression \n per tissue")
setnames(plotData2,"Fstat2","F-stat")

p2 <- forest(plotData2[,c(28,31,32,33)],
             est = plotData2$beta_IV,
             lower = plotData2$lowerCI95, 
             upper = plotData2$upperCI95,
             sizes = 0.5,
             ci_column = 2,
             ref_line = 0,
             theme = tm2,
             xlab = "Causal effect estimate",
             title = "MR ratio estimate of PCSK9 rs562556 on breast cancer survival")

plot(p2)

filename = paste0("../results/ForestPlots/MR_ratio_GeneExpression.png")
png(filename = filename,width = 1800, height = 700, res=200)
plot(p2)
dev.off()

#' # Table 1: Overview of MR-IVW results
#' ***
#' I want Exposure,	F-stat, beta, SE, pval per source
tab2 = cbind(MR_ratio[1:5,c(7,25,20,22,24)],MR_ratio[6:10,c(20,22,24)])
names(tab2)[3:5] = paste0("Morra_",names(tab2)[3:5]) 
names(tab2)[6:8] = paste0("Mei_",names(tab2)[6:8]) 
names(tab2) = gsub("_IV2","",names(tab2))
names(tab2) = gsub("_IV","",names(tab2))
names(tab2) = gsub("_p","_pval",names(tab2))

tab2[,phenotype := c("GE - Brain (Cb)", "GE - Esophagus (Musc)","GE - Spleen","PE - Women", "PE - Statin-free")]
tab2[,Fstat := round(Fstat,1)]

tab2[,Morra_beta := round(Morra_beta,3)]
tab2[,Morra_se := round(Morra_se,3)]
tab2[,Mei_beta := round(Mei_beta,3)]
tab2[,Mei_se := round(Mei_se,3)]
tab2[,Morra_pval := signif(Morra_pval,3)]
tab2[,Mei_pval := signif(Mei_pval,3)]

setnames(tab2,"phenotype","Exposure subgroups")
setnames(tab2,"Fstat","F-stat")

tab2

WriteXLS(x = "tab2", 
         ExcelFileName=paste0("../results/MainTable2.xlsx"), 
         SheetNames="Table2", 
         AutoFilter=T, 
         BoldHeaderRow=T,
         FreezeRow=1)

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
