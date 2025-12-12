#' ---
#' title: "Supplemental Figures: Forest Plots of MR-ratio"
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
#' All other outcomes
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

#' # Load data ####
#' ***
load("../results/04_MR.RData")
MRTab[,pval_adj := p.adjust(p=pval,method = "fdr")]

#' # Data prep ####
#' ***
plotData1 = copy(MRTab)
plotData1[,unique(exposure)]
plotData1[,unique(outcome)]
plotData1 = plotData1[!grepl("BCS",outcome),]
plotData1 = plotData1[!grepl("PLD",outcome),]
plotData1 = plotData1[grepl("rs562556",exposure),]

#' add columns necessary for plotting
dummy = unlist(strsplit(plotData1$exposure," - "))
plotData1[,setting := dummy[seq(2,length(dummy),4)]]
plotData1[,setting := gsub("Muscularis","Musc.",setting)]
dummy = unlist(strsplit(plotData1$outcome," - "))
plotData1[,outcomeSource := dummy[seq(3,length(dummy),4)]]
plotData1[,Outcome := dummy[seq(1,length(dummy),4)]]
plotData1[,lowerCI95 := beta-1.96*se]
plotData1[,upperCI95 := beta+1.96*se]

plotData1[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData1[signif==1, y_ast := ceiling(upperCI95)]

setorder(plotData1,outcomeSource,setting,beta)
plotData1$ES <- paste(rep(" ", 30), collapse = " ")
plotData1$CI <- ifelse(is.na(plotData1$se), "",sprintf("%.2f [%.2f, %.2f]",plotData1$beta, plotData1$lowerCI95, plotData1$upperCI95))
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]

plotData1[,mean := beta]
plotData1[,lower := lowerCI95]
plotData1[,upper := upperCI95]
setorder(plotData1,outcome,exposure)
plotData1[,maxSamples := rep(c("419433","1165700","333473"),each=6)]
plotData1[,maxCases := rep(c("30593","181522","NA"),each=6)]
plotData1[7:12,outcomeSource := paste(outcomeSource,"all")]
plotData1[13:18,outcomeSource := paste(outcomeSource,"fem.")]

plotData2 = copy(plotData1)
plotData2 = plotData2[grepl("PE",exposure),]

p2a = plotData2 |>
  forestplot(labeltext = c(Outcome, outcomeSource, setting, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -2, to = 2, by = 0.5),
             xlab = "logOR (95% CI) for BC and CAD Risk per 1-SD increment in PCSK9 levels",
             clip = c(-2,2),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\n\n"),
                outcomeSource = c("Outcome\nsource\n"),
                setting = c("Exposure\nsex\n"),
                CI = c("log(OR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a)

filename = paste0("../results/SupplementalFigures/MRratio_otherOutcomes_PE.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2a)
dev.off()

plotData3 = copy(plotData1)
plotData3 = plotData3[grepl("GE",exposure),]

p2b = plotData3 |>
  forestplot(labeltext = c(Outcome, outcomeSource, setting, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -0.2, to = 0.2, by = 0.05),
             xlab = "logOR (95% CI) for BC and CAD Risk per 1-SD increment in PCSK9 GE levels",
             clip = c(-0.2,0.2),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\n\n"),
                outcomeSource = c("Outcome\nsource\n"),
                setting = c("Exposure\ntissue\n"),
                CI = c("log(OR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b)

filename = paste0("../results/SupplementalFigures/MRratio_otherOutcomes_GE.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2b)
dev.off()

plotData4 = copy(plotData1)
plotData4 = plotData4[grepl("LDL",exposure),]

p2c = plotData4 |>
  forestplot(labeltext = c(Outcome, outcomeSource, setting, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -1, to = 1, by = 0.25),
             xlab = "logOR (95% CI) for BC and CAD Risk per 1-SD increment in LDL-C levels",
             clip = c(-1.1,1.1),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\n\n"),
                outcomeSource = c("Outcome\nsource\n"),
                setting = c("Exposure\nsex\n"),
                CI = c("log(OR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2c)

filename = paste0("../results/SupplementalFigures/MRratio_otherOutcomes_LDLC.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2c)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
