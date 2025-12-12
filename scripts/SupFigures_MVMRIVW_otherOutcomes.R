#' ---
#' title: "Main Figures: Forest Plots of MR-IVW"
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

#' # Load data ####
#' ***
load("../results/05_MVMR.RData")
MVMRTab[,pval_adj1 := p.adjust(p=pval_exp1,method = "fdr")]
MVMRTab[,pval_adj2 := p.adjust(p=pval_exp1,method = "fdr")]

#' # Plot 3: MR-IVW ####
#' ***
#' 
plotData1 = copy(MVMRTab)
plotData1[,unique(exposure1)]
plotData1[,unique(outcome)]
plotData1 = plotData1[!grepl("BCS",outcome) & !grepl("PLD",outcome),]

#' add columns necessary for plotting
dummy = unlist(strsplit(plotData1$outcome," - "))
plotData1[,Outcome := dummy[seq(1,length(dummy),4)]]
plotData1[,outcomeSource := dummy[seq(3,length(dummy),4)]]
plotData1[grepl("CAD - all",outcome),outcomeSource := paste(outcomeSource,"all")]
plotData1[grepl("CAD - fem",outcome),outcomeSource := paste(outcomeSource,"fem.")]

#' Plot for PCSK9
plotData2 = copy(plotData1)
plotData2[,lower := beta_exp1-1.96*SE_exp1]
plotData2[,upper := beta_exp1+1.96*SE_exp1]
plotData2[,mean := beta_exp1]

plotData2$ES <- paste(rep(" ", 30), collapse = " ")
plotData2$CI <- ifelse(is.na(plotData2$SE_exp1), "",sprintf("%.2f [%.2f, %.2f]",plotData2$beta_exp1, plotData2$lower, plotData2$upper))
plotData2[,FStat2 := round(condFstat_exp1,1)]
plotData2[,FStat2 := as.character(FStat2)]
plotData2[,SNPset := gsub("_"," & ",SNPset)]

setorder(plotData2,Outcome,SNPset,sex)

p2b = plotData2 |>
  forestplot(labeltext = c(sex, SNPset, Outcome, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             #xticks = seq(from = -4, to = 4, by = 2),
             xlab = "logOR (95% CI) for BC and CAD Risk \nper 1-SD increment in PCSK9 levels (corrected for LDL-C)",
             #clip = c(-4,4),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1),
               label = list(
                 grid::gpar(),
                 grid::gpar(fontface = "italic"),
                 grid::gpar(),
                 grid::gpar()))) |>
  fp_add_header(Outcome = c("Outcome \n\n"),
                outcomeSource = c("Outcome \nsource\n"),
                sex = c("Exposure\nsex\n"),
                SNPset = c("IV \nregions\n"),
                CI = c("log(OR) (95% CI)\n"),
                FStat2 = c("cond. \nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 5)

plot(p2b)

filename = paste0("../results/SupplementalFigures/MVMRIVW_otherOutcomes_PCSK9.png")
png(filename = filename,width = 2000, height = 900, res=200)
plot(p2b)
dev.off()

#' Same for LDL-C
#' 
plotData3 = copy(plotData1)
plotData3[,lower := beta_exp2-1.96*SE_exp2]
plotData3[,upper := beta_exp2+1.96*SE_exp2]
plotData3[,mean := beta_exp2]

plotData3$ES <- paste(rep(" ", 30), collapse = " ")
plotData3$CI <- ifelse(is.na(plotData3$SE_exp2), "",sprintf("%.2f [%.2f, %.2f]",plotData3$beta_exp2, plotData3$lower, plotData3$upper))
plotData3[,FStat2 := round(condFstat_exp2,1)]
plotData3[,FStat2 := as.character(FStat2)]
plotData3[,SNPset := gsub("_"," & ",SNPset)]

setorder(plotData3,Outcome,SNPset,sex)

p2c = plotData3 |>
  forestplot(labeltext = c(sex, SNPset, Outcome, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             #xticks = seq(from = -4, to = 4, by = 2),
             xlab = "logOR (95% CI) for BC and CAD Risk \nper 1-SD increment in LDL-C levels (corrected for PCSK9)",
             #clip = c(-4,4),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1),
               label = list(
                 grid::gpar(),
                 grid::gpar(fontface = "italic"),
                 grid::gpar(),
                 grid::gpar()))) |>
  fp_add_header(Outcome = c("Outcome \n\n"),
                outcomeSource = c("Outcome \nsource\n"),
                sex = c("Exposure\nsex\n"),
                SNPset = c("IV \nregions\n"),
                CI = c("log(OR) (95% CI)\n"),
                FStat2 = c("cond. \nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 5)

plot(p2c)

filename = paste0("../results/SupplementalFigures/MVMRIVW_otherOutcomes_LDLC.png")
png(filename = filename,width = 2000, height = 900, res=200)
plot(p2c)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
