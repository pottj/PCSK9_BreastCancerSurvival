#' ---
#' title: "Sup Figures: Forest Plots of MVMR-IVW"
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
#' What kind of plot do I want: same structure as for Figure 3, but now with the direct independent effects 
#' 
#' - Exposures: PCSK9 protein levels and LDL-C levels in females and sex-combined
#' - Outcome: Breast Cancer Survival in FinnGen and Morra et al.
#' 
#' Necessary columns: 
#' 
#' - Exposure
#' - Exposure sex
#' - Outcome Source (cohort name or publication first author)
#' - estimate + CI
#' - F-Statistic
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

load("../results/SupTables.RData")

#' # Data prep ####
#' ***
plotData1 = copy(stab6)

#' Filter data 
plotData1 = plotData1[grepl("Coronary",outcome),]
plotData1 = plotData1[grepl("all",outcomeSex),]

#' Get long format 
plotData1a=copy(plotData1)
plotData1b=copy(plotData1)
plotData1b[,exposure1 := exposure2]
plotData1b[,beta_exp1 := beta_exp2]
plotData1b[,se_exp1 := se_exp2]
plotData1b[,pval_exp1 := pval_exp2]
plotData1b[,pval_adj_exp1 := pval_adj_exp2]
plotData1b[,condFStat_exp1 := condFStat_exp2]

plotData1 = rbind(plotData1a,plotData1b)

#' Add columns necessary for plotting
plotData1[,mean := beta_exp1]
plotData1[,lower := beta_exp1-1.96*se_exp1]
plotData1[,upper := beta_exp1+1.96*se_exp1]
plotData1[,CI := ifelse(is.na(se_exp1), "",sprintf("%.2f [%.2f, %.2f]",beta_exp1, lower, upper))]
plotData1[,exposure1 := gsub("LDLC","LDL-C",exposure1)]
plotData1[,FStat2 := round(condFStat_exp1,1)]
plotData1[,FStat2 := as.character(FStat2)]

#' Set order: outcome source, exposure sex, exposure
setorder(plotData1,outcomeSource, -exposureSex,-exposure1)

#' ## Plot 
p4a = plotData1 |>
  filter(instrumentLocus == "PCSK9") |>
  forestplot(labeltext = c(exposure1, exposureSex, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -0.5, to = 1, by = 0.25),
             xlab = "logOR (95% CI) for CAD risk per 1-SD increment \nin PCSK9 or LDL-C levels (direct effects)",
             #clip = c(-3,3),
             #xlog = TRUE,
             align = c("l","l","l","r","r"),
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(exposure1 = c("\nExposure\n"),
                exposureSex = c("Exposure\nsex\n"),
                CI = c("\nlog(OR) [95% CI]\n"),
                FStat2 = c("cond.\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p4a)

filename = paste0("../results/SupplementalFigures/SupFig11A_MVMRIVW_CAD.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p4a)
dev.off()

#' ## Plot 
p4b = plotData1 |>
  filter(instrumentLocus != "PCSK9") |>
  forestplot(labeltext = c(exposure1, exposureSex, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 0.6, by = 0.1),
             xlab = "logOR (95% CI) for CAD risk per 1-SD increment \nin PCSK9 or LDL-C levels (direct effects)",
             #clip = c(-3,3),
             #xlog = TRUE,
             align = c("l","l","l","r","r"),
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(exposure1 = c("\nExposure\n"),
                exposureSex = c("Exposure\nsex\n"),
                CI = c("\nlog(OR) [95% CI]\n"),
                FStat2 = c("cond.\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p4b)

filename = paste0("../results/SupplementalFigures/SupFig11B_MVMRIVW_CAD.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p4b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
