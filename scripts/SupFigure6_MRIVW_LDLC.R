#' ---
#' title: "Supplemental Figures: Forest Plots of MR-IVW"
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
#' What kind of plot do I want: 
#' 
#' - Exposures: LDL-C levels in females and sex-combined with different instrument selection
#' - Outcome: Breast Cancer Survival in FinnGen and Morra et al.
#' 
#' Necessary columns: 
#' 
#' - Gene region
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
plotData1 = copy(stab5)

#' Filter data 
plotData1 = plotData1[grepl("Survival",outcome),]
plotData1 = plotData1[grepl("LDL",exposure),]
plotData1 = plotData1[exposureInstruments != "PCSK9 and HMGCR",]

#' Add columns necessary for plotting
plotData1[,mean := beta]
plotData1[,lower := beta-1.96*se]
plotData1[,upper := beta+1.96*se]
plotData1[,CI := ifelse(is.na(se), "",sprintf("%.2f [%.2f, %.2f]",beta, lower, upper))]
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]
plotData1[exposureSetting=="all",exposureSetting := "sex-combined"]

#' Set order: outcome source, exposure sex, exposure
setorder(plotData1,-exposureInstruments,outcomeSource, exposureSetting)

#' ## Plot 
p3 = plotData1 |>
  forestplot(labeltext = c(exposureSetting, exposureInstruments, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -1.5, to = 1, by = 0.5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment \nin LDL-C levels (total effects)",
             clip = c(-1.5,1),
             #xlog = TRUE,
             align = c("l","l","l","r","r"),
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1),
               label = list(
                 grid::gpar(),
                 grid::gpar(fontface = "italic"),
                 grid::gpar(),
                 grid::gpar(),
                 grid::gpar(),
                 grid::gpar()))) |>
  fp_add_header(exposureInstruments = c("Instrument\nselection\n"),
                exposureSetting = c("Exposure\nsex\n"),
                outcomeSource = c("Outcome\nsource\n"),
                CI = c("\nlog(HR) [95% CI]\n"),
                FStat2 = c("\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p3)

filename = paste0("../results/SupplementalFigures/SupFig6_MRIVW_LDLC.png")
png(filename = filename,width = 2000, height = 900, res=200)
plot(p3)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
