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
#' What kind of plot do I want: 
#' 
#' - Exposures: PCSK9 gene expression levels in the various tissues
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
plotData1 = copy(stab5)

#' Filter data 
plotData1 = plotData1[grepl("Survival",outcome),]
plotData1 = plotData1[grepl("GE",exposure),]
plotData1 = plotData1[exposureInstruments == "PCSK9",]

#' Add columns necessary for plotting
plotData1[,mean := beta]
plotData1[,lower := beta-1.96*se]
plotData1[,upper := beta+1.96*se]
plotData1[,CI := ifelse(is.na(se), "",sprintf("%.2f [%.2f, %.2f]",beta, lower, upper))]
plotData1[,exposure := gsub(" levels","",exposure)]
plotData1[,exposure := gsub(" GE","",exposure)]
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]

#' Set order: outcome source, exposure sex, exposure
setorder(plotData1, exposureSetting,outcomeSource)

#' ## Plot 
p3 = plotData1 |>
  forestplot(labeltext = c(exposureSetting, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -0.75, to = 0.75, by = 0.25),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment \nin PCSK9 gene expression levels (total effects)",
             #clip = c(-1.5,1),
             #xlog = TRUE,
             align = c("l","l","r","r"),
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(exposureSetting = c("Exposure\ntissue\n"),
                outcomeSource = c("Outcome\nsource\n"),
                CI = c("\nlog(HR) [95% CI]\n"),
                FStat2 = c("\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p3)

filename = paste0("../results/SupplementalFigures/SupFig5_MRIVW_GE.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p3)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
