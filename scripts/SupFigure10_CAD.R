#' ---
#' title: "Supp Figures: Forest Plots of other outcome"
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
#' - Exposure: PCSK9 protein levels in females
#' - Outcome: CAD
#' 
#' A) Original + Pooled analysis: data from Mei et al, in the MR-ratio approach
#' B) Replication: FinnGen & Morra et al. data in the MR-ratio approach
#' 
#' Necessary columns: 
#' 
#' - Outcome Source (cohort name or publication first author)
#' - Genetic model
#' - Sample size
#' - estimate + CI
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
plotData1 = copy(stab4)
plotData2 = copy(stab5)

#' Filter data 
plotData1 = plotData1[outcome == "Coronary Artery Disease",]
plotData1 = plotData1[outcomeSex == "all",]
plotData1 = plotData1[!grepl("GE",exposure),]
plotData1[,method := "MR-ratio"]
plotData2 = plotData2[outcome == "Coronary Artery Disease",]
plotData2 = plotData2[outcomeSex == "all",]
plotData2 = plotData2[exposureInstruments == "PCSK9",]
plotData2 = plotData2[!grepl("GE",exposure),]
plotData2[,method := "MR-IVW"]
plotData3 = rbind(plotData1,plotData2,fill=T)

#' Add columns necessary for plotting
plotData3[,mean := beta]
plotData3[,lower := beta-1.96*se]
plotData3[,upper := beta+1.96*se]
plotData3[,CI := ifelse(is.na(se), "",sprintf("%.2f [%.2f, %.2f]",beta, lower, upper))]
plotData3[,FStat2 := round(FStat,1)]
plotData3[,FStat2 := as.character(FStat2)]
plotData3[,exposure := gsub(" levels","",exposure)]
plotData3[,exposure := gsub(" PE","",exposure)]
plotData3[exposureSetting=="all",exposureSetting := "sex-combined"]

setorder(plotData3,-exposure,exposureSetting)

#' ## Plot A) single variant
p2a = plotData3 |>
  filter(method == "MR-ratio") |>
  forestplot(labeltext = c(exposure, exposureSetting, CI,FStat2),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 2, by = 0.25),
             xlab = "logOR (95% CI) for CAD risk per 1-SD increment in PCSK9 or LDL-C levels",
             #clip = c(-10,100),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(exposure = c("\nExposure\n"),
                exposureSetting = c("Exposure\nsex\n"),
                CI = c("\nlog(OR) [95% CI]\n"),
                FStat2 = c("\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p2a)

filename = paste0("../results/SupplementalFigures/SupFig10A_MRratio_CAD.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2a)
dev.off()

#' ## Plot B) multiple variants
p2b = plotData3 |>
  filter(method != "MR-ratio") |>
  forestplot(labeltext = c(exposure, exposureSetting, CI,FStat2),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 1, by = 0.25),
             xlab = "logOR (95% CI) for CAD risk per 1-SD increment in PCSK9 or LDL-C levels",
             #clip = c(-10,100),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(exposure = c("\nExposure\n"),
                exposureSetting = c("Exposure\nsex\n"),
                CI = c("\nlog(OR) [95% CI]\n"),
                FStat2 = c("\nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p2b)

filename = paste0("../results/SupplementalFigures/SupFig10B_MRIVW_CAD.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
