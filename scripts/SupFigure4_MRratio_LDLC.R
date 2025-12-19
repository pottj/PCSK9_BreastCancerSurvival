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
#' What kind of plot do I want: 
#' 
#' - Exposure: LDL-C levels in females
#' - Outcome: Breast Cancer Survival 
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

#' Filter data 
plotData1 = plotData1[grepl("Survival",outcome),]
plotData1 = plotData1[grepl("LDL",exposure),]
plotData1 = plotData1[grepl("fem",exposureSetting),]

#' Add columns necessary for plotting
plotData1[,mean := beta]
plotData1[,lower := beta-1.96*se]
plotData1[,upper := beta+1.96*se]
plotData1[,CI := ifelse(is.na(se), "",sprintf("%.2f [%.2f, %.2f]",beta, lower, upper))]

#' Set order: same as in Table 1: HMF, Bertucci, Nik-Zainal, TGCA, pooled, FinnGen recessive, FinnGen additive, Morra 
plotData1[,outcomeSource]
plotData1[,outcomeModel]
plotData1[,ranking := c(8,5,4,2,1,3,7,6)]
setorder(plotData1,ranking)

#' Add maximal sample sizes
plotData1[,maxSamples := c("<283","<553","<101","<519","<1,456","4,648","4,648","91,686")]

#' ## Plot A) Original studies
plotData2 = copy(plotData1)
plotData2 = plotData2[ranking<=5,]
#plotData2[Outcome == "Mei et al.", Outcome := "  meta-analysis"]

p2a = plotData2 |>
  filter(outcomeSource != "Mei et al. - pooled") |>
  forestplot(labeltext = c(outcomeSource, outcomeModel, maxSamples, CI),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 70, by = 10),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in LDL-C levels",
             clip = c(-1,70),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_append_row(mean  = filter(plotData2, outcomeSource == "Mei et al. - pooled")$mean,
                lower = filter(plotData2, outcomeSource == "Mei et al. - pooled")$lower,
                upper = filter(plotData2, outcomeSource == "Mei et al. - pooled")$upper,
                outcomeSource = filter(plotData2, outcomeSource == "Mei et al. - pooled")$outcomeSource,
                outcomeModel = filter(plotData2, outcomeSource == "Mei et al. - pooled")$outcomeModel,
                maxSamples = filter(plotData2, outcomeSource == "Mei et al. - pooled")$maxSamples,
                CI = filter(plotData2, outcomeSource == "Mei et al. - pooled")$CI,
                is.summary = TRUE) |>
  fp_add_header(outcomeSource = c("Outcome\nsource\n"),
                outcomeModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) [95% CI]\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a)

filename = paste0("../results/SupplementalFigures/SupFig4A_MRratio_LDLC.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2a)
dev.off()

#' ## Plot B) Replication cohorts
plotData3 = copy(plotData1)
plotData3 = plotData3[ranking>5,]

p2b = plotData3 |>
  forestplot(labeltext = c(outcomeSource, outcomeModel, maxSamples, CI),
             boxsize = 0.25,
             xticks = seq(from = -5, to = 10, by = 2.5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in LDL-C levels",
             clip = c(-6,11),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(outcomeSource = c("Outcome\nsource\n"),
                outcomeModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) [95% CI]\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b)

filename = paste0("../results/SupplementalFigures/SupFig4B_MRratio_LDLC.png")
png(filename = filename,width = 2000, height = 500, res=200)
plot(p2b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
