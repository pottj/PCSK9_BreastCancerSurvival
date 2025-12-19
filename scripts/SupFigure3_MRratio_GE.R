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
#' - Exposure: PCSK9 gene expression levels in two tissues
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
plotData1 = plotData1[grepl("GE",exposure),]

#' Add columns necessary for plotting
plotData1[,mean := beta]
plotData1[,lower := beta-1.96*se]
plotData1[,upper := beta+1.96*se]
plotData1[,CI := ifelse(is.na(se), "",sprintf("%.2f [%.2f, %.2f]",beta, lower, upper))]

#' Set order: same as in Table 1: HMF, Bertucci, Nik-Zainal, TGCA, pooled, FinnGen recessive, FinnGen additive, Morra 
plotData1[,outcomeSource]
plotData1[,outcomeModel]
plotData1[,ranking := rep(c(8,5,4,2,1,3,7,6),2)]
setorder(plotData1,ranking)

#' Add maximal sample sizes
plotData1[,maxSamples := rep(c("<283","<553","<101","<519","<1,456","4,648","4,648","91,686"),each=2)]

#' ## Plot A) Original studies
plotData2 = copy(plotData1)
plotData2 = plotData2[ranking<=5,]

p2a = plotData2 |>
  filter(outcomeSource != "Mei et al. - pooled") |>
  forestplot(labeltext = c(outcomeSource, exposureSetting, maxSamples, CI),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 15, by = 2.5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in PCSK9 GE levels",
             #clip = c(-10,100),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_append_row(mean  = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$mean,
                lower = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$lower,
                upper = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$upper,
                outcomeSource = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$outcomeSource,
                exposureSetting = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$exposureSetting,
                maxSamples = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$maxSamples,
                CI = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting == "Spleen")$CI,
                is.summary = TRUE) |>
  fp_append_row(mean  = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$mean,
                lower = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$lower,
                upper = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$upper,
                outcomeSource = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$outcomeSource,
                exposureSetting = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$exposureSetting,
                maxSamples = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$maxSamples,
                CI = filter(plotData2, outcomeSource == "Mei et al. - pooled" & exposureSetting != "Spleen")$CI,
                is.summary = TRUE) |>
  fp_add_header(outcomeSource = c("Outcome\nsource\n"),
                exposureSetting = c("Exposure\ntissue\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) [95% CI]\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a)

filename = paste0("../results/SupplementalFigures/SupFig3A_MRratio_GE.png")
png(filename = filename,width = 2000, height = 900, res=200)
plot(p2a)
dev.off()

#' ## Plot B) Replication cohorts
plotData3 = copy(plotData1)
plotData3 = plotData3[ranking>5,]

p2b = plotData3 |>
  forestplot(labeltext = c(outcomeSource, exposureSetting, maxSamples, CI),
             boxsize = 0.25,
             #xticks = seq(from = -10, to = 15, by = 5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in PCSK9 GE levels",
             #clip = c(-10,15),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(outcomeSource = c("Outcome\nsource\n"),
                exposureSetting = c("Exposure\ntissue\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) [95% CI]\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b)

filename = paste0("../results/SupplementalFigures/SupFig3B_MRratio_GE.png")
png(filename = filename,width = 2000, height = 700, res=200)
plot(p2b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
