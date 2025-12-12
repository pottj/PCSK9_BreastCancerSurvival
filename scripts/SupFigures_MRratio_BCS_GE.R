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
#' PCSK9 gene expression on BCS using rs562556
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
plotData1 = plotData1[grepl("BCS",outcome),]
plotData1 = plotData1[grepl("rs562556",exposure),]
plotData1 = plotData1[grepl("GE",exposure),]
setorder(plotData1,beta)

#' add columns necessary for plotting
plotData1[,tissue := gsub("PCSK9 GE levels - ","",exposure)]
plotData1[,tissue := gsub(" - .*","",tissue)]
plotData1[,geneticModel := gsub(".* - ","",outcome)]
dummy = unlist(strsplit(plotData1$outcome," - "))
plotData1[,outcomeSource := dummy[seq(3,length(dummy),4)]]
plotData1[,lowerCI95 := beta-1.96*se]
plotData1[,upperCI95 := beta+1.96*se]

plotData1[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData1[signif==1, y_ast := ceiling(upperCI95)]

setorder(plotData1,outcomeSource,-geneticModel,tissue,beta)
plotData1$ES <- paste(rep(" ", 30), collapse = " ")
plotData1$CI <- ifelse(is.na(plotData1$se), "",sprintf("%.2f [%.2f, %.2f]",plotData1$beta, plotData1$lowerCI95, plotData1$upperCI95))
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]
plotData1[,Outcome := gsub("BCS - females - ","",outcome)]
plotData1[,Outcome := gsub(" - .*","",Outcome)]

plotData1[,mean := beta]
plotData1[,lower := lowerCI95]
plotData1[,upper := upperCI95]
plotData1[tissue!="Spleen",tissue := "Esophagus musc."]

plotData2 = copy(plotData1)
plotData2 = plotData2[Outcome %in% c("Bertucci et al.","HMF Cohort","Nik-Zainal et al.","TCGA-BRCA","Mei et al."),]
plotData2[Outcome == "Mei et al.", Outcome := "  meta-analysis"]
plotData2[,maxSamples := rep(c("<553","<283","<1,456","<101","<519"),each=2)]

p2a = plotData2 |>
  filter(Outcome != "  meta-analysis") |>
  forestplot(labeltext = c(Outcome, geneticModel, maxSamples, tissue, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = 0, to = 10, by = 2),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 GE levels",
             clip = c(-1,10),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_append_row(mean  = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$mean,
                lower = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$lower,
                upper = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$upper,
                Outcome = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$Outcome,
                geneticModel = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$geneticModel,
                maxSamples = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$maxSamples,
                tissue = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$tissue,
                CI = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$CI,
                FStat2 = filter(plotData2, Outcome == "  meta-analysis" & tissue != "Spleen")$FStat2,
                is.summary = TRUE) |>
  fp_append_row(mean  = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$mean,
                lower = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$lower,
                upper = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$upper,
                Outcome = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$Outcome,
                geneticModel = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$geneticModel,
                maxSamples = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$maxSamples,
                tissue = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$tissue,
                CI = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$CI,
                FStat2 = filter(plotData2, Outcome == "  meta-analysis" & tissue == "Spleen")$FStat2,
                is.summary = TRUE) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                tissue = c("Exposure\ntissue\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 5)

plot(p2a)

filename = paste0("../results/SupplementalFigures/MRratio_BCS_GE_Mei.png")
png(filename = filename,width = 2200, height = 1000, res=200)
plot(p2a)
dev.off()

plotData3 = copy(plotData1)
plotData3 = plotData3[!(Outcome %in% c("Bertucci et al.","HMF Cohort","Nik-Zainal et al.","TCGA-BRCA","Mei et al.")),]
plotData3[,maxSamples := rep(c("4,648","4,648","91,686"),each=2)]

p2b = plotData3 |>
  forestplot(labeltext = c(Outcome, geneticModel, maxSamples, tissue, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -2, to = 2, by = 1),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 GE levels",
             clip = c(-2,2),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                tissue = c("Exposure\ntissue\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 5)

plot(p2b)

filename = paste0("../results/SupplementalFigures/MRratio_BCS_GE_BCAC.png")
png(filename = filename,width = 2200, height = 700, res=200)
plot(p2b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
