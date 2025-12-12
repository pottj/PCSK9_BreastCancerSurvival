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
load("../results/04_MR.RData")
MRTab[,pval_adj := p.adjust(p=pval,method = "fdr")]

#' # Plot 2: MR-IVW ####
#' ***
#' 
plotData1 = copy(MRTab)
plotData1[,unique(exposure)]
plotData1[,unique(outcome)]
plotData1 = plotData1[grepl("BCS",outcome),]
plotData1 = plotData1[!grepl("rs562556",exposure),]

#' add columns necessary for plotting
dummy = unlist(strsplit(plotData1$exposure," - "))
plotData1[,setting := dummy[seq(2,length(dummy),4)]]
plotData1[,setting := gsub("all","sex-combined",setting)]
plotData1[,setting := gsub("Adipose Visceral Omentum","Adipose VO",setting)]
plotData1[,setting := gsub("Brain Cerebellar Hemisphere","Brain CH",setting)]
plotData1[,setting := gsub("Brain Cerebellum","Brain Cb",setting)]
plotData1[,type := dummy[seq(1,length(dummy),4)]]
plotData1[,gene := dummy[seq(3,length(dummy),4)]]

dummy = unlist(strsplit(plotData1$outcome," - "))
plotData1[,outcomeSource := dummy[seq(3,length(dummy),4)]]
plotData1[,lowerCI95 := beta-1.96*se]
plotData1[,upperCI95 := beta+1.96*se]

setorder(plotData1,setting,-type,gene,outcomeSource)
plotData1$ES <- paste(rep(" ", 30), collapse = " ")
plotData1$CI <- ifelse(is.na(plotData1$se), "",sprintf("%.2f [%.2f, %.2f]",plotData1$beta, plotData1$lowerCI95, plotData1$upperCI95))
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]

plotData1[,mean := beta]
plotData1[,lower := lowerCI95]
plotData1[,upper := upperCI95]

plotData2 = copy(plotData1)
plotData2 = plotData2[gene == "PCSK9" & !grepl("GE",type),]
plotData2[,type2 := gsub(" .*","",type)]
setorder(plotData2,outcomeSource,-type,setting)

p2a = plotData2 |>
  forestplot(labeltext = c(type2, setting, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -1.5, to = 1, by = 0.5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in PCSK9 or LDL-C levels",
             #clip = c(-1.7,1),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(type2 = c("Exposure\n\n"),
                outcomeSource = c("Outcome\nsource\n"),
                setting = c("Exposure\nsex\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a)

filename = paste0("../results/MainFigures/MRIVW_BCS_PE.png")
png(filename = filename,width = 2000, height = 700, res=200)
plot(p2a)
dev.off()

plotData2 = plotData2[setting == "females",]

p2b = plotData2 |>
  forestplot(labeltext = c(type2, outcomeSource, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -1.5, to = 1, by = 0.5),
             xlab = "logHR (95% CI) for BC Survival per 1-SD increment in PCSK9 or LDL-C levels",
             #clip = c(-1.7,1),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(type2 = c("Exposure\n\n"),
                outcomeSource = c("Outcome\nsource\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p2b)

filename = paste0("../results/MainFigures/MRIVW_BCS_PE_females.png")
png(filename = filename,width = 1800, height = 500, res=200)
plot(p2b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
