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
MVMRTab[,pval_adj1 := p.adjust(p=pval_exp2,method = "fdr")]

#' # Plot 3: MR-IVW ####
#' ***
#' 
plotData1 = copy(MVMRTab)
plotData1[,unique(exposure1)]
plotData1[,unique(outcome)]
plotData1 = plotData1[grepl("BCS",outcome),]

#' add columns necessary for plotting
plotData1[,outcomeSource := gsub("BCS - females - ","",outcome)]
plotData1[,outcomeSource := gsub(" - additive","",outcomeSource)]
plotData1[,sex := gsub("all","sex-combined",sex)]
plotData1[,lowerCI95 := beta_exp2-1.96*SE_exp2]
plotData1[,upperCI95 := beta_exp2+1.96*SE_exp2]

plotData1[,signif := ifelse(pval_adj1 < 0.05,1,0)]
plotData1[signif==1, y_ast := ceiling(upperCI95)]

plotData2 = copy(plotData1)
setorder(plotData2,outcomeSource,-SNPset,sex,beta_exp2)
plotData2$ES <- paste(rep(" ", 30), collapse = " ")
plotData2$CI <- ifelse(is.na(plotData2$SE_exp2), "",sprintf("%.2f [%.2f, %.2f]",plotData2$beta_exp2, plotData2$lowerCI95, plotData2$upperCI95))
plotData2[,FStat2 := round(condFstat_exp2,1)]
plotData2[,FStat2 := as.character(FStat2)]
plotData2[,Outcome := outcomeSource]
plotData2[,SNPset := gsub("_"," & ",SNPset)]

plotData2[,mean := beta_exp2]
plotData2[,lower := lowerCI95]
plotData2[,upper := upperCI95]

p2b = plotData2 |>
  forestplot(labeltext = c(sex, SNPset, Outcome,CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -4, to = 4, by = 2),
             xlab = "logHR (95% CI) for BC Survival \nper 1-SD increment in LDL-C levels (corrected for PCSK9)",
             clip = c(-4,4),
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
  fp_add_header(Outcome = c("Outcome \nsource\n"),
                sex = c("Exposure\nsex\n"),
                SNPset = c("IV \nregions\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("cond. \nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b)

filename = paste0("../results/SupplementalFigures//MVMRIVW_BCS_LDLC.png")
png(filename = filename,width = 2000, height = 800, res=200)
plot(p2b)
dev.off()

plotData3 = copy(plotData2)
plotData3 = plotData3[sex=="females",]
p3b = plotData3 |>
  forestplot(labeltext = c(Outcome, SNPset, CI,FStat2),
             boxsize = 0.25,
             xticks = seq(from = -4, to = 4, by = 2),
             xlab = "logHR for BC Survival (95% CI) \nper 1-SD increment in LDL-C levels (corrected for PCSK9)",
             clip = c(-4,4),
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
  fp_add_header(Outcome = c("Outcome \nsource\n"),
                SNPset = c("IV\nregion\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("cond. \nF-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 3)

plot(p3b)

filename = paste0("../results/SupplementalFigures/MVMRIVW_BCS_LDLC_females.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p3b)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
