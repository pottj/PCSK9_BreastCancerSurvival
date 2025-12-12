#' ---
#' title: "Main Figures: Forest Plots of MR-ratio"
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

#' # Data prep ####
#' ***
plotData1 = copy(MRTab)
plotData1[,unique(exposure)]
plotData1[,unique(outcome)]
plotData1 = plotData1[grepl("BCS",outcome),]
plotData1 = plotData1[grepl("rs562556",exposure),]
plotData1 = plotData1[grepl("PE",exposure),]
setorder(plotData1,beta)

#' add columns necessary for plotting
plotData1[,sex := "sex-combined"]
plotData1[grepl("females",exposure),sex := "females"]
plotData1[,geneticModel := gsub(".* - ","",outcome)]
dummy = unlist(strsplit(plotData1$outcome," - "))
plotData1[,outcomeSource := dummy[seq(3,length(dummy),4)]]
plotData1[,lowerCI95 := beta-1.96*se]
plotData1[,upperCI95 := beta+1.96*se]

plotData1[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData1[signif==1, y_ast := ceiling(upperCI95)]

setorder(plotData1,outcomeSource,-geneticModel,sex,beta)
plotData1$ES <- paste(rep(" ", 30), collapse = " ")
plotData1$CI <- ifelse(is.na(plotData1$se), "",sprintf("%.2f [%.2f, %.2f]",plotData1$beta, plotData1$lowerCI95, plotData1$upperCI95))
plotData1[,FStat2 := round(FStat,1)]
plotData1[,FStat2 := as.character(FStat2)]
plotData1[,Outcome := gsub("BCS - females - ","",outcome)]
plotData1[,Outcome := gsub(" - .*","",Outcome)]

plotData1[,mean := beta]
plotData1[,lower := lowerCI95]
plotData1[,upper := upperCI95]

plotData2 = copy(plotData1)
plotData2 = plotData2[Outcome %in% c("Bertucci et al.","HMF Cohort","Nik-Zainal et al.","TCGA-BRCA","Mei et al."),]
plotData2[Outcome == "Mei et al.", Outcome := "  meta-analysis"]

p2a = plotData2 |>
  filter(Outcome != "  meta-analysis") |>
  forestplot(labeltext = c(Outcome, geneticModel, sex, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -10, to = 100, by = 10),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 levels",
             clip = c(-10,100),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_append_row(mean  = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$mean,
                lower = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$lower,
                upper = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$upper,
                Outcome = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$Outcome,
                geneticModel = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$geneticModel,
                sex = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$sex,
                CI = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$CI,
                FStat2 = filter(plotData2, Outcome == "  meta-analysis" & sex == "sex-combined")$FStat2,
                is.summary = TRUE) |>
  fp_append_row(mean  = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$mean,
                lower = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$lower,
                upper = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$upper,
                Outcome = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$Outcome,
                geneticModel = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$geneticModel,
                sex = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$sex,
                CI = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$CI,
                FStat2 = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$FStat2,
                is.summary = TRUE) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                sex = c("Exposure\nsex\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a)

filename = paste0("../results/MainFigures/MRratio_BCS_PE_Mei.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p2a)
dev.off()

plotData3 = copy(plotData1)
plotData3 = plotData3[!(Outcome %in% c("Bertucci et al.","HMF Cohort","Nik-Zainal et al.","TCGA-BRCA","Mei et al.")),]

p2b = plotData3 |>
  forestplot(labeltext = c(Outcome, geneticModel, sex, CI, FStat2),
             boxsize = 0.25,
             xticks = seq(from = -10, to = 15, by = 5),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 levels",
             clip = c(-10,15),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                sex = c("Exposure\nsex\n"),
                CI = c("log(HR) (95% CI)\n"),
                FStat2 = c("F-Stat\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b)

filename = paste0("../results/MainFigures/MRratio_BCS_PE_BCAC.png")
png(filename = filename,width = 2000, height = 700, res=200)
plot(p2b)
dev.off()

#' # Females only ####
#' ***
#' I want to use females only. That means I can remove the sex and F-stat column (because it is the same all the time). Instead I add the max sample size for these studies

plotData2[,maxSamples := rep(c("<553","<283","<1,456","<101","<519"),each=2)]

p2a2 = plotData2 |>
  filter(Outcome != "  meta-analysis") |>
  filter(sex == "females") |>
  forestplot(labeltext = c(Outcome, geneticModel,maxSamples, CI),
             boxsize = 0.25,
             xticks = seq(from = -10, to = 100, by = 10),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 levels",
             clip = c(-10,100),
             #xlog = TRUE,
             align = c("l","l","r","r"),
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_append_row(mean  = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$mean,
                lower = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$lower,
                upper = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$upper,
                Outcome = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$Outcome,
                geneticModel = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$geneticModel,
                maxSamples = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$maxSamples,
                CI = filter(plotData2, Outcome == "  meta-analysis" & sex != "sex-combined")$CI,
                is.summary = TRUE) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) (95% CI)\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2a2)

filename = paste0("../results/MainFigures/MRratio_BCS_PE_Mei_females.png")
png(filename = filename,width = 2000, height = 600, res=200)
plot(p2a2)
dev.off()

plotData3[,maxSamples := rep(c("4,648","4,648","91,686"),each=2)]

p2b2 = plotData3 |>
  filter(sex == "females") |>
  forestplot(labeltext = c(Outcome, geneticModel, maxSamples, CI),
             boxsize = 0.25,
             xticks = seq(from = -10, to = 15, by = 5),
             xlab = "logHR for BC Survival (95% CI) per 1-SD increment in PCSK9 levels",
             clip = c(-10,15),
             align = c("l","l","r","r"),
             #xlog = TRUE,
             vertices = TRUE,
             txt_gp = fpTxtGp(
               ticks = gpar(cex = 1),
               xlab = gpar(cex = 1))) |>
  fp_add_header(Outcome = c("Outcome\nsource\n"),
                geneticModel = c("Genetic\nmodel\n"),
                maxSamples = c("Sample\nsize\n"),
                CI = c("log(HR) (95% CI)\n")) |>
  fp_set_zebra_style("#EFEFEF") |>
  fp_add_lines() |>
  fp_decorate_graph(graph.pos = 4)

plot(p2b2)

filename = paste0("../results/MainFigures/MRratio_BCS_PE_BCAC_females.png")
png(filename = filename,width = 2000, height = 500, res=200)
plot(p2b2)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
