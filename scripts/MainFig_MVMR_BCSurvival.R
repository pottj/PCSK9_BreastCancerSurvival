#' ---
#' title: "Create Main Figure: Forest Plots of MVMR-IVW"
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
#' What to show in the main figure? 
#' 
#' - Forest Plot of the MVMR results
#' - For exposures PCSK9 protein levels and LDL-C
#' - For outcome Breast Cancer Survival and CAD (females)
#' - Only MVMR in which the conditional F-Statistic is greater than 8
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

#' # Get MR results
#' ***
load("../results/05_MVMR.RData")

MVMRTab[,pval_adj1 := p.adjust(p=pval_exp1,method = "fdr")]
MVMRTab[,pval_adj2 := p.adjust(p=pval_exp2,method = "fdr")]

MVMRTab = MVMRTab[outcome %in% c("BCS - BCAC","CAD - females","CAD - all","BC - FinnGen + UKB")]
MVMRTab = MVMRTab[!(grepl("all",sex) & grepl("female",outcome)),]
MVMRTab = MVMRTab[!(grepl("female",sex) & grepl("all",outcome)),]
MVMRTab = MVMRTab[condFstat_exp1 >9 & condFstat_exp2 >9,]

MVMRTab[SNPset=="PCSK9_v1",SNPset := "PCSK9"]
MVMRTab[SNPset=="PCSK9_HMGCR",SNPset := "PCSK9 + HMGCR"]
MVMRTab[grepl("CAD",outcome),outcome := "CAD risk[Aragam et al.]"]
MVMRTab[grepl("BCS",outcome),outcome := "BC Survival [BCAC]"]
MVMRTab[grepl("FinnGen",outcome),outcome := "BC risk [FinnGen + UKB]"]

MVMRTab[,SNPset2 := paste(c(8,8,8,4,4,4,6,6,6,2,2,2),SNPset)]

MVMRTab1 = copy(MVMRTab)
MVMRTab1 = MVMRTab1[,c(1:4,5:9,17,19)]
MVMRTab2 = copy(MVMRTab)
MVMRTab2 = MVMRTab2[,c(1:4,10:14,18,19)]
setnames(MVMRTab1,"exposure1","exposure")
setnames(MVMRTab1,"pval_adj1","pval_adj")
names(MVMRTab1) = gsub("_exp1","",names(MVMRTab1))
setnames(MVMRTab2,"exposure2","exposure")
setnames(MVMRTab2,"pval_adj2","pval_adj")
names(MVMRTab2) = gsub("_exp2","",names(MVMRTab2))

MVMR = rbind(MVMRTab1,MVMRTab2)

#' Add columns necessary for plotting
MVMR[,lowerCI95 := beta-1.96*SE]
MVMR[,upperCI95 := beta+1.96*SE]

MVMR$signif <- ifelse(MVMR$pval_adj < 0.05,1,0)
MVMR[signif==1, y_ast:=upperCI95 + 0.4]
MVMR[signif==1 & grepl("CAD",outcome), y_ast:=upperCI95 + 0.2]

MVMR[,colGroup2 := paste0(exposure,sex,sep="::")]
# MVMR[exposure == "LDLC",SNPset2 := gsub("[2,4,6,8] ","",SNPset2)]
# MVMR[exposure == "LDLC",SNPset2 := paste(c(7,7,3,3,5,5,1,1),SNPset2)]

MVMR[exposure == "PCSK9", exposure := "exposure 1: PCSK9 levels"]
MVMR[exposure == "LDLC", exposure := "exposure 2: LDL-C levels"]

#' # Forest Plot for BC Survival ####
#' ***
p1 = ggplot(data=MVMR[grepl("BC S",outcome),], aes(x=SNPset2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=colGroup2)) +
  # Make range for ggplot values based on the data and AES specified in first line:
  geom_pointrange()+ 
  # add a dotted line at x=0 after flip:
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  # Makes whiskers on the range (more aesthetically pleasing):
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  # Makes DV header (Can handle multiple DVs):
  facet_wrap(~exposure,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MVMR[MVMR$signif ==1 & grepl("BC S",outcome), ],aes(x=SNPset2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logHR for breast cancer survival (95% CI) \nper 1-SD increment in exposure levels") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  scale_colour_manual(values = c("darkorange","darkorange4","dodgerblue","dodgerblue3")) +
  # My personal theme for GGplots:
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="gray90"),
        strip.text.x = element_text(size = 15),
        legend.position ="none", 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_blank(), 
        panel.border= element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.title=element_text(size=15),
        #axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black",face = "italic"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=10), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))
p1

filename = paste0("../results/Main_Sup_Figs/Main_ForestPlot_MVMR_BCS.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p1)
dev.off()

#' # Forest Plot for CAD in all ####
#' ***
p2 = ggplot(data=MVMR[grepl("CAD",outcome),], aes(x=SNPset2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=colGroup2)) +
  # Make range for ggplot values based on the data and AES specified in first line:
  geom_pointrange()+ 
  # add a dotted line at x=0 after flip:
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  # Makes whiskers on the range (more aesthetically pleasing):
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  # Makes DV header (Can handle multiple DVs):
  facet_wrap(~exposure,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MVMR[MVMR$signif ==1 & grepl("CAD",outcome), ],aes(x=SNPset2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logOR for CAD risk (95% CI) \nper 1-SD increment in exposure levels") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  scale_colour_manual(values = c("darkorange","darkorange4","dodgerblue","dodgerblue3")) +
  # My personal theme for GGplots:
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="gray90"),
        strip.text.x = element_text(size = 15),
        legend.position ="none", 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_blank(), 
        panel.border= element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.title=element_text(size=15),
        #axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black",face = "italic"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=10), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))
p2

filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MVMR_CAD.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p2)
dev.off()

#' # Forest Plot for BC ris in all ####
#' ***
p3 = ggplot(data=MVMR[grepl("BC risk",outcome),], aes(x=SNPset2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=colGroup2)) +
  # Make range for ggplot values based on the data and AES specified in first line:
  geom_pointrange()+ 
  # add a dotted line at x=0 after flip:
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  # Makes whiskers on the range (more aesthetically pleasing):
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  # Makes DV header (Can handle multiple DVs):
  facet_wrap(~exposure,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MVMR[MVMR$signif ==1 & grepl("BC risk",outcome), ],aes(x=SNPset2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logOR for BC risk (95% CI) \nper 1-SD increment in exposure levels") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  scale_colour_manual(values = c("darkorange","darkorange4","dodgerblue","dodgerblue3")) +
  # My personal theme for GGplots:
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="gray90"),
        strip.text.x = element_text(size = 15),
        legend.position ="none", 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_blank(), 
        panel.border= element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.title=element_text(size=15),
        #axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        #axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black",face = "italic"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=10), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))
p3

filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MVMR_BC.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p3)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
