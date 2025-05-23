#' ---
#' title: "Create Main Figure: Forest Plots of MR-IVW"
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
#' - Forest Plot of the MR results (both IVW and ratio)
#' - For exposures PCSK9 protein levels and LDL-C
#' - For outcome Breast Cancer Survival 
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

#' # Filter ####
#' ***
MRTab[,unique(outcome)]
MRTab = MRTab[grepl("BCS",outcome),]
MRTab = MRTab[!grepl("gene expression",exposure),]
MRTab[grepl("BCAC",outcome), outcome := "BC Survival [BCAC]"]
MRTab[grepl("Mei ",outcome), outcome := "BC Survival [Mei et al.]"]

#' add nice name for phenotype
unique(MRTab$exposure)
MRTab[,exposure := gsub(" protein levels","",exposure)]
MRTab[,exposure := gsub("LDL-C levels","LDLC",exposure)]
MRTab[,exposure := gsub(" ratio","",exposure)]
MRTab[,exposure := gsub(" PCSK9 SNPs","",exposure)]
MRTab[,exposure := gsub("HMGCR SNPs","(HMGCR)",exposure)]
MRTab[,exposure := gsub("gw SNPs","(any)",exposure)]
MRTab[,exposure := gsub("free","all",exposure)]
MRTab[!grepl("[()]",exposure),exposure := paste(exposure, "(PCSK9)")]
MRTab[,exposure := paste(c(7,7,7,8,8,8,6,6,5,5,2,6,4,1,5,3),exposure)]

MRTab[,outcome2 := gsub("BC Survival","MR-IVW",outcome)]
MRTab[method=="MR-ratio",outcome2 := gsub("MR-IVW","MR-ratio",outcome2)]
MRTab[method=="MR-ratio" & !grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn<1,456\nrs562556")]
MRTab[method=="MR-ratio" & grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn=91,868\nrs562556")]
MRTab[method!="MR-ratio",outcome2 := paste0(outcome2,"\nn=91,868\nmultiple SNPs")]

#' Add columns necessary for plotting
MRTab[,lowerCI95 := beta-1.96*se]
MRTab[,upperCI95 := beta+1.96*se]

MRTab[grepl("PCSK9",exposure),colGroup := "PCSK9"]
MRTab[grepl("LDLC",exposure),colGroup := "LDLC"]

MRTab$signif <- ifelse(MRTab$pval_adj < 0.05,1,0)
MRTab[signif==1, y_ast:=upperCI95 + 2]
MRTab[signif==1 & method=="MR-IVW", y_ast:=-0.1]

MRTab[,sex := "females"]
MRTab[grepl("all",exposure),sex := "all"]

MRTab[,alphaGroup := 1]
MRTab[grepl("all",exposure),alphaGroup := 0.95]

MRTab[,colGroup2 := paste0(colGroup,sex,sep="::")]
MRTab[,exposure2 := gsub("PCSK9 - females [(]","",exposure)]
MRTab[,exposure2 := gsub("LDLC - females [(]","",exposure2)]
MRTab[,exposure2 := gsub("PCSK9 - all [(]","",exposure2)]
MRTab[,exposure2 := gsub("LDLC - all [(]","",exposure2)]
MRTab[,exposure2 := gsub("[)]","",exposure2)]

#' # Forest Plot ####
#' ***

p1 = ggplot(data=MRTab, aes(x=exposure2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=colGroup2)) +
  # Make range for ggplot values based on the data and AES specified in first line:
  geom_pointrange()+ 
  # add a dotted line at x=0 after flip:
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  # Makes whiskers on the range (more aesthetically pleasing):
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  # Makes DV header (Can handle multiple DVs):
  facet_wrap(~outcome2,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MRTab[MRTab$signif ==1, ],aes(x=exposure2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
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
        legend.position ="none", 
        axis.line.x = element_line(colour = "black"), 
        axis.line.y = element_blank(), 
        panel.border= element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(),
        axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black",face = "italic"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))

p1
filename = paste0("../results/Main_Sup_Figs/Main_ForestPlot_MR_BCS.png")
png(filename = filename,width = 2500, height = 1500, res=200)
plot(p1)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
