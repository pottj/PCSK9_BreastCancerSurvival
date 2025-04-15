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
#' In this script, I load the supplemental tables and restrict S2 and S3 to the BCS outcome. 
#' 
#' Then I create a forest plot with facets for each outcome source (Mei vs Morra) and MR method (IVW vs ratio). 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
load("../results/SupTables.RData")

MR_IVW = copy(stab2)
MR_ratio = copy(stab3)

#' # Filter ####
#' ***
MR_IVW[,unique(outcome)]
MR_IVW = MR_IVW[grepl("Survival",outcome),]
MR_IVW[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]

MR_ratio[,unique(outcome)]
MR_ratio = MR_ratio[grepl("Survival",outcome),]
MR_ratio = MR_ratio[exposure_pval<0.05,]
MR_ratio[outcome=="BCSurvival [BCAC]", outcome := "BC Survival [BCAC]"]
MR_ratio[outcome=="BCSurvival [Mei et al.]", outcome := "BC Survival [Mei et al.]"]

#' Merge data sets - what columns do I really need?
MR_ratio[,nSNPs := 1]
setnames(MR_ratio,"beta_ratio","beta_IVW")
setnames(MR_ratio,"se_ratio","se_IVW")
setnames(MR_ratio,"pval_ratio","pval_IVW")
setnames(MR_ratio,"Fstat","FStat")
MR_merged = rbind(MR_IVW[,1:7],MR_ratio[,c(7,13,24,20,21,22,23)])
MR_merged[,method := c(rep("MR-IVW",8),rep("MR-ratio",22))]

#' add nice name for phenotype
unique(MR_merged$exposure)
MR_merged[exposure=="Adipose_Visceral_Omentum", exposure:="Adipose VO"]
MR_merged[exposure=="Brain_Cerebellum", exposure:="Brain Cb"]
MR_merged[exposure=="Brain_Cerebellar_Hemisphere", exposure:="Brain CH"]
MR_merged[exposure=="Nerve_Tibial", exposure:="Nerve tib"]
MR_merged[exposure=="Whole_Blood", exposure:="Whole blood"]
MR_merged[exposure=="Esophagus_Muscularis", exposure:="Esophagus musc"]
MR_merged[exposure=="Artery_Aorta", exposure:="Artery aorta"]
MR_merged[exposure=="Esophagus_Gastroesophageal_Junction", exposure:="Esophagus GJ"]
MR_merged[exposure=="Skin_Not_Sun_Exposed_Suprapubic", exposure:="Skin sun not exp"]
MR_merged[exposure=="Skin_Sun_Exposed_Lower_leg", exposure:="Skin sun exp"]

MR_merged[!grepl("PCSK9",exposure), exposure := paste0("GE - ",exposure)]
MR_merged[grepl("PCSK9",exposure), exposure := gsub("PCSK9_","PE - ",exposure)]
MR_merged[exposure=="PE - free", exposure := "PE - statin-free individuals"]

MR_merged[,outcome2 := gsub("BC Survival","MR-IVW",outcome)]
MR_merged[method=="MR-ratio",outcome2 := gsub("MR-IVW","MR-ratio",outcome2)]
MR_merged[method=="MR-ratio" & !grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn<1,456\nrs562556")]
MR_merged[method=="MR-ratio" & grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn=91,868\nrs562556")]
MR_merged[method!="MR-ratio",outcome2 := paste0(outcome2,"\nn=91,868\nmultiple SNPs")]

#' Add columns necessary for plotting
MR_merged[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_merged[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_merged[grepl("PE",exposure),colGroup := "pQTL"]
MR_merged[grepl("GE",exposure),colGroup := "eQTL"]

MR_merged$signif <- ifelse(MR_merged$pval_IVW < 0.0045,1,0)
MR_merged[signif==1, y_ast:=upperCI95 + 2]
MR_merged[signif==1, y_ast2:=max(upperCI95, na.rm=TRUE) + 2]

#' # Forest Plot ####
#' ***

p1 = ggplot(data=MR_merged, aes(x=exposure, y=beta_IVW, ymin=lowerCI95, ymax=upperCI95,color=colGroup)) +
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
  geom_point(data = MR_merged[MR_merged$signif ==1, ],aes(x=exposure, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logHR for breast cancer survival (95% CI) \nper 1-SD increment in PCSK9 levels") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  scale_colour_manual(values = c("darkred","steelblue"), labels = c("GE", "PE")) +
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
        axis.text=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))

filename = paste0("../results/Main_Sup_Figs/Main_ForestPlot_MR_BCS.png")
png(filename = filename,width = 2500, height = 1500, res=200)
plot(p1)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
