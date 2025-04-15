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

#' # Get MR results
#' ***
load("../results/SupTables.RData")

#' What do I want? 
#' 
#' - BCS and CAD 
#' - PCSK9 (females) and LDLC (females)
#' - MR and MVMR 
#' 
mrtab1 = copy(stab2)
mrtab1 = mrtab1[exposure == "PCSK9_free"]
mrtab1 = mrtab1[grepl("Survival",outcome) | grepl("CAD",outcome) ]

mrtab2 = copy(stab5)
mrtab2 = mrtab2[flag %in% c("flag4")]
mrtab2 = mrtab2[grepl("Survival",outcome) | grepl("CAD",outcome) ]

mrtab3 = copy(stab7)
mrtab3 = mrtab3[exposure1 == "PCSK9_free"]
mrtab3 = mrtab3[grepl("Survival",outcome) | grepl("CAD",outcome) ]
mrtab3 = mrtab3[flag == "flag1"]

MR_merged = rbind(mrtab1,mrtab2,fill=T)
MR_merged[,flag := "MR"]

mrtab4 = copy(mrtab3)
mrtab4 = mrtab4[,c(3,1,2,4,5,6,7,13,14,15)]
mrtab5 = copy(mrtab3)
mrtab5 = mrtab5[,c(8,1,2,9,10,11,12,13,14,15)]
names(mrtab4) = names(MR_merged)
names(mrtab5) = names(MR_merged)
mrtab4[,flag := "MVMR"]
mrtab5[,flag := "MVMR"]

MR_merged = rbind(MR_merged,mrtab4,mrtab5)
MR_merged[, outcome := gsub("BCSurvival","Breast Cancer Survival",outcome)]
MR_merged[, outcome := gsub("CAD","Coronary Atherosclerosis",outcome)]
MR_merged[, exposure := gsub("_females","",exposure)]
MR_merged[, exposure := gsub("_free","",exposure)]

#' Add columns necessary for plotting
MR_merged[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_merged[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_merged[!grepl("LDLC",exposure),colGroup := "pQTL"]
MR_merged[grepl("LDLC",exposure),colGroup := "eQTL"]

MR_merged$signif <- ifelse(MR_merged$pval_IVW < 0.05,1,0)
MR_merged[signif==1, y_ast:=upperCI95 + 0.1]
MR_merged[signif==1, y_ast2:=max(upperCI95, na.rm=TRUE) + 0.1,by=outcome]

MR_merged[,exposure2 := paste(flag,exposure,sep = " - ")]

#' # Forest Plot ####
#' ***

p1 = ggplot(data=MR_merged, aes(x=exposure2, y=beta_IVW, ymin=lowerCI95, ymax=upperCI95,color=colGroup)) +
  # Make range for ggplot values based on the data and AES specified in first line:
  geom_pointrange()+ 
  # add a dotted line at x=0 after flip:
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  # Makes whiskers on the range (more aesthetically pleasing):
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  # Makes DV header (Can handle multiple DVs):
  facet_wrap(~outcome,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MR_merged[MR_merged$signif ==1, ],aes(x=exposure2, y=y_ast2),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logHR for breast cancer survival (95% CI)                    logOR for coronary atherosclerosis risk (95% CI) \nper 1-SD increment in exposure levels") + 
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
        text=element_text(family="Times New Roman",size=10), plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "cm"))

filename = paste0("../results/Main_Sup_Figs/Main_ForestPlot_MVMR_BCS_CAD.png")
png(filename = filename,width = 2000, height = 1000, res=200)
plot(p1)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
