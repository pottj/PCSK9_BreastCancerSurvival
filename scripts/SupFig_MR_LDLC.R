#' ---
#' title: "Create Sup Figure: Forest Plots of MR-IVW"
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
.libPaths()

#' # Load data ####
#' ***
load("../results/SupTables.RData")

MR_IVW_PCSK9 = copy(stab2)
MR_IVW_LDLC = copy(stab5)

#' # Filter ####
#' ***
MR_IVW_PCSK9[,unique(exposure)]
MR_IVW_PCSK9 = MR_IVW_PCSK9[grepl("PCSK9",exposure),]
MR_IVW_PCSK9[,exposure := gsub("_"," - ",exposure)]
MR_IVW_PCSK9[,exposure := gsub("free","statin-free",exposure)]
MR_IVW_PCSK9[,exposure := paste(exposure,"PCSK9 SNPs",sep = " - ")]

MR_IVW_LDLC[,unique(exposure)]
MR_IVW_LDLC[,exposure := gsub("_"," - ",exposure)]
MR_IVW_LDLC[flag == "flag1",exposure := paste(exposure,"pruned",sep = " - ")]
MR_IVW_LDLC[flag == "flag2",exposure := paste(exposure,"PCSK9 SNPs",sep = " - ")]
MR_IVW_LDLC[flag == "flag3",exposure := paste(exposure,"HMGCR SNPs",sep = " - ")]
MR_IVW_LDLC[flag == "flag4",exposure := paste(exposure,"PCSK9 & HMGCR SNPs",sep = " - ")]

#' Merge data sets - what columns do I really need?
MR_merged = rbind(MR_IVW_LDLC,MR_IVW_PCSK9,fill=T)

#' Add columns necessary for plotting
MR_merged[,lowerCI95 := beta_IVW-1.96*se_IVW]
MR_merged[,upperCI95 := beta_IVW+1.96*se_IVW]

MR_merged[!grepl("LDLC",exposure),colGroup := "pQTL"]
MR_merged[grepl("LDLC",exposure),colGroup := "eQTL"]

MR_merged$signif <- ifelse(MR_merged$pval_IVW < 0.05,1,0)
MR_merged[signif==1, y_ast:=upperCI95 + 0.1,by=outcome]
MR_merged[signif==1, y_ast2:=max(upperCI95, na.rm=TRUE) + 0.1,by=outcome]
# MR_merged[signif==1, y_ast3:=min(lowerCI95, na.rm=TRUE) - 0.5]

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
  facet_wrap(~outcome,scales = "free_x")+ 
  # flip coordinates (puts labels on y axis):
  coord_flip() + 
  # specifies the size and shape of the geompoint: 
  geom_point(shape = 15, size = 2) + 
  # add asteristik for the significant results: 
  geom_point(data = MR_merged[MR_merged$signif ==1, ],aes(x=exposure, y=y_ast2),
             shape = "*", size=8, show.legend = FALSE, color="black") +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logHR, logOR, or 1-SD increment in outcome (95% CI) \nper 1-SD increment in exposure levels") + 
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

filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MR_LDLC.png")
png(filename = filename,width = 3500, height = 2000, res=200)
plot(p1)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
