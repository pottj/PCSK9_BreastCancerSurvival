#' ---
#' title: "Create Sup Figures: Forest Plots of MR-IVW"
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
#' What to show in the supplemental figures? 
#' 
#' - Forest Plot of the MR results (both IVW and ratio)
#' - Per exposures (females & all together)
#' - For all outcomes (one each row) 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"
load_meta = F

source("../SourceFile.R")

#' # No gene expression ####
#' ***
load("../results/04_MR.RData")

#' Add columns necessary for plotting
MRTab[,pval_adj := p.adjust(p=pval,method = "fdr")]
MRTab[,lowerCI95 := beta-1.96*se]
MRTab[,upperCI95 := beta+1.96*se]

MRTab[, signif := ifelse(pval_adj < 0.05,1,0)]
MRTab[signif==1, y_ast:=upperCI95 + 1]
MRTab[signif==1 & method=="MR-IVW", y_ast:=upperCI95 + 0.1]

#' Remove Mei et al. (not relevant for the remaining outcomes)
MRTab = MRTab[!grepl("Mei ",outcome),]

#' Also remove Gene expression - they will get their own plot!
MRTab = MRTab[!grepl("expression",exposure),]

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

#' Get color group
MRTab[grepl("PCSK9",exposure),colGroup := "PCSK9"]
MRTab[grepl("LDLC",exposure),colGroup := "LDLC"]
MRTab[,sex := "females"]
MRTab[grepl("all",exposure),sex := "all"]
MRTab[,colGroup2 := paste0(colGroup,sex,sep="::")]

#' For CAD, I only want to compare all with all and females with females
filt1 = grepl("females",MRTab$sex) & MRTab$outcome == "CAD - all"
table(filt1)
filt2 = grepl("all",MRTab$sex) & MRTab$outcome == "CAD - females"
table(filt2)
filt = filt1 | filt2
MRTab = MRTab[!filt,]
MRTab[grepl("CAD",outcome),outcome := "CAD - Aragam et al."]

#' Loop 
myOutcomes = unique(MRTab$outcome)
myOutcomes2 = c("BCS","BC","PAaD","CAD")

for(i in 1:length(myOutcomes)){
  #i=4
  plotData = copy(MRTab)
  plotData = plotData[outcome == myOutcomes[i],]
  setorder(plotData,-exposure)
  plotData[,exposure := paste(c(8,8,7,7,2,6,6,4,1,5,5,3),exposure)]
  
  # get y-axis label
  if(i==1) myYlab = "logHR for BC Survival (95% CI) \nper 1-SD increment in exposure levels"
  if(i==2) myYlab = "logOR for BC risk (95% CI) \nper 1-SD increment in exposure levels"
  if(i==3) myYlab = "1-SD increment for Parental Age at Death (95% CI) \nper 1-SD increment in exposure levels"
  if(i==4) myYlab = "logOR for CAD risk (95% CI) \nper 1-SD increment in exposure levels"
  
  # plot
  p1 = ggplot(data=plotData, aes(x=exposure, y=beta, ymin=lowerCI95, ymax=upperCI95,color=colGroup2)) +
    # Make range for ggplot values based on the data and AES specified in first line:
    geom_pointrange()+ 
    # add a dotted line at x=0 after flip:
    geom_hline(yintercept=0, lty=2, linewidth =1) +  
    # Makes whiskers on the range (more aesthetically pleasing):
    geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
    # Makes DV header (Can handle multiple DVs):
    facet_wrap(~method,scales = "free_x")+ 
    # flip coordinates (puts labels on y axis):
    coord_flip() + 
    # specifies the size and shape of the geompoint: 
    geom_point(shape = 15, size = 2) + 
    # add asteristik for the significant results: 
    geom_point(data = plotData[signif ==1, ],aes(x=exposure, y=y_ast),
               shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
    # Blank Title for the Graph:
    ggtitle("")+
    # Label on the Y axis (flipped specification do to coord_flip):
    xlab("") +
    # Label on the X axis (flipped specification do to coord_flip):
    ylab(myYlab) + 
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
          axis.title=element_text(size=15),
          axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
          plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
          axis.text.y=element_text(family="Times New Roman",size=10, color = "Black",face = "italic"), 
          axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
          text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))
  
  p1
  filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MR_",myOutcomes2[i],".png")
  png(filename = filename,width = 2500, height = 1500, res=200)
  plot(p1)
  dev.off() 
}

#' # Only gene expression ####
#' ***
load("../results/04_MR.RData")

#' Add columns necessary for plotting
MRTab[,pval_adj := p.adjust(p=pval,method = "fdr")]
MRTab[,lowerCI95 := beta-1.96*se]
MRTab[,upperCI95 := beta+1.96*se]

MRTab[, signif := ifelse(pval_adj < 0.05,1,0)]
MRTab[signif==1 , y_ast:= 0.3]
MRTab[signif==1 & outcome == "PLD - UKB", y_ast:= 0.05]

#' Remove Mei et al. (not relevant for the remaining outcomes)
MRTab = MRTab[!grepl("Mei ",outcome),]

#' Restrict for gene expression
MRTab = MRTab[grepl("expression",exposure),]

#' add nice name for phenotype
unique(MRTab$exposure)
unique(MRTab$outcome)

MRTab[grepl("BCAC",outcome), outcome := "BC Survival [BCAC]"]
MRTab[grepl("FinnGen",outcome), outcome := "BC risk [FinnGen + UKB]"]
MRTab[grepl("PLD ",outcome), outcome := "Parental Age at Death [UKB]"]
MRTab[grepl("CAD - fe",outcome), outcome := "CAD risk in females [Aragam et al.]"]
MRTab[grepl("CAD - all",outcome), outcome := "CAD risk [Aragam et al.]"]

MRTab[,exposure := gsub("ratio","(rs562556)",exposure)]
MRTab[,exposure := gsub("PCSK9 gene expression - ","",exposure)]

setorder(MRTab,exposure)
MRTab[,exposure := paste(c(1,1,1,1,1,
                           2,2,2,2,2,
                           3,3,3,3,3,
                           9,9,9,9,9,
                           4,4,4,4,4,
                           5,5,5,5,5,
                           6,6,6,6,6,
                           7,7,7,7,7,
                           99,99,99,99,99,
                           8,8,8,8,8),exposure)]

# plot for BC risk and BC survival 
p1 = ggplot(data=MRTab[grepl("BC",outcome),], aes(x=exposure, y=beta, ymin=lowerCI95, ymax=upperCI95,col=signif)) +
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
  geom_point(data = MRTab[signif ==1 &grepl("BC",outcome), ],aes(x=exposure, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logHR for breast cancer survival (95% CI) \nlogOR for breast cancer risk (95% CI) \nper 1-SD increment in PCSK9 gene expression") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  # scale_colour_manual(values = c("darkorange","darkorange4","dodgerblue","dodgerblue3")) +
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
        axis.title=element_text(size=15),
        axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))

p1

# plot for CAD risk and parental age at death
p2 = ggplot(data=MRTab[!grepl("BC",outcome),], aes(x=exposure, y=beta, ymin=lowerCI95, ymax=upperCI95,col=signif)) +
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
  geom_point(data = MRTab[signif ==1 & !grepl("BC",outcome), ],aes(x=exposure, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
  # Blank Title for the Graph:
  ggtitle("")+
  # Label on the Y axis (flipped specification do to coord_flip):
  xlab("") +
  # Label on the X axis (flipped specification do to coord_flip):
  ylab("logOR for CAD risk (95% CI) \n1-SD increment for Parental Age at Death (95% CI) \nper 1-SD increment in PCSK9 gene expression") + 
  # limits and tic marks on X axis (flipped specification do to coord_flip):
  #scale_y_continuous(limits = c(-.50,.50), breaks = c(-.50,-.25,0,.25,.50))+
  # add color coding: 
  # scale_colour_manual(values = c("darkorange","darkorange4","dodgerblue","dodgerblue3")) +
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
        axis.title=element_text(size=15),
        axis.title.x = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 20, r = 0, b = 0, l = 0)),
        axis.title.y = element_text(family="Times New Roman",colour = "Black", margin = margin(t = 0, r = 20, b = 0, l = 0)),
        plot.title = element_text(family="Times New Roman", colour = "Black", margin = margin(t = 0, r = 0, b = 20, l = 0)),
        axis.text.y=element_text(family="Times New Roman",size=10, color = "Black"), 
        axis.text.x=element_text(family="Times New Roman",size=10, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))

p2

filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MR_GE_BC_BCS.png")
png(filename = filename,width = 2500, height = 1500, res=200)
plot(p1)
dev.off() 

filename = paste0("../results/Main_Sup_Figs/Sup_ForestPlot_MR_GE_CAD_PLD.png")
png(filename = filename,width = 2500, height = 1500, res=200)
plot(p2)
dev.off() 

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
