#' ---
#' title: "Create Poster Figure: Forest Plots of MR-IVW"
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
#' What to show on the poster? just one plot, just females, facet per method, colour per outcome, shape per exposure, fill per drug target
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
MRTab[grepl("FinnGen",outcome), outcome := "BC Survival [FinnGen]"]
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
MRTab[method=="MR-ratio" & grepl("Mei",outcome),outcome2 := paste0(outcome2,"\nn<1,456\nrs562556")]
MRTab[method=="MR-ratio" & grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn=91,868\nrs562556")]
MRTab[method=="MR-ratio" & grepl("FinnGen",outcome),outcome2 := paste0(outcome2,"\nn=4648\nrs562556")]
MRTab[method!="MR-ratio" & grepl("BCAC",outcome),outcome2 := paste0(outcome2,"\nn=91,868\nmultiple SNPs")]
MRTab[method!="MR-ratio" & grepl("FinnGen",outcome),outcome2 := paste0(outcome2,"\nn=4648\nmultiple SNPs")]

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

#' # Forest Plots ####
#' ***
#' 
#' ## Test Plot ####
plotData = copy(MRTab)
plotData = plotData[grepl("females",exposure),]
plotData = plotData[!grepl("any",exposure),]
plotData[,exposure := c(rep("PCSK9",5),rep("LDL-C",7))]
plotData[,locus := c(rep("PCSK9",10),rep("HMGCR",2))]
plotData[,outcome := gsub("BC Survival ","",outcome)]
plotData[,dummy := paste(outcome,exposure,locus)]
setorder(plotData,dummy)
plotData[1:9,dummy2 := paste0("res 0",1:9)]
plotData[10:12,dummy2 := paste0("res 1",0:2)]
plotData[dummy2 == "res 10",dummy2 := "res 09"]
plotData[dummy2 == "res 08",dummy2 := "res 07"]
plotData[dummy2 == "res 05",dummy2 := "res 04"]
plotData[dummy2 == "res 03",dummy2 := "res 02"]
plotData[y_ast<0.5,y_ast := 0.5]

plotData[method=="MR-ratio",method := paste(method,"SNP rs562556",sep=" - ")]
plotData[method=="MR-IVW",  method := paste(method,"multiple SNPs",sep=" - ")]
# plotData[,beta := beta*(-1)]
# plotData[,lowerCI95 := beta-1.96*se]
# plotData[,upperCI95 := beta+1.96*se]

p = ggplot(data=plotData, aes(x=dummy2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=outcome,shape=exposure,fill=locus)) +
  geom_pointrange()+ 
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  facet_wrap(~method,scales = "free_x")+ 
  coord_flip() + 
  geom_point(size = 3) + 
  geom_point(data = plotData[plotData$signif ==1, ],aes(x=dummy2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
  ggtitle("")+
  xlab("") +
  ylab("logHR for breast cancer survival (95% CI) \nper 1-SD increment in exposure levels") + 
  scale_colour_manual(values = c("darkorange","dodgerblue3","darkgrey")) +
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values = c("lightgrey","black"),)+
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="gray90"), 
        #legend.position ="none", 
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
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))+ 
  labs(fill = "Drug target",color = "Outcome source",shape = "Exposure")

p
filename = paste0("../results/PosterFigures/MR_BCS_females.png")
png(filename = filename,width = 1800, height = 900, res=200)
plot(p)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
