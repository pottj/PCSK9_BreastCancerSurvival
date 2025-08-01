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

MVMRTab = MVMRTab[outcome %in% c("BCS - BCAC","BCS - FinnGen")]
MVMRTab = MVMRTab[grepl("female",sex),]
MVMRTab = MVMRTab[condFstat_exp1 >9 & condFstat_exp2 >9,]

MVMRTab[SNPset=="PCSK9_v1",SNPset := "PCSK9"]
MVMRTab[SNPset=="PCSK9_HMGCR",SNPset := "PCSK9 + HMGCR"]

MVMRTab1 = copy(MVMRTab)
MVMRTab1 = MVMRTab1[,c(1:4,5:9,17)]
MVMRTab2 = copy(MVMRTab)
MVMRTab2 = MVMRTab2[,c(1:4,10:14,18)]
setnames(MVMRTab1,"exposure1","exposure")
setnames(MVMRTab1,"pval_adj1","pval_adj")
names(MVMRTab1) = gsub("_exp1","",names(MVMRTab1))
setnames(MVMRTab2,"exposure2","exposure")
setnames(MVMRTab2,"pval_adj2","pval_adj")
names(MVMRTab2) = gsub("_exp2","",names(MVMRTab2))

MVMR = rbind(MVMRTab1,MVMRTab2)
setorder(MVMR,SNPset,outcome,exposure)
MVMR[,outcome := gsub("BCS - ","[",outcome)]
MVMR[,outcome := paste0(outcome,"]")]

#' Add columns necessary for plotting
MVMR[,lowerCI95 := beta-1.96*SE]
MVMR[,upperCI95 := beta+1.96*SE]

MVMR$signif <- ifelse(MVMR$pval_adj < 0.05,1,0)
MVMR[signif==1, y_ast:=upperCI95 + 0.4]
MVMR[signif==1 & grepl("CAD",outcome), y_ast:=upperCI95 + 0.2]

MVMR[,dummy := paste(SNPset,outcome,exposure)]
setorder(MVMR,dummy)
MVMR[,dummy2 := paste0("res 0",1:8)]
MVMR[,method := "MVMR-IVW"]

p = ggplot(data=MVMR, aes(x=dummy2, y=beta, ymin=lowerCI95, ymax=upperCI95,color=outcome,shape=exposure,fill=SNPset)) +
  geom_pointrange()+ 
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_errorbar(aes(ymin=lowerCI95, ymax=upperCI95), width=0.5, cex=1)+ 
  geom_segment(aes(x = 8, y = 0, xend = 8, yend = 3.9), size=1,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 8, y = 0, xend = 8, yend = -3.9), size=1,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 7, y = 0, xend = 7, yend = 3.9), size=1,
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 7, y = 0, xend = 7, yend = -3.9), size=1,
               arrow = arrow(length = unit(0.25, "cm"))) +
  facet_wrap(~method,scales = "free_x")+ 
  coord_flip() + 
  geom_point(size = 3) + 
  geom_point(data = MVMR[MVMR$signif ==1, ],aes(x=dummy2, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, color="black",alpha=1) +
  ggtitle("")+
  xlab("") +
  scale_y_continuous(limits = c(-4, 4)) + 
  
  ylab("logHR for breast cancer survival (95% CI) \nper 1-SD increment in exposure levels") + 
  scale_colour_manual(values = c("darkorange","dodgerblue3")) +
  scale_shape_manual(values = c(24,21)) +
  scale_fill_manual(values = c("black","white"),)+
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
  labs(fill = "Drug targets",color = "Outcome source",shape = "Exposure")

p
filename = paste0("../results/PosterFigures/MVMR_BCS_females.png")
png(filename = filename,width = 1500, height = 900, res=200)
plot(p)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
