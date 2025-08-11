#' ---
#' title: "Supplemental Figure: Forest Plots of MR-IVW"
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

load("../results/05_MVMR.RData")
MVMRTab[,pval_adj1 := p.adjust(p=pval_exp1,method = "fdr")]
MVMRTab[,pval_adj2 := p.adjust(p=pval_exp2,method = "fdr")]

#' # Plot 1: MR-ratio for rs562556 ####
#' ***
plotData1 = copy(MRTab)
plotData1[,unique(outcome)]
plotData1 = plotData1[!grepl("BCS",outcome),]
plotData1 = plotData1[!grepl("gene expression",exposure),]
plotData1 = plotData1[method == "MR-ratio",]
plotData1[,sex := "both sexes"]
plotData1[grepl("females",exposure),sex := "females"]

#' update outcome and exposure names
plotData1[outcome == "BC - FinnGen + UKB", outcome := "BC [FinnGen + UKB]"]
plotData1[outcome == "CAD - all", outcome := "CAD [Aragam et al.]"]
plotData1[outcome == "CAD - females", outcome := "CAD F [Aragam et al.]"]
plotData1[outcome == "PLD - UKB", outcome := "PAD [UKB]"]

plotData1[grepl("PCSK9",exposure), exposure := "PCSK9"]
plotData1[grepl("LDL",exposure), exposure := "LDL-C"]

#' add columns necessary for plotting
plotData1[,lowerCI95 := beta-1.96*se]
plotData1[,upperCI95 := beta+1.96*se]

plotData1[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData1[signif==1, y_ast := 2]

plotData1[,dummy := paste(outcome,exposure,method)]
plotData1[,dummy2 := "PCSK9 function on other outcomes"]
plotData1 = plotData1[exposure=="PCSK9"]

p1 = ggplot(data=plotData1, aes(x=dummy, y=beta,
                                ymin=lowerCI95, ymax=upperCI95,
                                color=outcome,shape=method,
                                alpha=sex,fill=outcome)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.4),
                  size=1,linewidth=1)+ 
  geom_point(aes(x=dummy, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, 
             position = position_dodge2(width = 0.4),
             color="black",alpha=1)+
  facet_wrap(~dummy2,scales = "free_x")+ 
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in PCSK9 levels") + 
  #scale_y_continuous(limits = c(-20, 31)) + 
  #scale_colour_manual(values = c("darkorange","dodgerblue3","grey40")) +
  #scale_fill_manual(values = c("darkorange","dodgerblue3","grey40")) +
  scale_shape_manual(values = c(24)) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"), 
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
  labs(color = "Outcome source",shape = "Outcome genetic\nmodel", alpha = "Exposure sex")+ guides(fill = "none")

p1
filename = paste0("../results/SupplementalFigures/Plot1_MRratio_PCSK9.png")
png(filename = filename,width = 2000, height = 900, res=200)
plot(p1)
dev.off()

#' # Plot 2: MR-ratio & MR-IVW of PCSK9 ####
#' ***
plotData2 = copy(MRTab)
plotData2[,unique(outcome)]
plotData2 = plotData2[!grepl("BCS",outcome),]
plotData2 = plotData2[!grepl("ratio",exposure),]
plotData2 = plotData2[!grepl("LDL",exposure),]

#' update outcome and exposure names
plotData2[outcome == "BC - FinnGen + UKB", outcome := "BC [FinnGen + UKB]"]
plotData2[outcome == "CAD - all", outcome := "CAD [Aragam et al.]"]
plotData2[outcome == "CAD - females", outcome := "CAD F [Aragam et al.]"]
plotData2[outcome == "PLD - UKB", outcome := "PAD [UKB]"]

plotData2[, exposure := gsub("PCSK9 gene expression - ","",exposure)]
plotData2[, exposure := gsub("PCSK9 protein levels - all","Whole Blood \nprotein levels",exposure)]
plotData2[, exposure := gsub("PCSK9 protein levels - females","Whole Blood \nprotein levels (females)",exposure)]
plotData2[exposure == "Whole Blood", exposure := "Whole Blood \ngene expression"]

#' add columns necessary for plotting
plotData2[,lowerCI95 := beta-1.96*se]
plotData2[,upperCI95 := beta+1.96*se]
plotData2[,dummy := "PCSK9 level on other outcome"]

plotData2[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData2[signif==1, y_ast := upperCI95 + 0.1]

p2 = ggplot(data=plotData2, aes(x=exposure, y=beta, 
                                ymin=lowerCI95, ymax=upperCI95,
                                color=outcome,shape=method)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.6),
                  size=1,linewidth=1)+ 
  facet_wrap(~dummy,scales = "free_x")+ 
  # geom_point(aes(x=exposure, y=y_ast),
  #            shape = "*", size=8, show.legend = FALSE, 
  #            position = position_dodge2(width = 0.4),
  #            color="black",alpha=1)+
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in PCSK9 levels") + 
  #scale_colour_manual(values = c("darkorange","dodgerblue3")) +
  scale_shape_manual(values = c(16,17)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"), 
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
  labs(color = "Outcome source",shape = "MR method")

p2
filename = paste0("../results/SupplementalFigures/Plot2_MRratio_IVW.png")
png(filename = filename,width = 1800, height = 1500, res=200)
plot(p2)
dev.off()

#' # Plot 3: MR-IVW of LDL-C ####
#' ***
plotData3 = copy(MRTab)
plotData3[,unique(outcome)]
plotData3 = plotData3[!grepl("BCS",outcome),]
plotData3 = plotData3[!grepl("ratio",exposure),]
plotData3 = plotData3[grepl("LDL",exposure),]
plotData3[,sex := "both sexes"]
plotData3[grepl("females",exposure),sex := "females"]

#' update outcome and exposure names
plotData3[outcome == "BC - FinnGen + UKB", outcome := "BC [FinnGen + UKB]"]
plotData3[outcome == "CAD - all", outcome := "CAD [Aragam et al.]"]
plotData3[outcome == "CAD - females", outcome := "CAD F [Aragam et al.]"]
plotData3[outcome == "PLD - UKB", outcome := "PAD [UKB]"]

plotData3[, exposure := gsub("LDL-C levels - all ","",exposure)]
plotData3[, exposure := gsub("LDL-C levels - females ","",exposure)]
plotData3[, exposure := gsub("PCSK9_HMGCR","PCSK9 and HMGCR",exposure)]

#' add columns necessary for plotting
plotData3[,lowerCI95 := beta-1.96*se]
plotData3[,upperCI95 := beta+1.96*se]
plotData3[,dummy := "LDL-C level on other outcome"]
plotData3[,exposure2 := c(rep(c(1,2,4),each=4),rep(c(1,2,4),each=4),rep(3,8))]
plotData3[,exposure3 := paste(exposure2,exposure)]
plotData3[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData3[signif==1, y_ast:=upperCI95 + 0.1]

p3 = ggplot(data=plotData3, aes(x=exposure3, y=beta, 
                                ymin=lowerCI95, ymax=upperCI95,
                                color=outcome,alpha=sex)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.8),
                  size=1,linewidth=1)+ 
  facet_wrap(~dummy,scales = "free_x")+ 
  coord_flip() + 
  # geom_point(aes(x=exposure3, y=y_ast),
  #            shape = "*", size=8, show.legend = FALSE,
  #            position = position_dodge2(width = 0.6),
  #            color="black",alpha=1)+
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in LDL-C levels") + 
  #  scale_colour_manual(values = c("darkorange","dodgerblue3")) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"), 
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
  labs(color = "Outcome source",alpha = "Exposure sex")

p3
filename = paste0("../results/SupplementalFigures/Plot3_MRIVW_LDLC.png")
png(filename = filename,width = 2000, height = 1500, res=200)
plot(p3)
dev.off()

#' # PLot 4: MR and MVMR ####
#' ***
plotData4a = copy(MVMRTab)
plotData4a = plotData4a[!(outcome %in% c("BCS - BCAC","BCS - FinnGen")),]
plotData4a = plotData4a[condFstat_exp1 >9 & condFstat_exp2 >9,]
plotData4a[,sex := gsub("all","both sexes",sex)]
plotData4a[outcome == "BC - FinnGen + UKB", outcome := "BC [FinnGen + UKB]"]
plotData4a[outcome == "CAD - all", outcome := "CAD [Aragam et al.]"]
plotData4a[outcome == "CAD - females", outcome := "CAD F [Aragam et al.]"]
plotData4a[outcome == "PLD - UKB", outcome := "PAD [UKB]"]

plotData4a[SNPset=="PCSK9_v1",SNPset := "MVMR - PCSK9"]
plotData4a[SNPset=="PCSK9_HMGCR",SNPset := "MVMR - PCSK9 and HMGCR"]
plotData4a[,method := "MVMR-IVW"]

plotData4a = plotData4a[,c(3,1,4,19,6,7,8,9,15,16,17,2)]
names(plotData4a)[1:11] = names(plotData2)[1:11]
plotData4b = copy(plotData2)
plotData4b = plotData4b[grepl("protein",exposure),1:11]
plotData4b[,sex := "both sexes"]
plotData4b[grepl("fem",exposure),sex := "females"]
plotData4b[,exposure := "MR - PCSK9"]

plotData4 = rbind(plotData4b,plotData4a)

#' Add columns necessary for plotting
plotData4[,lowerCI95 := beta-1.96*se]
plotData4[,upperCI95 := beta+1.96*se]
plotData4[,dummy := "PCSK9 levels on other outcome, cond. on LDL-C"]

p4 = ggplot(data=plotData4, aes(x=exposure, y=beta, 
                                ymin=lowerCI95, ymax=upperCI95,
                                color=outcome,alpha=sex,shape=method)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.8),
                  size=1,linewidth=1)+ 
  facet_wrap(~dummy,scales = "free_x")+ 
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in PCSK9 levels") + 
  #scale_colour_manual(values = c("darkorange","dodgerblue3")) +
  scale_alpha_manual(values = c(0.5,1)) +
  scale_shape_manual(values = c(16,15)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"), 
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
  labs(color = "Outcome source",alpha = "Exposure sex",shape = "MR method")

p4
filename = paste0("../results/SupplementalFigures/Plot4_MVMRIVW_PCSK9.png")
png(filename = filename,width = 2000, height = 1500, res=200)
plot(p4)
dev.off()

#' # Plot 5: Combine 1 & 2 ###
#' ***
plotData12 = rbind(plotData1,plotData2,fill=T)

plotData12[,table(method)]

plotData12[is.na(sex),sex := "both sexes"]
plotData12[grepl("fem",exposure),sex := "females"]
plotData12[,table(sex)]

plotData12[,signif := ifelse(pval_adj < 0.05,1,0)]
plotData12[signif==1, y_ast:=upperCI95 + 0.1]
plotData12[,table(signif)]

plotData12[is.na(dummy2),dummy2 := dummy]
plotData12[exposure != "PCSK9",dummy := exposure]

plotData12[,table(dummy2)]
plotData12[grepl("function",dummy2),dummy2 := "PCSK9 function (rs562556)"]
plotData12[grepl("level",dummy2),dummy2 := "PCSK9 levels (QTLs)"]

plotData12[1:4,dummy]
plotData12[1:4,dummy := "WB PE"]
plotData12[5:8,dummy := "WB PE F"]
plotData12[9:12,dummy := "WB PE"]
plotData12[13:16,dummy := "WB PE F"]
plotData12[45:48,dummy := "WB GE"]

plotData12[,dummy := gsub("Brain - Cerebellar Hemisphere","Brain - CH",dummy)]
plotData12[,dummy := gsub("Brain - Cerebellum","Brain - Cb",dummy)]
plotData12[17:20,dummy := "Adipose VO"]

p12 = ggplot(data=plotData12, aes(x=dummy, y=beta,
                                  ymin=lowerCI95, ymax=upperCI95,
                                  color=outcome,shape=method,
                                  alpha=sex,fill=outcome)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.8),
                  size=1,linewidth=1)+ 
  geom_point(aes(x=dummy, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, 
             position = position_dodge2(width = 0.8),
             color="black",alpha=1)+
  facet_wrap(~dummy2,scales = "free")+ 
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in PCSK9 levels") + 
  #scale_y_continuous(limits = c(-20, 31)) + 
  # scale_colour_manual(values = c("darkorange","dodgerblue3","grey40")) +
  # scale_fill_manual(values = c("darkorange","dodgerblue3","grey40")) +
  scale_shape_manual(values = c(16,17)) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"),
        strip.text = element_text(size = 15),
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
        axis.text.y=element_text(family="Times New Roman",size=15, color = "Black"), 
        axis.text.x=element_text(family="Times New Roman",size=15, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))+ 
  labs(color = "Outcome \nsource",shape = "MR method", alpha = "Exposure sex")+ guides(fill = "none")

p12
filename = paste0("../results/SupplementalFigures/Plot12_PCSK9.png")
png(filename = filename,width = 2400, height = 2000, res=200)
plot(p12)
dev.off()


p12 = ggplot(data=plotData12[c(1:16,21:24,29:32,45:48)], aes(x=dummy, y=beta,
                                  ymin=lowerCI95, ymax=upperCI95,
                                  color=outcome,shape=method,
                                  alpha=sex,fill=outcome)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.6),
                  size=1,linewidth=1)+ 
  geom_point(aes(x=dummy, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, 
             position = position_dodge2(width = 0.6),
             color="black",alpha=1)+
  facet_wrap(~dummy2,scales = "free")+ 
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in PCSK9 levels") + 
  #scale_y_continuous(limits = c(-20, 31)) + 
  # scale_colour_manual(values = c("darkorange","dodgerblue3","grey40")) +
  # scale_fill_manual(values = c("darkorange","dodgerblue3","grey40")) +
  scale_shape_manual(values = c(16,17)) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"),
        strip.text = element_text(size = 15),
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
        axis.text.y=element_text(family="Times New Roman",size=15, color = "Black"), 
        axis.text.x=element_text(family="Times New Roman",size=15, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))+ 
  labs(color = "Outcome \nsource",shape = "MR method", alpha = "Exposure sex")+ guides(fill = "none")

p12
filename = paste0("../results/SupplementalFigures/Plot12_PCSK9_filtered.png")
png(filename = filename,width = 2400, height = 1500, res=200)
plot(p12)
dev.off()

#' # Plot 6: 3 & 4 combined ####
#' ***
plotData34 = rbind(plotData3,plotData4,fill=T)

plotData34[,table(method)]
plotData34[,table(sex)]

plotData34[,signif := ifelse(pval_adj < 0.05,1,0)]

plotData34[is.na(exposure3),exposure3 := exposure]
plotData34[,table(exposure3)]
plotData34[exposure3 == "MVMR - PCSK9 and HMGCR",exposure3 := "3 PCSK9 and HMGCR SNPs"]
plotData34[exposure3 == "MVMR - PCSK9",exposure3 := "2 PCSK9 SNPs"]
plotData34 = plotData34[exposure3 != "MR - PCSK9"]

plotData34[,table(dummy)]
plotData34[,dummy := gsub("cond","\ncond",dummy)]
plotData34[grepl("LDL-C level",dummy),dummy := paste0(dummy,"\n using multiple QTLs")]

p34 = ggplot(data=plotData34, aes(x=exposure3, y=beta,
                                  ymin=lowerCI95, ymax=upperCI95,
                                  color=outcome,shape=method,
                                  alpha=sex,fill=outcome)) +
  geom_hline(yintercept=0, lty=2, linewidth =1) +  
  geom_pointrange(position = position_dodge2(width = 0.8),
                  size=1,linewidth=1)+ 
  geom_point(aes(x=exposure3, y=y_ast),
             shape = "*", size=8, show.legend = FALSE, 
             position = position_dodge2(width = 0.8),
             color="black",alpha=1)+
  facet_wrap(~dummy,scales = "free_x")+ 
  coord_flip() + 
  ggtitle("")+
  xlab("") +
  ylab("logOR for BC or CAD Risk or 1-SD increment for PAD (95% CI)\nper 1-SD increment in exposure levels") + 
  #scale_y_continuous(limits = c(-20, 31)) + 
  # scale_colour_manual(values = c("darkorange","dodgerblue3","grey40")) +
  # scale_fill_manual(values = c("darkorange","dodgerblue3","grey40")) +
  scale_shape_manual(values = c(16,15)) +
  scale_alpha_manual(values = c(0.5,1)) +
  theme(line = element_line(colour = "black", linewidth = 1), 
        strip.background = element_rect(fill="grey80"), 
        strip.text = element_text(size = 15),
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
        axis.text.y=element_text(family="Times New Roman",size=15, color = "Black"), 
        axis.text.x=element_text(family="Times New Roman",size=15, color = "Black"), 
        text=element_text(family="Times New Roman",size=15), plot.margin = margin(t = 0.25, r = 1, b = 0.5, l = 0.25, unit = "cm"))+ 
  labs(color = "Outcome source",shape = "MR method", alpha = "Exposure sex")+ guides(fill = "none")

p34
filename = paste0("../results/SupplementalFigures/Plot34_LDLC_PCSK9.png")
png(filename = filename,width = 2800, height = 1500, res=200)
plot(p34)
dev.off()

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")
