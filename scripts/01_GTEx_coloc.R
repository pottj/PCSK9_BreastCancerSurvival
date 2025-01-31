#' ---
#' title: "Check GTEx v8 data"
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
#' In this script, I load the GTEx v8 data as generated in the [sex- and statin-stratified meta-GWAS project](https://github.com/GenStatLeipzig/GWAMA_PCSK9_strat) (script *05_1_coloc_get_eQTL_data.R*) 
#' 
#' I will first load each data set, filter for PCSK9 gene, get the number of genome-wide significant SNPs (potential IVs), and the perform one hypercoloc approach, to check how similar the signals are. Finally, I save the genome-wide significant eQTLs for later use. 
#' 
#' # Initialize ####
#' ***
rm(list = ls())
time0<-Sys.time()
server = "laptop_BSU"

source("../SourceFile.R")
.libPaths()

#' # Load data ####
#' ***
myFiles = list.files(path = GTEx)
myFiles

dumTab1 = foreach(i = 1:length(myFiles))%do%{
  #i=1
  load(paste0(GTEx,myFiles[i]))
  data0 = data0[gene == "PCSK9",]
  data0
}
GTExData = rbindlist(dumTab1)
GTExData[,table(tissue,pval<5e-8)]
GTExData[,table(tissue,pval<1e-6)]
GTExData[,min(pval),by=tissue]

#' Although expressed in these tissues, there is not always a strong instrument. Colon (transverse), esophagus (mucosa), and small intestine (terminal ileum) have no eQTLs with at least suggestive significance. Liver and pancreas lack of genome-wide significant SNPs (p<5e-8), but have at least suggestive significant eQTLs (p<1e-6).
#' 
#' # Colocalization ####
#' ***
#' I will test all combination of the 6 tissues with at least suggestive significant SNPs.
#' 
myTissues = GTExData[pval<1e-6,unique(tissue)]

dumTab2 = foreach(i = 1:(length(myTissues)-1))%do%{
  #i=1
  data1 = copy(GTExData)
  data1 = data1[tissue == myTissues[i],]
  
  dumTab3 = foreach(j = (i+1):length(myTissues))%do%{
    #j=2
    data2 = copy(GTExData)
    data2 = data2[tissue == myTissues[j],]
    
    # filter for overlapping SNPs
    data1 = data1[variant_id_b38 %in% data2$variant_id_b38,]
    data2 = data2[variant_id_b38 %in% data1$variant_id_b38,]
    stopifnot(data1$variant_id_b38 == data2$variant_id_b38)
    
    # run coloc
    my_res<- coloc::coloc.abf(dataset1=list(beta=data1$beta,
                                            varbeta=(data1$se)^2,
                                            N=data1$n_samples,
                                            snp=data1$variant_id_b38,
                                            MAF=data1$maf,
                                            type="quant"),
                              dataset2=list(beta=data2$beta,
                                            varbeta=(data2$se)^2,
                                            N=data2$n_samples,
                                            snp=data2$variant_id_b38,
                                            MAF=data2$maf,
                                            type="quant"))
    my_res2<-my_res[[1]]
    my_res3<-as.data.table(my_res2)
    my_res3<-t(my_res3)
    my_res3<-as.data.table(my_res3)
    names(my_res3)<-names(my_res2)
    
    my_res3[,trait1:=myTissues[i]]
    my_res3[,trait2:=myTissues[j]]
    my_res3
    
  }
  ColocResults1 = rbindlist(dumTab3)
  ColocResults1
}
ColocResults = rbindlist(dumTab2)
ColocResults[PP.H4.abf>0.8,]
ColocResults[PP.H3.abf>0.8,]
ColocResults[PP.H4.abf>0.5,]
ColocResults[PP.H3.abf>0.5,]

#' There are three cluster: 
#' 
#' 1. trio of brain (cerebellar hemisphere), brain (cerebellum), and pancreas, which share the same signal, 
#' 2. pair of lung and nerve (tibial), which share the same signal
#' 3. liver, which does not share its signal with the other five tissues. 
#'  
#' # Save data ####
#' ***
save(GTExData,file = "../temp/GTExV8_selectedTissues_PCSK9.RData")

GTExData = GTExData[pval<1e-6,]
GTExData[,table(tissue)]
GTExData[,length(unique(variant_id_b38))]
save(GTExData,file = "../temp/GTExV8_potentialIVs.RData")

#' # Session Info ####
#' ***
sessionInfo()
message("\nTOTAL TIME : " ,round(difftime(Sys.time(),time0,units = "mins"),3)," minutes")

