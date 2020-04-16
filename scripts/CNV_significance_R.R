#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(rtracklayer))

## Runs on single samples

# cnvs_file = read.table(snakemake@input[["cnp"]], header = T) %>% as.tbl()
# ratio_file = snakemake@params[["withcontrol"]]

basepath = "/Users/diazmiyar/mnt/20200103-WGS/"
setwd (basepath)
#for (sample in c('WT_Ctrl','TP53_Ctrl','TP53_48Rec1','TP53_48Rec2','TP53_Aneu1','TP53_Aneu2')){
for (sample in c('K562_Control','K562_Treated1','K562_Treated2')){
  for (control in c("freec_single","freec_control")){
    #if (!(sample == "WT_Ctrl" & control == "freec_control")){ # skip that sample
    if (!(sample == "K562_Control" & control == "freec_control")){ # skip that sample
      message ("assessing significance in ", sample, " ",control)
      filepath = paste0 ("derived_data/", sample, "/", control, "/")
      
      dataTable <-read.table(paste0(filepath, "freec_ratio.txt"), header=TRUE);
      ratio<-data.frame(dataTable)
      
      dataTable <- read.table(paste0(filepath, "freec_CNVs"), header=FALSE)
      cnvs<- data.frame(dataTable)
      
      ratio$Ratio[which(ratio$Ratio==-1)]=NA
      
      cnvs.bed=GRanges(cnvs[,1],IRanges(cnvs[,2],cnvs[,3]))  
      ratio.bed=GRanges(ratio$Chromosome,IRanges(ratio$Start,ratio$Start),score=ratio$Ratio)
      
      overlaps <- subsetByOverlaps(ratio.bed,cnvs.bed)
      normals <- setdiff(ratio.bed,cnvs.bed)
      normals <- subsetByOverlaps(ratio.bed,normals)
      
      #mu <- mean(score(normals),na.rm=TRUE)
      #sigma<- sd(score(normals),na.rm=TRUE)
      
      #hist(score(normals),n=500,xlim=c(0,2))
      #hist(log(score(normals)),n=500,xlim=c(-1,1))
      
      #shapiro.test(score(normals)[which(!is.na(score(normals)))][5001:10000])
      #qqnorm (score(normals)[which(!is.na(score(normals)))],ylim=(c(0,10)))
      #qqline(score(normals)[which(!is.na(score(normals)))], col = 2)
      
      #shapiro.test(log(score(normals))[which(!is.na(score(normals)))][5001:10000])
      #qqnorm (log(score(normals))[which(!is.na(score(normals)))],ylim=(c(-6,10)))
      #qqline(log(score(normals))[which(!is.na(score(normals)))], col = 2)
      
      numberOfCol=length(cnvs)
      
      for (i in c(1:length(cnvs[,1]))) {
        values <- score(subsetByOverlaps(ratio.bed,cnvs.bed[i]))
        #wilcox.test(values,mu=mu)
        W <- function(values,normals){resultw <- try(wilcox.test(values,score(normals)), silent = TRUE)
        if(class(resultw)=="try-error") return(list("statistic"=NA,"parameter"=NA,"p.value"=NA,"null.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultw}
        KS <- function(values,normals){resultks <- try(ks.test(values,score(normals)), silent = TRUE)
        # In ks.test(values, score(normals)) :
        #   p-value will be approximate in the presence of ties
        if(class(resultks)=="try-error") return(list("statistic"=NA,"p.value"=NA,"alternative"=NA,"method"=NA,"data.name"=NA)) else resultks}
        #resultks <- try(KS <- ks.test(values,score(normals)), silent = TRUE)
        #	if(class(resultks)=="try-error") NA) else resultks
        cnvs[i,numberOfCol+1]=W(values,normals)$p.value
        cnvs[i,numberOfCol+2]=KS(values,normals)$p.value
      }
      
      if (numberOfCol==5) {
        names(cnvs)=c("chr","start","end","copy number","status","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
      }
      if (numberOfCol==7) {
        names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
      }
      if (numberOfCol==9) {
        names(cnvs)=c("chr","start","end","copy number","status","genotype","uncertainty","somatic/germline","precentageOfGermline","WilcoxonRankSumTestPvalue","KolmogorovSmirnovPvalue")  
      }
      
      write.table(cnvs, file=paste0(filepath, "freec_CNVs.p.value.txt"),sep="\t",quote=F,row.names=F)
      #write.table(cnvs, file=snakemake@output,sep="\t",quote=F,row.names=F)
    }
  }
}


