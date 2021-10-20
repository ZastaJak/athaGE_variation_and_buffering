#!/usr/bin/Rscript

library(DMRcaller)

bare_grange_csv<-read.csv("bare_atha_grange.csv")

bare_atha_grange<-GRanges(seqnames = bare_grange_csv[,2], ranges = IRanges(start = as.numeric(bare_grange_csv[,3]), end = as.numeric(bare_grange_csv[,4]), names = bare_grange_csv[,1]), strand = bare_grange_csv[,6])


setwd("WT/SRR4999961/")

WT<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")

chloroplast_data<-WT[seqnames(WT)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
WT_corrected<-WT
WT_corrected$readsM<-round((WT$readsM)-(WT$readsN)*(1-rate))
WT_corrected$readsM[WT_corrected$readsM<0]<-0
WT_corrected$readsN<-round(WT$readsN*rate)

rm(WT)

setwd("../../met1-1/SRR4999962")

met1_1<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")

chloroplast_data<-met1_1[seqnames(met1_1)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_1_corrected<-met1_1
met1_1_corrected$readsM<-round((met1_1$readsM)-(met1_1$readsN)*(1-rate))
met1_1_corrected$readsM[met1_1_corrected$readsM<0]<-0
met1_1_corrected$readsN<-round(met1_1$readsN*rate)

rm(met1_1)

setwd("../../met1-3/SRR4999963")

met1_3<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")

chloroplast_data<-met1_3[seqnames(met1_3)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_3_corrected<-met1_3
met1_3_corrected$readsM<-round((met1_3$readsM)-(met1_3$readsN)*(1-rate))
met1_3_corrected$readsM[met1_3_corrected$readsM<0]<-0
met1_3_corrected$readsN<-round(met1_3$readsN*rate)

rm(met1_3)

setwd("../../")

seqlevels(WT_corrected)<-c("1", "2", "3", "4", "5", "Pt", "Mt")
seqlevels(met1_1_corrected)<-c("1", "2", "3", "4", "5", "Pt", "Mt")
seqlevels(met1_3_corrected)<-c("1", "2", "3", "4", "5", "Pt", "Mt")

DMRsCG_WT_to_met1_1 <- computeDMRs(WT_corrected,
                          met1_1_corrected,
                          regions = NULL,
                          context = "CG",
                          method = "bins",
                          binSize = 100,
                          test = "score",
                          pValueThreshold = 0.01,
                          minCytosinesCount = 4,
                          minProportionDifference = 0.4,
                          minGap = 200,
                          minSize = 50,
                          minReadsPerCytosine = 4,
                          cores = 1)

DMRsCG_WT_to_met1_3 <- computeDMRs(WT_corrected,
                          met1_3_corrected,
                          regions = NULL,
                          context = "CG",
                          method = "bins",
                          binSize = 100,
                          test = "score",
                          pValueThreshold = 0.01,
                          minCytosinesCount = 4,
                          minProportionDifference = 0.4,
                          minGap = 200,
                          minSize = 50,
                          minReadsPerCytosine = 4,
                          cores = 1)

DMRsCG_met1_1_to_met1_3 <- computeDMRs(met1_1_corrected,
                          met1_3_corrected,
                          regions = NULL,
                          context = "CG",
                          method = "bins",
                          binSize = 100,
                          test = "score",
                          pValueThreshold = 0.01,
                          minCytosinesCount = 4,
                          minProportionDifference = 0.4,
                          minGap = 200,
                          minSize = 50,
                          minReadsPerCytosine = 4,
                          cores = 1)

save(DMRsCG_WT_to_met1_1, file = "CG-WT-met1_1.RData")
save(DMRsCG_WT_to_met1_3, file = "CG-WT-met1_3.RData")
save(DMRsCG_met1_1_to_met1_3, file = "CG-met1_1-met1_3.RData")