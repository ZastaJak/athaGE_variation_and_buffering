#!/usr/bin/Rscript

library(DMRcaller)

setwd("../WT/SRR4999961/")
WT<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../../met1-1/SRR4999962/")
met1_1<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../../met1-3/SRR4999963/")
met1_3<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../../comp_old/")

chloroplast_data<-WT[seqnames(WT)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
WT_corrected<-WT
WT_corrected$readsM<-round((WT$readsM)-(WT$readsN)*(1-rate))
WT_corrected$readsM[WT_corrected$readsM<0]<-0
WT_corrected$readsN<-round(WT$readsN*rate)

WT<-WT_corrected
rm(WT_corrected)

chloroplast_data<-met1_1[seqnames(met1_1)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_1_corrected<-met1_1
met1_1_corrected$readsM<-round((met1_1$readsM)-(met1_1$readsN)*(1-rate))
met1_1_corrected$readsM[met1_1_corrected$readsM<0]<-0
met1_1_corrected$readsN<-round(met1_1$readsN*rate)

met11<-met1_1_corrected
rm(met1_1_corrected)
rm(met1_1)

chloroplast_data<-met1_3[seqnames(met1_3)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_3_corrected<-met1_3
met1_3_corrected$readsM<-round((met1_3$readsM)-(met1_3$readsN)*(1-rate))
met1_3_corrected$readsM[met1_3_corrected$readsM<0]<-0
met1_3_corrected$readsN<-round(met1_3$readsN*rate)

met13<-met1_3_corrected
rm(met1_3_corrected)
rm(met1_3)


bare_grange_csv<-read.csv("1k_50_frame.csv")

bare_grange<-GRanges(seqnames = bare_grange_csv[,2], ranges = IRanges(start = as.numeric(bare_grange_csv[,3]), end = as.numeric(bare_grange_csv[,4]), names = bare_grange_csv[,1]), strand = bare_grange_csv[,6])

WT_CG_P<-analyseReadsInsideRegionsForCondition(bare_grange, WT, "CG")
met11_CG_P<-analyseReadsInsideRegionsForCondition(bare_grange, met11, "CG")
met13_CG_P<-analyseReadsInsideRegionsForCondition(bare_grange, met13, "CG")

save(WT_CG_P, file = "WT_CG_P.RData")
save(met11_CG_P, file = "met11_CG_P.RData")
save(met13_CG_P, file = "met13_CG_P.RData")