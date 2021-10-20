#!/usr/bin/Rscript

library(DMRcaller)

epi1<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")

bare_grange_csv<-read.csv("bare_atha_grange.csv")

bare_grange<-GRanges(seqnames = bare_grange_csv[,2], ranges = IRanges(start = as.numeric(bare_grange_csv[,3]), end = as.numeric(bare_grange_csv[,4]), names = bare_grange_csv[,1]), strand = bare_grange_csv[,6])
chloroplast_data<-epi1[seqnames(epi1)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
epi1_corrected<-epi1
epi1_corrected$readsM<-round((epi1$readsM)-(epi1$readsN)*(1-rate))
epi1_corrected$readsM[epi1_corrected$readsM<0]<-0
epi1_corrected$readsN<-round(epi1$readsN*rate)
epi1<-epi1_corrected
out_grange_CHG_A<-analyseReadsInsideRegionsForCondition(bare_grange, epi1, "CHG")

save(out_grange_CHG_A, file = "CHG-Atha.RData")