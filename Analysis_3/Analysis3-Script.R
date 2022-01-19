##
## Run file - Analysis 3
## Requires :
## - athtiling1.0rcdf
## - Gene expression data
## - Data from separate DMRcaller analysis
## - TAIR changelog
##  
## Requires 2 folders: base_data_entry, containing methylation data and other, and expression_calculation with, rcdf file, and subfolder "both" containing expression data
##
## IMPORTANT:
## This analysis *MUST* be run with empty environment
## 

library(tidyverse)
library(ggplot2)
library(gridExtra)
library(DMRcaller)
library(beepr)
library(patchwork)
library(limma)
library(affy)
library(TxDb.Athaliana.BioMart.plantsmart22)
library(viridis)
library(ggExtra)
library(ggVennDiagram)
library(affyPLM)
library(ggvenn)
library(VennDiagram)

#####
## Step 1: Expression data extraction
## Extract expression results from bio-replicates, and process
#####

set.seed(101)
setwd("expression_calculation/")

## Load the rcdf object
load("athtiling1.0rcdf.rda")

## Load all the reads
setwd("both")
## Read into R
allreads_1<-ReadAffy(cdfname = "athtiling1.0rcdf")
## Generate expressionset through RMA, but without normalisation
all_rma_1_nonorm<-rma(allreads_1, normalize = FALSE)
## Perform normalisation
all_rma_all<-affyPLM::normalize.ExpressionSet.loess(all_rma_1_nonorm, transfn = "antilog")
## Extract expression data
expression_all<-data.frame(exprs(all_rma_all))
setwd("../")

## Generate the matrix with contents of each column
design_all <- model.matrix(~ 0+factor(c(1,1,1,2,2,2,3,3,3,4,4,4)))
colnames(design_all)<-c("WT_met11", "met11", "WT_met13", "met13")

## Fit the data to a linear model using the matrix
fit1_all <- lmFit(all_rma_all, design_all)

## Construct contrast matrix with WT1 to met1-1 and WT3 to met1-3
contrast.matrix1_all <- makeContrasts(met11-WT_met11, met13-WT_met13, levels=design_all)

## Fit the data to contrasts
fit2_all <- contrasts.fit(fit1_all, contrast.matrix1_all)
## Perform empirical bayes analysis
fit2_all <- eBayes(fit2_all)

bayes_fit<-fit2_all

setwd("../")

#####
## Step 2: Gene data extraction
## Generate gene and promoter data for separate DMRcaller analysis
#####

## Extract data for all genes
all_genes<-genes(TxDb.Athaliana.BioMart.plantsmart22)

## Identify which genes are present in expression data
analysed_genes<-rownames(fit2_all$coefficients)


setwd("base_data_entry/")
## Load TAIR changelog
Tair_changelog<-read.csv("TAIR10_locushistory.csv.gz", header = FALSE)
## Fixing first line
Tair_changelog[1,1]<-"AT1G00001"
setwd("../")

## Identify genes present in analysed genes, but not present in annotation
different_genes<-setdiff(analysed_genes, names(all_genes))

## Finds where obsolete genes were merged with annotated ones
merge_obsolete<-as.data.frame(tibble(Tair_changelog) %>%
                                filter(V3 == "mergeobsolete"))

## Here, generate replacement gene assignments
replacements<-c() ## Replacement gene assignments 
obsolete_removed<-c() ## Obsolete genes without replacements
replaced<-c() ## Genes that have been replaced
for(i in 1:length(different_genes)){
  grepped_MO<-grep(different_genes[i], merge_obsolete[,1], ignore.case = TRUE)
  if(length(grepped_MO)==1){
    replacement<-merge_obsolete[grepped_MO,5]
    replacements<-c(replacements, replacement)
    replaced<-c(replaced, different_genes[i])
  }
  else{
    obsolete_removed<-c(obsolete_removed, different_genes[i])
  }
}

## Remove spaces present in corrections file
replacements<-gsub(" ", "", replacements)


## Fix analysed_genes (ie: replace genes with replacements)
for(i in 1:length(replaced)){
  analysed_genes[replaced[i]==analysed_genes]<-replacements[i]
}

## Remove duplicate genes
included_dupes<-analysed_genes
analysed_genes<-unique(analysed_genes)


## Extract data for all genes
## Also, remove obsolete genes
all_genes_loop<-c()
for(i in 1:length(analysed_genes)){
  
  all_genes_loop<-c(all_genes_loop,(1:length(all_genes))[names(all_genes)==analysed_genes[i]])
  
}

## Extract location data for genes
new_atha_grange<-(all_genes[all_genes_loop])

## Generate location data for promoters
## 1000 downstream, 50 upstream
prom_DF<-data.frame((new_atha_grange))
prom_DF_plus<-prom_DF[prom_DF[,5]=="+",]
prom_DF_minus<-prom_DF[prom_DF[,5]=="-",]
prom_DF_plus_processing<-data.frame((prom_DF_plus[,2]-1000), (prom_DF_plus[,2]+50))
prom_DF_minus_processing<-data.frame((prom_DF_minus[,3]+1000), (prom_DF_minus[,3]-50))
prom_range_plus<-GRanges(IRanges(start = prom_DF_plus_processing[,1], end = prom_DF_plus_processing[,2], names = (prom_DF_plus[,6])), strand = "+", seqnames = prom_DF_plus[,1])
prom_range_minus<-GRanges(IRanges(start = prom_DF_minus_processing[,2], end = prom_DF_minus_processing[,1], names = (prom_DF_minus[,6])), strand = "-", seqnames = prom_DF_minus[,1])
prom_range_unsorted<-c(prom_range_plus, prom_range_minus)
prom_range_unsorted_names<-names(prom_range_unsorted)

## Sort promoter data, to have the same order as gene data
prom_range_unsorted_names_loop<-c()
for(i in 1:length(analysed_genes)){
  prom_range_unsorted_names_loop<-c(prom_range_unsorted_names_loop, (1:length(prom_range_unsorted_names))[prom_range_unsorted_names==analysed_genes[i]])
}
prom_grange<-prom_range_unsorted[prom_range_unsorted_names_loop]

## Rename the object
bare_atha_grange<-new_atha_grange

## Prepare data frames for export
DF_AG<-data.frame(new_atha_grange)
rownames(DF_AG)<-names(new_atha_grange)

## Export the files for DMRcaller separate analysis
write.csv(as.data.frame(prom_grange), "1k_50_frame.csv")
write.csv(as.data.frame(new_atha_grange), "bare_atha_grange.csv")

#####
## Step 3: Data entry
## Enter data from separate DMRcaller analysis
#####

setwd("base_data_entry/")

load("WT_CG_A.RData")
load("WT_CHG_A.RData")
load("WT_CHH_A.RData")

load("met11_CG_A.RData")
load("met11_CHG_A.RData")
load("met11_CHH_A.RData")

load("met13_CG_A.RData")
load("met13_CHG_A.RData")
load("met13_CHH_A.RData")

load("WT_CG_P.RData")
load("WT_CHG_P.RData")
load("WT_CHH_P.RData")

load("met11_CG_P.RData")
load("met11_CHG_P.RData")
load("met11_CHH_P.RData")

load("met13_CG_P.RData")
load("met13_CHG_P.RData")
load("met13_CHH_P.RData")

load("CG-WT-met1_1.RData")
load("CG-WT-met1_3.RData")
load("CG-met1_1-met1_3.RData")

load("CHG-WT-met1_1.RData")
load("CHG-WT-met1_3.RData")
load("CHG-met1_1-met1_3.RData")

load("CHH-WT-met1_1.RData")
load("CHH-WT-met1_3.RData")
load("CHH-met1_1-met1_3.RData")

setwd("../")

#####
## Step 4: Data coercion
## Ensure the objects are ready to be merged are ready to be merged
#####

## Identify differences between analysed genes and the GRange objects (obsolete and duplicate genes)
differences<-unique(c(setdiff(names(bare_atha_grange), rownames(expression_all)), setdiff(rownames(expression_all), names(bare_atha_grange))))

## Adjust atha_grange
bare_atha_grange_names<-names(bare_atha_grange)
bare_adjustment<-c()
for(i in 1:length(differences)){
  
  bare_adjustment<-c(bare_adjustment, (1:29576)[bare_atha_grange_names==differences[i]])
  
}

adj_bare_atha_grange<-bare_atha_grange[-bare_adjustment]
## Now use this to adjust other data
## bare_adjustment also applies to bound_changes_2

## DMRcaller data:
corrected_names<-names(adj_bare_atha_grange)

unordered_names<-names(WT_CG_A)
ordering_vec_ranges_a<-c()
for(i in corrected_names){
  
  if(sum(unordered_names==i)>0){
    loopnum<-(1:length(unordered_names))[unordered_names==i]
    ordering_vec_ranges_a<-c(ordering_vec_ranges_a, loopnum)
  }
}

unordered_names<-names(WT_CG_P)
ordering_vec_ranges_p<-c()
for(i in corrected_names){
  
  if(sum(unordered_names==i)>0){
    loopnum<-(1:length(unordered_names))[unordered_names==i]
    ordering_vec_ranges_p<-c(ordering_vec_ranges_p, loopnum)
  }
}

## Expression_all:
unordered_names<-rownames(expression_all)
ordering_vec_exp_new<-c()
for(i in corrected_names){
  
  if(sum(unordered_names==i)>0){
    loopnum<-(1:length(unordered_names))[unordered_names==i]
    ordering_vec_exp_new<-c(ordering_vec_exp_new, loopnum)
  }
}


## Ordering
## WARNING: From this point on, objects are coerced to use the same gene list, in the same order.

WT_CG_A<-WT_CG_A[ordering_vec_ranges_a]
WT_CHG_A<-WT_CHG_A[ordering_vec_ranges_a]
WT_CHH_A<-WT_CHH_A[ordering_vec_ranges_a]
met11_CG_A<-met11_CG_A[ordering_vec_ranges_a]
met11_CHG_A<-met11_CHG_A[ordering_vec_ranges_a]
met11_CHH_A<-met11_CHH_A[ordering_vec_ranges_a]
met13_CG_A<-met13_CG_A[ordering_vec_ranges_a]
met13_CHG_A<-met13_CHG_A[ordering_vec_ranges_a]
met13_CHH_A<-met13_CHH_A[ordering_vec_ranges_a]

WT_CG_P<-WT_CG_P[ordering_vec_ranges_p]
WT_CHG_P<-WT_CHG_P[ordering_vec_ranges_p]
WT_CHH_P<-WT_CHH_P[ordering_vec_ranges_p]
met11_CG_P<-met11_CG_P[ordering_vec_ranges_p]
met11_CHG_P<-met11_CHG_P[ordering_vec_ranges_p]
met11_CHH_P<-met11_CHH_P[ordering_vec_ranges_p]
met13_CG_P<-met13_CG_P[ordering_vec_ranges_p]
met13_CHG_P<-met13_CHG_P[ordering_vec_ranges_p]
met13_CHH_P<-met13_CHH_P[ordering_vec_ranges_p]

## bayes_fit uses ordering_vec_exp_new as well, but because it is a list it needs to be applied manually when used
expression_all<-expression_all[ordering_vec_exp_new,]

#####
## Step 5: Methylation change extraction
## Identify presence of differentially methylated regions within genes (generate bound_changes_2 and its prom version)
#####

## For genes
## For each context
for(context in c("CG", "CHG", "CHH")){
  ## For three types
  for(met_type in c("WT_to_met1_1", "WT_to_met1_3", "met1_1_to_met1_3")){
    ## Set the range that's analysed
    input_range<-get(paste0("DMRs", context, "_", met_type))
    ## Find overlaps between DMRs and ranges of a gene, of select size
    overlap_tibble<-tibble(data.frame(findOverlaps(input_range, bare_atha_grange, minoverlap = 50)))
    
    ## Unless DMR is located, assume methylation is maintained
    all_met<-rep("maintained", length(bare_atha_grange))
    
    numeric_names_overlap_table<-as.numeric(names(table(overlap_tibble[,2])))
    
    ## For each overlap:
    for(i in 1:length(table(overlap_tibble[,2]))){
      requested_DF<-(input_range[as.numeric(unlist(overlap_tibble[numeric_names_overlap_table[i] == overlap_tibble[,2],1]))])
      met_num<-unique(as.numeric(unlist(overlap_tibble[numeric_names_overlap_table[i] == overlap_tibble[,2],2])))
      ## If overlapping only gained or lost DMRs, write that
      if(length(unique(mcols(requested_DF)[,11]))==1){
        all_met[met_num]<-unique(mcols(requested_DF)[,11])
      }
      ## If overlapping both gain AND loss DMRs, write "both@
      if(length(unique(mcols(requested_DF)[,11]))!=1){
        all_met[met_num]<-"both"
      }
    }
    ## Assign output
    assign(paste0("bare_ag_", context, "_", met_type), all_met)
    beep()
  }
}

## For promoters
## As before
for(context in c("CG", "CHG", "CHH")){
  for(met_type in c("WT_to_met1_1", "WT_to_met1_3", "met1_1_to_met1_3")){
    input_range<-get(paste0("DMRs", context, "_", met_type))
    overlap_tibble<-tibble(data.frame(findOverlaps(input_range, prom_grange, minoverlap = 50)))
    
    all_met<-rep("maintained", length(prom_grange))
    
    numeric_names_overlap_table<-as.numeric(names(table(overlap_tibble[,2])))
    for(i in 1:length(table(overlap_tibble[,2]))){
      requested_DF<-(input_range[as.numeric(unlist(overlap_tibble[numeric_names_overlap_table[i] == overlap_tibble[,2],1]))])
      met_num<-unique(as.numeric(unlist(overlap_tibble[numeric_names_overlap_table[i] == overlap_tibble[,2],2])))
      if(length(unique(mcols(requested_DF)[,11]))==1){
        all_met[met_num]<-unique(mcols(requested_DF)[,11])
      }
      if(length(unique(mcols(requested_DF)[,11]))!=1){
        all_met[met_num]<-"both"
      }
    }
    
    assign(paste0("prom_", context, "_", met_type), all_met)
    beep()
  }
}

## Merge these objects into a single data frame
## Separate for genes and promoters
bound_changes<-data.frame(bare_ag_CG_WT_to_met1_1, bare_ag_CG_WT_to_met1_3, bare_ag_CG_met1_1_to_met1_3, bare_ag_CHG_WT_to_met1_1, bare_ag_CHG_WT_to_met1_3, bare_ag_CHG_met1_1_to_met1_3, bare_ag_CHH_WT_to_met1_1, bare_ag_CHH_WT_to_met1_3, bare_ag_CHH_met1_1_to_met1_3)
bound_changes_prom<-data.frame(prom_CG_WT_to_met1_1, prom_CG_WT_to_met1_3, prom_CG_met1_1_to_met1_3, prom_CHG_WT_to_met1_1, prom_CHG_WT_to_met1_3, prom_CHG_met1_1_to_met1_3, prom_CHH_WT_to_met1_1, prom_CHH_WT_to_met1_3, prom_CHH_met1_1_to_met1_3)

## Change column names
colnames(bound_changes)<-c("CG_WT_met1_1", "CG_WT_met1_3", "CG_met1_1_met1_3", "CHG_WT_met1_1", "CHG_WT_met1_3", "CHG_met1_1_met1_3", "CHH_WT_met1_1", "CHH_WT_met1_3", "CHH_met1_1_met1_3")
colnames(bound_changes_prom)<-c("CG_WT_met1_1", "CG_WT_met1_3", "CG_met1_1_met1_3", "CHG_WT_met1_1", "CHG_WT_met1_3", "CHG_met1_1_met1_3", "CHH_WT_met1_1", "CHH_WT_met1_3", "CHH_met1_1_met1_3")

## Turn to tibble, and add names
bound_changes_2<-tibble(names(bare_atha_grange), bound_changes)
colnames(bound_changes_2)<-c("genename",colnames(bound_changes_2)[-1])

bound_changes_prom_2<-tibble(names(prom_grange), bound_changes_prom)
colnames(bound_changes_prom_2)<-c("genename",colnames(bound_changes_prom_2)[-1])

## Sound alert that it's finished
## This step is time-consuming
beep()

#####
## Step 6: FDR calculation
## Calculate the FDR for differential expression between WT1 and met1_1 and WT3 and met1_3
#####

bayes_WT_met1_1_FDR<-p.adjust(bayes_fit$p.value[,1], method = "fdr")
bayes_WT_met1_3_FDR<-p.adjust(bayes_fit$p.value[,2], method = "fdr")

#####
## Step 7: Gene type assignment
## Assign each gene its methylation type
#####

met_ID_loop<-c()
for(i in 1:length(WT_CG_A)){
  
  ## Handle errors
  if(unlist(mcols(WT_CHG_A[i])[3])=="NaN" | unlist(mcols(WT_CHH_A[i])[3])=="NaN" | unlist(mcols(WT_CG_A[i])[3])=="NaN"){
    
    met_ID_loop<-c(met_ID_loop, "ERROR:Undefined") ## No mutation error
    
  }
  else{
    ## Check for CHG and CHH methylation
    if(unlist(mcols(WT_CHG_A[i])[3])>0.05 | unlist(mcols(WT_CHH_A[i])[3])>0.05){
      met_ID_loop<-c(met_ID_loop, "Transposable-like")
    }
    else{
      ## Otherwise, either gbM, if proportion is above 0.1 in CG
      if(unlist(mcols(WT_CG_A[i])[3])>0.10){
        met_ID_loop<-c(met_ID_loop, "gbM")
      }
      ## Or non-methylated, if it's equal to or lower than 0.1
      else{
        met_ID_loop<-c(met_ID_loop, "non-methylated")
      }
    }
    
    ## Print a progress-bar into the console, as this is time consiming
    if((i/500)%%1==0){
      print(i/(length(WT_CG_A)))
    }
  }
}

## Generate the data frame
bound_met_type<-data.frame(cbind(names(WT_CG_A[,1]), met_ID_loop))

## Make a sound
beep()

#####
## Step 8: Gene expression calculations
## Calculate the expression mean, SD, and CV
#####

## Calculate mean between bio-reps in their gorups
exp_mean<-function(x){
  return(c(mean(x[1:3]), mean(x[4:6]), mean(x[7:9]), mean(x[10:12])))
}

## Calculate SD
exp_stdev<-function(x){
  return(c(sd(x[1:3]), sd(x[4:6]), sd(x[7:9]), sd(x[10:12])))
}

## calculate CV from mean and SD
mean_expression<-t(apply(expression_all, 1, exp_mean))
sd_expression<-t(apply(expression_all, 1, exp_stdev))
CV_expression<-sd_expression/mean_expression

## Update colnames
colnames(mean_expression)<-c("WT_met1_1", "met1_1", "WT_met1_3", "met1_3")
colnames(CV_expression)<-c("WT_met1_1", "met1_1", "WT_met1_3", "met1_3")

#####
## Step 9: Reads calculation
## Extract how many reads for methylation, and determine if number is high enough for analysis
#####

## WT to met1_1
WT_met1_1_difference_CG_A<-(mcols(met11_CG_A)[,3]-mcols(WT_CG_A)[,3])
WT_met1_1_difference_CHG_A<-(mcols(met11_CHG_A)[,3]-mcols(WT_CHG_A)[,3])
WT_met1_1_difference_CHH_A<-(mcols(met11_CHH_A)[,3]-mcols(WT_CHH_A)[,3])

WT_met1_1_difference_CG_P<-(mcols(met11_CG_P)[,3]-mcols(WT_CG_P)[,3])
WT_met1_1_difference_CHG_P<-(mcols(met11_CHG_P)[,3]-mcols(WT_CHG_P)[,3])
WT_met1_1_difference_CHH_P<-(mcols(met11_CHH_P)[,3]-mcols(WT_CHH_P)[,3])

## Assume that enough reads by default
WT_met1_1_sigdiff_CG_A<-tibble(names(met11_CG_A), rep(">=25_reads", 29575))
WT_met1_1_sigdiff_CHG_A<-tibble(names(met11_CG_A), rep(">=25_reads", 29575))
WT_met1_1_sigdiff_CHH_A<-tibble(names(met11_CG_A), rep(">=25_reads", 29575))
WT_met1_1_sigdiff_CG_P<-tibble(names(met11_CG_P), rep(">=25_reads", 29575))
WT_met1_1_sigdiff_CHG_P<-tibble(names(met11_CG_P), rep(">=25_reads", 29575))
WT_met1_1_sigdiff_CHH_P<-tibble(names(met11_CG_P), rep(">=25_reads", 29575))

## If either WT or mutant has less than 25 reads total, mark it as such
WT_met1_1_sigdiff_CG_A[((mcols(WT_CG_A)[,2]+mcols(WT_CG_A)[,1]) < 25 | (mcols(met11_CG_A)[,2]+mcols(met11_CG_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_1_sigdiff_CHG_A[((mcols(WT_CHG_A)[,2]+mcols(WT_CHG_A)[,1]) < 25 | (mcols(met11_CHG_A)[,2]+mcols(met11_CHG_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_1_sigdiff_CHH_A[((mcols(WT_CHH_A)[,2]+mcols(WT_CHH_A)[,1]) < 25 | (mcols(met11_CHH_A)[,2]+mcols(met11_CHH_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_1_sigdiff_CG_P[((mcols(WT_CG_P)[,2]+mcols(WT_CG_P)[,1]) < 25 | (mcols(met11_CG_P)[,2]+mcols(met11_CG_P)[,1]) < 25),2]<-"<25_reads"
WT_met1_1_sigdiff_CHG_P[((mcols(WT_CHG_P)[,2]+mcols(WT_CHG_P)[,1]) < 25 | (mcols(met11_CHG_P)[,2]+mcols(met11_CHG_P)[,1]) < 25),2]<-"<25_reads"
WT_met1_1_sigdiff_CHH_P[((mcols(WT_CHH_P)[,2]+mcols(WT_CHH_P)[,1]) < 25 | (mcols(met11_CHH_P)[,2]+mcols(met11_CHH_P)[,1]) < 25),2]<-"<25_reads"


## WT to met1_3
WT_met1_3_difference_CG_A<-(mcols(met13_CG_A)[,3]-mcols(WT_CG_A)[,3])
WT_met1_3_difference_CHG_A<-(mcols(met13_CHG_A)[,3]-mcols(WT_CHG_A)[,3])
WT_met1_3_difference_CHH_A<-(mcols(met13_CHH_A)[,3]-mcols(WT_CHH_A)[,3])

WT_met1_3_difference_CG_P<-(mcols(met13_CG_P)[,3]-mcols(WT_CG_P)[,3])
WT_met1_3_difference_CHG_P<-(mcols(met13_CHG_P)[,3]-mcols(WT_CHG_P)[,3])
WT_met1_3_difference_CHH_P<-(mcols(met13_CHH_P)[,3]-mcols(WT_CHH_P)[,3])


WT_met1_3_sigdiff_CG_A<-tibble(names(met13_CG_A), rep(">=25_reads", 29575))
WT_met1_3_sigdiff_CHG_A<-tibble(names(met13_CG_A), rep(">=25_reads", 29575))
WT_met1_3_sigdiff_CHH_A<-tibble(names(met13_CG_A), rep(">=25_reads", 29575))
WT_met1_3_sigdiff_CG_P<-tibble(names(met13_CG_P), rep(">=25_reads", 29575))
WT_met1_3_sigdiff_CHG_P<-tibble(names(met13_CG_P), rep(">=25_reads", 29575))
WT_met1_3_sigdiff_CHH_P<-tibble(names(met13_CG_P), rep(">=25_reads", 29575))

## If either WT or mutant has less than 25 reads total, mark it as such
WT_met1_3_sigdiff_CG_A[((mcols(WT_CG_A)[,2]+mcols(WT_CG_A)[,1]) < 25 | (mcols(met13_CG_A)[,2]+mcols(met13_CG_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_3_sigdiff_CHG_A[((mcols(WT_CHG_A)[,2]+mcols(WT_CHG_A)[,1]) < 25 | (mcols(met13_CHG_A)[,2]+mcols(met13_CHG_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_3_sigdiff_CHH_A[((mcols(WT_CHH_A)[,2]+mcols(WT_CHH_A)[,1]) < 25 | (mcols(met13_CHH_A)[,2]+mcols(met13_CHH_A)[,1]) < 25),2]<-"<25_reads"
WT_met1_3_sigdiff_CG_P[((mcols(WT_CG_P)[,2]+mcols(WT_CG_P)[,1]) < 25 | (mcols(met13_CG_P)[,2]+mcols(met13_CG_P)[,1]) < 25),2]<-"<25_reads"
WT_met1_3_sigdiff_CHG_P[((mcols(WT_CHG_P)[,2]+mcols(WT_CHG_P)[,1]) < 25 | (mcols(met13_CHG_P)[,2]+mcols(met13_CHG_P)[,1]) < 25),2]<-"<25_reads"
WT_met1_3_sigdiff_CHH_P[((mcols(WT_CHH_P)[,2]+mcols(WT_CHH_P)[,1]) < 25 | (mcols(met13_CHH_P)[,2]+mcols(met13_CHH_P)[,1]) < 25),2]<-"<25_reads"

#####
## Step 10: Data merging
## Merge all of the data calculated/extracted before into a pair of objects!
#####

##
## Met1_1
##

all_tibble_met1_1_expression<-tibble(data.frame(
  rownames(mean_expression), ##genenames
  mean_expression[,1:2], ## mean expression
  CV_expression[,1:2], ##CV of expression
  bayes_fit$coefficients[ordering_vec_exp_new,1], ## expression fold change
  bayes_fit$p.value[ordering_vec_exp_new,1], ## p-value
  bayes_WT_met1_1_FDR[ordering_vec_exp_new], ## FDR
  bound_met_type[,2] ##type of methylation
))
colnames(all_tibble_met1_1_expression)<-c("genename", "WT_expression", "met1_1_expression", "WT_expression_CV", "met1_1_expression_CV", "expression_fc_new", "p_value_new", "FDR_new", "methylation_type")

## First genes
## Extract relevant data
component1<-data.frame(all_tibble_met1_1_expression, mcols(WT_CG_A[,3]), mcols(met11_CG_A[,3]), bound_changes_2[-bare_adjustment,c(2)], "CG", "gene", WT_met1_1_difference_CG_A, WT_met1_1_sigdiff_CG_A[,2])
component2<-data.frame(all_tibble_met1_1_expression, mcols(WT_CHG_A[,3]), mcols(met11_CHG_A[,3]), bound_changes_2[-bare_adjustment,c(5)], "CHG", "gene", WT_met1_1_difference_CHG_A, WT_met1_1_sigdiff_CHG_A[,2])
component3<-data.frame(all_tibble_met1_1_expression, mcols(WT_CHH_A[,3]), mcols(met11_CHH_A[,3]), bound_changes_2[-bare_adjustment,c(8)], "CHH", "gene", WT_met1_1_difference_CHH_A, WT_met1_1_sigdiff_CHH_A[,2])
## Ensure colnames match
colnames(component1)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component2)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component3)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")

## Then promoters
component4<-data.frame(all_tibble_met1_1_expression, mcols(WT_CG_P[,3]), mcols(met11_CG_P[,3]), bound_changes_prom_2[-bare_adjustment,c(2)], "CG", "promoter", WT_met1_1_difference_CG_P, WT_met1_1_sigdiff_CG_P[,2])
component5<-data.frame(all_tibble_met1_1_expression, mcols(WT_CHG_P[,3]), mcols(met11_CHG_P[,3]), bound_changes_prom_2[-bare_adjustment,c(5)], "CHG", "promoter", WT_met1_1_difference_CHG_P, WT_met1_1_sigdiff_CHG_P[,2])
component6<-data.frame(all_tibble_met1_1_expression, mcols(WT_CHH_P[,3]), mcols(met11_CHH_P[,3]), bound_changes_prom_2[-bare_adjustment,c(8)], "CHH", "promoter", WT_met1_1_difference_CHH_P, WT_met1_1_sigdiff_CHH_P[,2])
colnames(component4)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component5)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component6)<-c(colnames(all_tibble_met1_1_expression), "proportion_WT", "proportion_met11", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")

## Merge data for genes and promoters
all_tibble_met1_1<-tibble(rbind(component1, component2, component3, component4, component5, component6))
all_tibble_met1_1<-tibble(all_tibble_met1_1, unlist(c(rep(bound_met_type[,2], 3), rep("promoter", (length(component4[,1])*3)))))
colnames(all_tibble_met1_1)<-c(colnames(all_tibble_met1_1)[-17], "methylation_type_loc")

##
##Met1_3
##

all_tibble_met1_3_expression<-tibble(data.frame(
  rownames(mean_expression), ##genenames
  mean_expression[,3:4], ##mean expression
  CV_expression[,3:4], ##CV of expression
  bayes_fit$coefficients[ordering_vec_exp_new,2], ##expression fold change
  bayes_fit$p.value[ordering_vec_exp_new,2], ##p-value
  bayes_WT_met1_3_FDR[ordering_vec_exp_new], ##FDR
  bound_met_type[,2] ##type of methylation
  
))
colnames(all_tibble_met1_3_expression)<-c("genename", "WT_expression", "met1_3_expression", "WT_expression_CV", "met1_3_expression_CV", "expression_fc_new", "p_value_new", "FDR_new", "methylation_type")

## Genes
component1<-data.frame(all_tibble_met1_3_expression, mcols(WT_CG_A[,3]), mcols(met13_CG_A[,3]), bound_changes_2[-bare_adjustment,c(3)], "CG", "gene", WT_met1_3_difference_CG_A, WT_met1_3_sigdiff_CG_A[,2])
component2<-data.frame(all_tibble_met1_3_expression, mcols(WT_CHG_A[,3]), mcols(met13_CHG_A[,3]), bound_changes_2[-bare_adjustment,c(6)], "CHG", "gene", WT_met1_3_difference_CHG_A, WT_met1_3_sigdiff_CHG_A[,2])
component3<-data.frame(all_tibble_met1_3_expression, mcols(WT_CHH_A[,3]), mcols(met13_CHH_A[,3]), bound_changes_2[-bare_adjustment,c(9)], "CHH", "gene", WT_met1_3_difference_CHH_A, WT_met1_3_sigdiff_CHH_A[,2])
colnames(component1)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component2)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component3)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")

## Promoters
component4<-data.frame(all_tibble_met1_3_expression, mcols(WT_CG_P[,3]), mcols(met13_CG_P[,3]), bound_changes_prom_2[-bare_adjustment,c(3)], "CG", "promoter", WT_met1_3_difference_CG_P, WT_met1_3_sigdiff_CG_P[,2])
component5<-data.frame(all_tibble_met1_3_expression, mcols(WT_CHG_P[,3]), mcols(met13_CHG_P[,3]), bound_changes_prom_2[-bare_adjustment,c(6)], "CHG", "promoter", WT_met1_3_difference_CHG_P, WT_met1_3_sigdiff_CHG_P[,2])
component6<-data.frame(all_tibble_met1_3_expression, mcols(WT_CHH_P[,3]), mcols(met13_CHH_P[,3]), bound_changes_prom_2[-bare_adjustment,c(9)], "CHH", "promoter", WT_met1_3_difference_CHH_P, WT_met1_3_sigdiff_CHH_P[,2])
colnames(component4)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component5)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")
colnames(component6)<-c(colnames(all_tibble_met1_3_expression), "proportion_WT", "proportion_met13", "methylation_state", "methylation_context", "methylated_body", "methylation_difference", "methylation_difference_reads")

## Merge data for genes and promoters
all_tibble_met1_3<-tibble(rbind(component1, component2, component3, component4, component5, component6))
all_tibble_met1_3<-tibble(all_tibble_met1_3, unlist(c(rep(bound_met_type[,2], 3), rep("promoter", (length(component4[,1])*3)))))
colnames(all_tibble_met1_3)<-c(colnames(all_tibble_met1_3)[-17], "methylation_type_loc")



## Separate genes without methylation data
all_tibble_met1_1_w_undefined<-all_tibble_met1_1
all_tibble_met1_1<-(all_tibble_met1_1 %>% filter(methylation_type_loc != "ERROR:Undefined"))

all_tibble_met1_3_w_undefined<-all_tibble_met1_3
all_tibble_met1_3<-(all_tibble_met1_3 %>% filter(methylation_type_loc != "ERROR:Undefined"))

## Separage genes with less than 25 reads
all_tibble_met1_1_w_low_reads<-all_tibble_met1_1
all_tibble_met1_3_w_low_reads<-all_tibble_met1_3

all_tibble_met1_1<-(all_tibble_met1_1 %>% filter(methylation_difference_reads==">=25_reads"))
all_tibble_met1_3<-(all_tibble_met1_3 %>% filter(methylation_difference_reads==">=25_reads"))

## Keep only methylation in CG context
all_tibble_met1_1_old<-all_tibble_met1_1
all_tibble_met1_3_old<-all_tibble_met1_3

all_tibble_met1_1<-(all_tibble_met1_1 %>% filter(methylation_context == "CG"))
all_tibble_met1_3<-(all_tibble_met1_3 %>% filter(methylation_context == "CG"))

## Remove genes that overlap with gain DMRs (both "gain" and "both")
## No genes are removed for met1_3
all_tibble_met1_1_w_other<-all_tibble_met1_1
all_tibble_met1_3_w_other<-all_tibble_met1_3

all_tibble_met1_1<-(all_tibble_met1_1 %>% filter(methylation_state == "loss" | methylation_state == "maintained"))
all_tibble_met1_3<-(all_tibble_met1_3 %>% filter(methylation_state == "loss" | methylation_state == "maintained"))

all_tibble_met1_1$methylation_type_loc<-factor(all_tibble_met1_1$methylation_type_loc, levels = c("gbM", "Transposable-like", "non-methylated", "promoter"))
all_tibble_met1_3$methylation_type_loc<-factor(all_tibble_met1_3$methylation_type_loc, levels = c("gbM", "Transposable-like", "non-methylated", "promoter"))

#####
## Step 11: Basic plot
## Draw the plot showing distribution of genes in 3 categories and promoters, expression in WT1 vs methylation proportion
#####

WT_vs_met_plot<-ggplot(all_tibble_met1_1, aes(WT_expression, proportion_WT, color = methylation_state)) + 
  geom_point() +
  facet_grid(cols = vars(methylation_type_loc), rows = vars(methylation_context)) +
  xlab("WT-1 expression") +
  ylab("CG methylation proportion") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_color_manual(values=rev(c("#9C9CA8", "#CA3535"))) +
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=17)) +
  labs(color='DMR overlap') +
  guides(color = guide_legend(override.aes = list(size = 5)))


dir.create("plots")
setwd("plots")
ggsave(WT_vs_met_plot, filename = "WT_vs_met_plot.png", width = 12, height = 5, device = "png", scale = 1.5)
setwd("../")

#####
## Step 12: Expression and methylation change plot
## Draw the plot showing expression and methylation data and their changes
#####

## Determine how the genes change between WT and met1-1
met1_1_list<-rep("Methylation_maintained_expression_changed", (dim(all_tibble_met1_1)[1]))
met1_1_list[all_tibble_met1_1$methylation_state == "loss"]<-"Methylation_loss_expression_changed" ## Because the latter supersedes this, only those where expression changed will stay
met1_1_list[all_tibble_met1_1$methylation_state == "loss" & ((all_tibble_met1_1$FDR_new>0.05) | (all_tibble_met1_1$expression_fc_new>=(-1) & all_tibble_met1_1$expression_fc_new<=(1)))]<-"Methylation_loss_expression_no_change"
met1_1_list[all_tibble_met1_1$methylation_state == "maintained" & ((all_tibble_met1_1$FDR_new>0.05) | (all_tibble_met1_1$expression_fc_new>=(-1) & all_tibble_met1_1$expression_fc_new<=(1)))]<-"Methylation_maintained_expression_no_change"

## Merge data
all_tibble_met1_1_quartered<-tibble(all_tibble_met1_1, met1_1_list)
colnames(all_tibble_met1_1_quartered)<-c(colnames(all_tibble_met1_1_quartered)[-18], "gene_state")

## Determine how the genes change between WT and met1-3
met1_3_list<-rep("Methylation_maintained_expression_changed", (dim(all_tibble_met1_3)[1]))
met1_3_list[all_tibble_met1_3$methylation_state == "loss"]<-"Methylation_loss_expression_changed" ## Because the latter supersedes this, only those where expression changed will stay
met1_3_list[all_tibble_met1_3$methylation_state == "loss" & ((all_tibble_met1_3$FDR_new>0.05) | (all_tibble_met1_3$expression_fc_new>=(-1) & all_tibble_met1_3$expression_fc_new<=(1)))]<-"Methylation_loss_expression_no_change"
met1_3_list[all_tibble_met1_3$methylation_state == "maintained" & ((all_tibble_met1_3$FDR_new>0.05) | (all_tibble_met1_3$expression_fc_new>=(-1) & all_tibble_met1_3$expression_fc_new<=(1)))]<-"Methylation_maintained_expression_no_change"

## Merge data
all_tibble_met1_3_quartered<-tibble(all_tibble_met1_3, met1_3_list)
colnames(all_tibble_met1_3_quartered)<-c(colnames(all_tibble_met1_3_quartered)[-18], "gene_state")


dir.create("plots")
setwd("plots")

## Generate the plot for WT1 and met1-1
met1_1_exp<-ggplot((all_tibble_met1_1_quartered%>% filter(methylated_body == "gene" & methylation_context == "CG")), aes(WT_expression, met1_1_expression, color = gene_state)) + 
  geom_point() +
  xlab("WT-1 expression") + 
  ylab("Met1-1 expression") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_color_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme(legend.position = "none") +
  xlim(c(2.45, 15.5)) +
  ylim(c(2.45, 15.5))

## Generate the plot for WT3 and met1-3
met1_3_exp<-ggplot((all_tibble_met1_3_quartered%>% filter(methylated_body == "gene" & methylation_context == "CG")), aes(WT_expression, met1_3_expression, color = gene_state)) + 
  geom_point() +
  xlab("WT-3 expression") + 
  ylab("Met1-3 expression") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_color_manual(values=c("#D81B60", "#1E88E5", "#FFC107", "#004D40")) +
  theme(legend.position = "none") +
  xlim(c(2.45, 15.5)) +
  ylim(c(2.45, 15.5))

## Merge the two plots
merged_plot<- (met1_1_exp | met1_3_exp)
ggsave(merged_plot, filename = paste0("gene_scatter_and_hist_quartered.png"), width = 16, height = 8, device = "png", scale = 1.15)
setwd("../")

#####
## Step 13: Gene filtering
## Generate object with genes that did not have significant expression change (FDR larger than 0.05, or fold change between -1 and 1)
## Calculate fold change in CV between WT and met1 mutants
## If FC above or below 1, CV either increased or decreased
#####

all_tibble_met1_1_nonsig<-(all_tibble_met1_1 %>% filter((FDR_new>0.05) | (expression_fc_new>=(-1) & expression_fc_new<=(1))))
all_tibble_met1_3_nonsig<-(all_tibble_met1_3 %>% filter((FDR_new>0.05) | (expression_fc_new>=(-1) & expression_fc_new<=(1))))

all_tibble_met1_1_nonsig_vec_1<-log2(all_tibble_met1_1_nonsig$met1_1_expression_CV/all_tibble_met1_1_nonsig$WT_expression_CV)
all_tibble_met1_3_nonsig_vec_1<-log2(all_tibble_met1_3_nonsig$met1_3_expression_CV/all_tibble_met1_3_nonsig$WT_expression_CV)
all_tibble_met1_1_nonsig_vec<-rep("No change", length(all_tibble_met1_1_nonsig_vec_1))
all_tibble_met1_3_nonsig_vec<-rep("No change", length(all_tibble_met1_3_nonsig_vec_1))


all_tibble_met1_1_nonsig_vec[all_tibble_met1_1_nonsig_vec_1>1]<-"increased CV"
all_tibble_met1_1_nonsig_vec[all_tibble_met1_1_nonsig_vec_1<(-1)]<-"decreased CV"

all_tibble_met1_3_nonsig_vec[all_tibble_met1_3_nonsig_vec_1>1]<-"increased CV"
all_tibble_met1_3_nonsig_vec[all_tibble_met1_3_nonsig_vec_1<(-1)]<-"decreased CV"

all_tibble_met1_1_nonsig_wchange<-tibble(all_tibble_met1_1_nonsig, all_tibble_met1_1_nonsig_vec)
all_tibble_met1_3_nonsig_wchange<-tibble(all_tibble_met1_3_nonsig, all_tibble_met1_3_nonsig_vec)

colnames(all_tibble_met1_1_nonsig_wchange)<-c(colnames(all_tibble_met1_1_nonsig_wchange[-18]), "CV_change_fold")
colnames(all_tibble_met1_3_nonsig_wchange)<-c(colnames(all_tibble_met1_3_nonsig_wchange[-18]), "CV_change_fold")

#####
## Step 14: CV change plots
## Draw plots depicting CV change
#####

dir.create("plots/")
setwd("plots/")

## Draw plot between WT1 and met1-1
plot_wt_met11_CV_nonsig_nomain<-ggplot((all_tibble_met1_1_nonsig_wchange %>% filter(methylation_state != "maintained")), aes(WT_expression_CV, met1_1_expression_CV, color = CV_change_fold)) + 
  geom_point() +
  facet_grid(cols = vars(methylation_type_loc), rows = vars(methylation_context)) +
  xlab("WT-1 CV of expression") +
  ylab("Met1_1 CV of expression") +
  geom_abline(intercept = 0, slope = 1, size = 3, alpha = 0.5) +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_color_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535"))) +
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=17)) +
  labs(color='CV change') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlim(c(0,0.4))+
  ylim(c(0,0.4))

## Draw plot between WT3 and met1-3
plot_wt_met13_CV_nonsig_nomain<-ggplot((all_tibble_met1_3_nonsig_wchange %>% filter(methylation_state != "maintained")), aes(WT_expression_CV, met1_3_expression_CV, color = CV_change_fold)) + 
  geom_point() +
  facet_grid(cols = vars(methylation_type_loc), rows = vars(methylation_context)) +
  xlab("WT-3 CV of expression") +
  ylab("Met1_3 CV of expression") +
  geom_abline(intercept = 0, slope = 1, size = 3, alpha = 0.5) +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_color_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535"))) +
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=17)) +
  labs(color='CV change') +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  xlim(c(0,0.4))+
  ylim(c(0,0.4))

## Save the plots individually
# ggsave(plot_wt_met11_CV_nonsig_nomain, filename = "plot_wt_met11_CV_nonsig_nomain.png", width = 12, height = 5, device = "png", scale = 1.5)
# ggsave(plot_wt_met13_CV_nonsig_nomain, filename = "plot_wt_met13_CV_nonsig_nomain.png", width = 12, height = 5, device = "png", scale = 1.5)

## But also merge them, and save together.
merged_nonsig_nomain<-(plot_wt_met11_CV_nonsig_nomain/plot_wt_met13_CV_nonsig_nomain  + plot_layout(guides = 'collect'))
ggsave(merged_nonsig_nomain, filename = "merged_CV_change_scatter.png", width = 16, height = 8, device = "png", scale = 1.15)


setwd("../")

#####
## Step 15: CV change barplots
## Draw plots depicting size of groups shown on previous scatterplot
#####

dir.create("plots/")
setwd("plots/")

## Count how many genes in each group for both comparisons
count_at1_s1<-all_tibble_met1_1_nonsig_wchange %>% count(CV_change_fold, methylation_state, methylation_type_loc)
count_at3_s1<-all_tibble_met1_3_nonsig_wchange %>% count(CV_change_fold, methylation_state, methylation_type_loc)

## Plot for WT1-met1_1
splitplot_1<-ggplot((count_at1_s1 %>% filter(methylation_state!="maintained")), aes(methylation_state, n, fill = CV_change_fold)) + 
  geom_col(position = "dodge") +
  facet_grid(cols = vars(methylation_type_loc)) +
  xlab("Overlap with DMRs met1_1") +
  ylab("Count") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_fill_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535"))) +
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=17)) +
  labs(fill='CV change') +
  guides(fill = guide_legend(override.aes = list(size = 5)))

## Plot for WT3-met1_3
splitplot_3<-ggplot((count_at3_s1 %>% filter(methylation_state!="maintained")), aes(methylation_state, n, fill = CV_change_fold)) + 
  geom_col(position = "dodge") +
  facet_grid(cols = vars(methylation_type_loc)) +
  xlab("Overlap with DMRs met1_3") +
  ylab("Count") +
  theme(axis.text.x = element_text(size = 18), 
        axis.text.y = element_text(size = 19),
        text = element_text(size = 17)) +
  scale_fill_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535"))) +
  theme(strip.text.x = element_text(size = 15, face = "bold"),
        strip.text.y = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=17)) +
  labs(fill='CV change') +
  guides(fill = guide_legend(override.aes = list(size = 5)))

## Merge the plots and save
splitplot_m<- splitplot_1 / splitplot_3 + plot_layout(guides = 'collect')
ggsave(splitplot_m, filename = "merged_CV_change_barplot.png", width = 16, height = 8, device = "png", scale = 1.15)

setwd("../")

#####
## Step 16: CV change venn diagrams
## Draws venn diagrams showing overlaps between gene groups
## They have to be merged together manually.
## Also draws the heatmap comparing the groups
#####

dir.create("venndiag")
setwd("venndiag")
for(CV_change in unique(all_tibble_met1_1_nonsig_wchange$CV_change_fold)){
  for(typeloc in unique(all_tibble_met1_1_nonsig_wchange$methylation_type_loc)){
    ## Extract genes in the category
    met1_1_vloop <- unlist((all_tibble_met1_1_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1]) 
    met1_3_vloop <- unlist((all_tibble_met1_3_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1])
    
    
    ## Draw the venn diagram
    venn.diagram(
      x = list(met1_1_vloop, met1_3_vloop),
      category.names = c("Met 1_1" , "Met 1_3"),
      filename = paste0(CV_change, "_", typeloc, ".png"),
      output=TRUE,
      
      imagetype = "png", height = 400, width = 400,
      resolution = 500, compression = "lzw",
      
      lwd = 2, lty = 'blank', fill = c("#f99595", "#95abf9"),
      
      cex = .3, fontface = "bold", fontfamily = "sans",
      
      cat.cex = 0.3, cat.fontface = "plain", cat.default.pos = "outer", 
      cat.pos = c(-2, 2), cat.dist = c(0.04, 0.04), cat.fontfamily = "sans"
    )
  }
  
}

## Draw the heatmap
for(typeloc in c("gbM", "Transposable-like", "non-methylated", "promoter")){
  loop2<-c()
  ## Extract data for each CV change
  for(CV_change in c("decreased CV", "increased CV", "No change")){
    met1_1_vloop <- unlist((all_tibble_met1_1_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1]) 
    met1_3_vloop <- unlist((all_tibble_met1_3_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1])
    CV_change_loop<-cbind(length(intersect(met1_1_vloop, met1_3_vloop)), sum(length(met1_1_vloop), length(met1_3_vloop))-2*length(intersect(met1_1_vloop, met1_3_vloop)))
    loop2<-rbind(loop2, CV_change_loop)
  }
  
  ## Calculate relationships using fisher's test
  yloop<-c()
  for(y in 1:3){
    xloop<-c()
    for(x in 1:3){
      xloop<-c(xloop, fisher.test(rbind(loop2[y,], loop2[x,]), alternative="t")$p.value)
    }
    yloop<-rbind(yloop, xloop)
  }
  rownames(yloop)<-c("decreased CV", "increased CV", "no CV change")
  colnames(yloop)<-c("decreased CV", "increased CV", "no CV change")
  
  ## Prepare data for the plot
  tloop<-tibble(c(rep("decreased CV", 3), rep("increased CV", 3), rep("no CV change", 3)) ,rep(c("decreased CV", "increased CV", "no CV change"), 3), c(yloop[,1], yloop[,2], yloop[,3]), log(c(yloop[,1], yloop[,2], yloop[,3]), 10))
  colnames(tloop)<-c("ro", "co", "val", "logval")
  
  ## Round the data
  roundloop<-c()
  for(x in tloop$val){
    if(x < 0.01){
      roundloop<-c(roundloop,(10^ceiling(log10(x))))
    }
    else{
      roundloop<-c(roundloop, (round(x, 4)))
    }
  }
  
  tloop<-tibble(tloop, roundloop)
  colnames(tloop)<-c("ro", "co", "val", "logval", "rounded_val")
  ## Ensure proper order
  tloop[,1]<-factor(unlist(tloop[,1]), levels = c('no CV change', 'increased CV', 'decreased CV'))
  ## Remove repetitious data
  tloop<-tloop[-c(4,7,8),]
  ## Draw the sub-plot
  assign(paste0("loopplot_",typeloc), ggplot(as_tibble(tloop), aes(ro, co, fill = logval))+
           geom_tile() + 
           scale_x_discrete(position = "top")+
           geom_text(aes(label = rounded_val), size = 8, colour = c("black", "black", "white", "black", "white", "black")) +
           scale_fill_gradient(low = "blue", high = "white") +
           labs(fill = "1e") +
           theme(axis.title=element_blank(), axis.text = element_text(size = 17)))
  
}
## Merge sub-plots into one and save
annoplot<-(`loopplot_gbM` | `loopplot_Transposable-like` | `loopplot_non-methylated` | loopplot_promoter)
ggsave(annoplot, filename = "annoplot.png", width = 10, height = 2, device = "png", scale = 3)

setwd("../")

#####
## Step 17: GO analysis data export
## Export data for GO analysis of increased CV overlaps of gbM genes, to be done with PANTHER
#####

typeloc<-"gbM"
CV_change<-"increased CV"
met1_1_genes_GO <- unlist((all_tibble_met1_1_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1]) 
met1_3_genes_GO <- unlist((all_tibble_met1_3_nonsig_wchange %>% filter(CV_change_fold == CV_change & methylation_state == "loss") %>% filter(methylation_type_loc==typeloc))[,1])

## Generate background from all analysed genes for met1_1
background<-unlist((all_tibble_met1_1 %>% filter(methylated_body == "gene"))[,1])
## Identify the overlap
inc_overlap<-intersect(met1_1_genes_GO, met1_3_genes_GO)

dir.create("GO")
setwd("GO")
## Save the txt files
write(background, "background.txt")
write(inc_overlap, "increased_CV_overlap.txt")

setwd("../")

#####
## Step 18: Base analysis
## Compute low resolution methylation profiles, coverage and correlation
#####

## Load data
setwd("raw_met_data/WT")
WT<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../Met1_1")
met1_1<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../Met1_3")
met1_3<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../")

## Correct data for chloroplast methylation

chloroplast_data<-WT[seqnames(WT)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
WT_corrected<-WT
WT_corrected$readsM<-round((WT$readsM)-(WT$readsN)*(1-rate))
WT_corrected$readsM[WT_corrected$readsM<0]<-0
WT_corrected$readsN<-round(WT$readsN*rate)

chloroplast_data<-met1_1[seqnames(met1_1)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_1_corrected<-met1_1
met1_1_corrected$readsM<-round((met1_1$readsM)-(met1_1$readsN)*(1-rate))
met1_1_corrected$readsM[met1_1_corrected$readsM<0]<-0
met1_1_corrected$readsN<-round(met1_1$readsN*rate)


chloroplast_data<-met1_3[seqnames(met1_3)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
met1_3_corrected<-met1_3
met1_3_corrected$readsM<-round((met1_3$readsM)-(met1_3$readsN)*(1-rate))
met1_3_corrected$readsM[met1_3_corrected$readsM<0]<-0
met1_3_corrected$readsN<-round(met1_3$readsN*rate)

## Remove the space-consuming cx report files
rm(WT, met1_1, met1_3)

## Make chromosome 1 object
regions_loop<-GRanges(seqnames = Rle(1), ranges = IRanges(1,30427670))

###########################################################
#####
#### CG

## Extract data from 3 bio-replicates. Low resolution profile, coverage, and spatial correlation 
profile_1_corr_CG <- computeMethylationProfile(WT_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_1_corr_CG <-computeMethylationDataCoverage(WT_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_1_corr_CG <-computeMethylationDataSpatialCorrelation(WT_corrected,
                                                             context="CG",
                                                             distances=c(1,10,100,1000,10000))

profile_2_corr_CG <- computeMethylationProfile(met1_1_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_2_corr_CG <-computeMethylationDataCoverage(met1_1_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_2_corr_CG <-computeMethylationDataSpatialCorrelation(met1_1_corrected,
                                                             context="CG",
                                                             distances=c(1,10,100,1000,10000))


profile_3_corr_CG <- computeMethylationProfile(met1_3_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_3_corr_CG <-computeMethylationDataCoverage(met1_3_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_3_corr_CG <-computeMethylationDataSpatialCorrelation(met1_3_corrected,
                                                             context="CG",
                                                             distances=c(1,10,100,1000,10000))

## Process the data to make the plots
profiledata_CG<- rbind(cbind(round((start(profile_1_corr_CG)+end(profile_1_corr_CG))/2), mcols(profile_1_corr_CG)[,3], "biorep-1"),
                       cbind(round((start(profile_2_corr_CG)+end(profile_2_corr_CG))/2), mcols(profile_2_corr_CG)[,3], "biorep-2"),
                       cbind(round((start(profile_3_corr_CG)+end(profile_3_corr_CG))/2), mcols(profile_3_corr_CG)[,3], "biorep-3"))
colnames(profiledata_CG)<-c("location_bp", "proportion", "population")
profiledata_CG<-as_tibble(profiledata_CG)
profiledata_CG$location_bp<-as.numeric(profiledata_CG$location_bp)
profiledata_CG$proportion<-as.numeric(profiledata_CG$proportion)
profiledata_CG$population<-factor(profiledata_CG$population, levels = c("biorep-1", "biorep-2", "biorep-3"))


coveragedata<- rbind(cbind(c(1,5,10,15, 20, 25), coverage_1_corr_CG,"biorep-1"),
                     cbind(c(1,5,10,15, 20, 25), coverage_2_corr_CG,"biorep-2"),
                     cbind(c(1,5,10,15, 20, 25), coverage_3_corr_CG, "biorep-3"))
colnames(coveragedata)<-c("min_reads", "coverage", "population")
coveragedata<-as_tibble(coveragedata)

coveragedata$min_reads<-as.numeric(coveragedata$min_reads)
coveragedata$coverage<-as.numeric(coveragedata$coverage)
coveragedata$population<-factor(coveragedata$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

correlationdata_CG<- rbind(cbind(names(spacorr_1_corr_CG), spacorr_1_corr_CG, "biorep-1"),
                           cbind(names(spacorr_2_corr_CG), spacorr_2_corr_CG,"biorep-2"),
                           cbind(names(spacorr_3_corr_CG), spacorr_3_corr_CG, "biorep-3"))
colnames(correlationdata_CG)<-c("distance", "correlation", "population")
correlationdata_CG<-as_tibble(correlationdata_CG)

correlationdata_CG$distance<-as.numeric(correlationdata_CG$distance)
correlationdata_CG$correlation<-as.numeric(correlationdata_CG$correlation)
correlationdata_CG$population<-factor(correlationdata_CG$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

## Generate the 3 plots

profileplot_CG<-ggplot(profiledata_CG, aes(location_bp, proportion, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  xlab("location (bp)") +
  ggtitle('Methylation profile (CG)')

coverageplot<-ggplot(coveragedata, aes(min_reads, coverage, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  xlab("minimum number of reads") +
  ggtitle('Coverage')

correlationplot_CG<-ggplot(correlationdata_CG, aes(distance, correlation, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10') +
  ggtitle('Correlation (CG)')

## Remove unnecessary files
rm(profile_1_corr_CG, profile_2_corr_CG, profile_3_corr_CG, coverage_1_corr_CG, coverage_2_corr_CG, coverage_3_corr_CG, spacorr_1_corr_CG, spacorr_2_corr_CG, spacorr_3_corr_CG,
   profiledata_CG, correlationdata_CG, coveragedata)
## Reclaim memory
gc()
####
#####
###########################################################
###########################################################
#####
#### CHG

## As before, except no spatial coverage

profile_1_corr_CHG <- computeMethylationProfile(WT_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_1_corr_CHG <-computeMethylationDataSpatialCorrelation(WT_corrected,
                                                              context="CHG",
                                                              distances=c(1,10,100,1000,10000))

profile_2_corr_CHG <- computeMethylationProfile(met1_1_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_2_corr_CHG <-computeMethylationDataSpatialCorrelation(met1_1_corrected,
                                                              context="CHG",
                                                              distances=c(1,10,100,1000,10000))


profile_3_corr_CHG <- computeMethylationProfile(met1_3_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_3_corr_CHG <-computeMethylationDataSpatialCorrelation(met1_3_corrected,
                                                              context="CHG",
                                                              distances=c(1,10,100,1000,10000))

profiledata_CHG<- rbind(cbind(round((start(profile_1_corr_CHG)+end(profile_1_corr_CHG))/2), mcols(profile_1_corr_CHG)[,3], "biorep-1"),
                        cbind(round((start(profile_2_corr_CHG)+end(profile_2_corr_CHG))/2), mcols(profile_2_corr_CHG)[,3], "biorep-2"),
                        cbind(round((start(profile_3_corr_CHG)+end(profile_3_corr_CHG))/2), mcols(profile_3_corr_CHG)[,3], "biorep-3"))
colnames(profiledata_CHG)<-c("location_bp", "proportion", "population")
profiledata_CHG<-as_tibble(profiledata_CHG)
profiledata_CHG$location_bp<-as.numeric(profiledata_CHG$location_bp)
profiledata_CHG$proportion<-as.numeric(profiledata_CHG$proportion)
profiledata_CHG$population<-factor(profiledata_CHG$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

correlationdata_CHG<- rbind(cbind(names(spacorr_1_corr_CHG), spacorr_1_corr_CHG, "biorep-1"),
                            cbind(names(spacorr_2_corr_CHG), spacorr_2_corr_CHG,"biorep-2"),
                            cbind(names(spacorr_3_corr_CHG), spacorr_3_corr_CHG, "biorep-3"))
colnames(correlationdata_CHG)<-c("distance", "correlation", "population")
correlationdata_CHG<-as_tibble(correlationdata_CHG)

correlationdata_CHG$distance<-as.numeric(correlationdata_CHG$distance)
correlationdata_CHG$correlation<-as.numeric(correlationdata_CHG$correlation)
correlationdata_CHG$population<-factor(correlationdata_CHG$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

profileplot_CHG<-ggplot(profiledata_CHG, aes(location_bp, proportion, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  xlab("location (bp)") +
  ggtitle('Methylation profile (CHG)')

correlationplot_CHG<-ggplot(correlationdata_CHG, aes(distance, correlation, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10') +
  ggtitle('Correlation (CHG)')

rm(profile_1_corr_CHG, profile_2_corr_CHG, profile_3_corr_CHG, spacorr_1_corr_CHG, spacorr_2_corr_CHG, spacorr_3_corr_CHG,
   profiledata_CHG, correlationdata_CHG)
gc()
####
#####
###########################################################
###########################################################
#####
#### CHH
profile_1_corr_CHH <- computeMethylationProfile(WT_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_1_corr_CHH <-computeMethylationDataSpatialCorrelation(WT_corrected,
                                                              context="CHH",
                                                              distances=c(1,10,100,1000,10000))

profile_2_corr_CHH <- computeMethylationProfile(met1_1_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_2_corr_CHH <-computeMethylationDataSpatialCorrelation(met1_1_corrected,
                                                              context="CHH",
                                                              distances=c(1,10,100,1000,10000))


profile_3_corr_CHH <- computeMethylationProfile(met1_3_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_3_corr_CHH <-computeMethylationDataSpatialCorrelation(met1_3_corrected,
                                                              context="CHH",
                                                              distances=c(1,10,100,1000,10000))

profiledata_CHH<- rbind(cbind(round((start(profile_1_corr_CHH)+end(profile_1_corr_CHH))/2), mcols(profile_1_corr_CHH)[,3], "biorep-1"),
                        cbind(round((start(profile_2_corr_CHH)+end(profile_2_corr_CHH))/2), mcols(profile_2_corr_CHH)[,3], "biorep-2"),
                        cbind(round((start(profile_3_corr_CHH)+end(profile_3_corr_CHH))/2), mcols(profile_3_corr_CHH)[,3], "biorep-3"))
colnames(profiledata_CHH)<-c("location_bp", "proportion", "population")
profiledata_CHH<-as_tibble(profiledata_CHH)
profiledata_CHH$location_bp<-as.numeric(profiledata_CHH$location_bp)
profiledata_CHH$proportion<-as.numeric(profiledata_CHH$proportion)
profiledata_CHH$population<-factor(profiledata_CHH$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

correlationdata_CHH<- rbind(cbind(names(spacorr_1_corr_CHH), spacorr_1_corr_CHH, "biorep-1"),
                            cbind(names(spacorr_2_corr_CHH), spacorr_2_corr_CHH,"biorep-2"),
                            cbind(names(spacorr_3_corr_CHH), spacorr_3_corr_CHH, "biorep-3"))
colnames(correlationdata_CHH)<-c("distance", "correlation", "population")
correlationdata_CHH<-as_tibble(correlationdata_CHH)

correlationdata_CHH$distance<-as.numeric(correlationdata_CHH$distance)
correlationdata_CHH$correlation<-as.numeric(correlationdata_CHH$correlation)
correlationdata_CHH$population<-factor(correlationdata_CHH$population, levels = c("biorep-1", "biorep-2", "biorep-3"))

profileplot_CHH<-ggplot(profiledata_CHH, aes(location_bp, proportion, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  xlab("location (bp)") +
  ggtitle('Methylation profile (CHH)')

correlationplot_CHH<-ggplot(correlationdata_CHH, aes(distance, correlation, color = population)) + 
  geom_point() +
  geom_line() +
  theme_bw() +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values = viridis(10)[c(1,5,8)]) +
  theme(legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4))) +
  ylim(c(0,1)) +
  scale_x_continuous(trans='log10') +
  ggtitle('Correlation (CHH)')

rm(profile_1_corr_CHH, profile_2_corr_CHH, profile_3_corr_CHH, spacorr_1_corr_CHH, spacorr_2_corr_CHH, spacorr_3_corr_CHH,
   profiledata_CHH, correlationdata_CHH)
gc()

rm(met1_1_corrected, met1_3_corrected, WT_corrected)
gc()
dir.create("met_profiles")
setwd("met_profiles")
## Merge the plots into a single figure
allplot<- ((profileplot_CG | profileplot_CHG) / (profileplot_CHH | (correlationplot_CG | correlationplot_CHG) / (correlationplot_CHH | coverageplot))) + plot_layout(guides = 'collect')
#allplot<-(((coverageplot | correlationplot_CG) / (correlationplot_CHG | correlationplot_CHH)) | profileplot_CG) / (profileplot_CHG | profileplot_CHH) + plot_layout(guides = 'collect')
ggsave(allplot, filename = "base_analysis_plots.png", width = 10, height = 10, device = "png", scale = 1.5)
setwd("../")

####
#####
###########################################################
