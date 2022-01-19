##
## Full run file - Analysis 1
## Requires pre-generated drought and hili mock files
## Requires 1 folder underneath wd: input, with drought, high.light, assignment data, tair changelog, and colour palette
##


library(tidyverse)
library(beepr)
library(limma)
library(TxDb.Athaliana.BioMart.plantsmart22)
library(ggplot2)
library(ggforce)
library(VennDiagram)
library(patchwork)
library(ggdendro)
library(reshape2)
library(grid)
library(viridis)
library(gridExtra)
library(GOfuncR)
library(GO.db)

#### 1. Data input, and GRange formation

## Step 1: Read prepared expression data, probe-to-gene assignments, and TAIR10 changelog
setwd("./input")

drought<-read.table("drought_mock.txt.gz")
high.light<-read.table("highlight_mock.txt.gz")
assignment_data<-read.csv("summary_by_probe.csv.gz")
Tair_changelog<-read.csv("TAIR10_locushistory.csv.gz", header = FALSE)

## First line of changelog file is incorrectly read - manually correcting here
Tair_changelog[1,1]<-"AT1G00001"

setwd("../")
dir.create("output")
setwd("output/")

## Step 2: Match high light and drought probes
## Order drought and high light
drought<-drought[order(drought[,1]),]
high.light<-high.light[order(high.light[,1]),]

## High light probeset is longer and contains all drought probes. 
## This gets coordinates of all drought probes within high light.
dro_hili_matches<-c()
for(i in 1:length(high.light[,1])){
  grepvar<-grep(high.light[i,1], drought[,1], ignore.case = TRUE)
  if(length(grepvar)!=0){
    dro_hili_matches<-c(dro_hili_matches, i)
  }
  if((i/100)%%1==0){
    print(i/325.78)
  }
}

## This merges matching probes between high light and drought. 
bound_experiments<-cbind(high.light[dro_hili_matches,], drought[-1])

## Here, merged data is normalized using cyclicloess method.
normalized_data<-data.frame(bound_experiments[,1], normalizeBetweenArrays(bound_experiments[,-1], method = "cyclicloess"))

## Step 3: Match each probe to a gene
## Create a GRange containing gene information to assign later
all_genes<-genes(TxDb.Athaliana.BioMart.plantsmart22)

## Change both assignment data and normalized data names to lowercase
## While grep can be instructed to ignore case, here it must be done this way
lowercase_normalizeddata_probes<-(tolower(normalized_data[,1]))
lowercase_assignmentdata_probes<-(tolower(assignment_data[,1]))

## Here we assign each CATMA probe a gene
## 3 possible outputs: Assigned gene, if single linked gene, OR
## "No linked genes", if no genes linked to probe, OR "Multiple linked genes", if probe linked to more than one gene.
non_assigned<-c() ## Debug: non_assigned should stay empty. If not empty, something went wrong!
assigned<-c() ## output
for(i in 1:length(normalized_data[,1])){
  grepped<-assignment_data[grep(paste("^", lowercase_normalizeddata_probes[i], "$", sep=""), lowercase_assignmentdata_probes),]
  if(length(grepped[,1])!=0){
    if(grepped[,3]==""){
      grepped[,3]<-"No linked genes"
    }
    if(grepped[,4]!=""){
      grepped[,3]<-"Multiple linked genes"
    } 
    output<-grepped[c(1,3)]
    assigned<-rbind(assigned, output)
  }
  else{
    non_assigned<-c(non_assigned, lowercase_normalizeddata_probes[i])
  }
}

## Step 4: Assign each probe-gene pair location
## Attempts to build a dataset to form the GRange object.
## However, a number of genes are absent from genetic location annotation
## This code identifies these genes
no_gene<-c() #Absent genes
GRange_compatible_data_tst<-c()
for(i in 1:length(normalized_data[,1])){
  gene_info<-as.data.frame(all_genes[names(all_genes)==assigned[i,2]])
  if(length(gene_info[,1])!=0){
    output<-data.frame(normalized_data[i,1], gene_info[6], gene_info[-6], normalized_data[i,-1])
    GRange_compatible_data_tst<-rbind(GRange_compatible_data_tst, output)
  }
  else{
    if(assigned[i,2]!="Multiple linked genes" & assigned[i,2]!="No linked genes"){
      no_gene<-c(no_gene, assigned[i,2])
    }
  }
  if((i/100)%%1==0){
    print(i/325.01)
  }
}

## Finds where obsolete genes were merged with annotated ones
merge_obsolete<-as.data.frame(tibble(Tair_changelog) %>%
                                filter(V3 == "mergeobsolete"))

## Here, generate replacement gene assignments
replacements<-c() ## Replacement gene assignments 
obsolete_removed<-c() ## Obsolete genes without replacements
replaced<-c() ## Genes that have been replaced
for(i in 1:length(no_gene)){
  grepped_MO<-grep(no_gene[i], merge_obsolete[,1], ignore.case = TRUE)
  if(length(grepped_MO)==1){
    replacement<-merge_obsolete[grepped_MO,5]
    replacements<-c(replacements, replacement)
    replaced<-c(replaced, no_gene[i])
  }
  else{
    obsolete_removed<-c(obsolete_removed, no_gene[i])
  }
}

## For purposes of Sankey plot, duplicates assigned
assigned2<-assigned

## Replaces obsolete assigned genes
for(i in 1:length(replaced)){
  rep_grep<-grep(replaced[i], assigned2[,2], ignore.case = TRUE)
  assigned2[rep_grep,2]<-replacements[i]
}

## Calculates data for the GRange object, now with replaced obsolete genes
no_gene2<-c() ## Debug: This should match obsolete_removed.
GRange_compatible_data<-c()
for(i in 1:length(normalized_data[,1])){
  gene_info<-as.data.frame(all_genes[names(all_genes)==assigned2[i,2]])
  if(length(gene_info[,1])!=0){
    output<-data.frame(normalized_data[i,1], gene_info[6], gene_info[-6], normalized_data[i,-1])
    GRange_compatible_data<-rbind(GRange_compatible_data, output)
  }
  else{
    if(assigned2[i,2]!="Multiple linked genes" & assigned2[i,2]!="No linked genes"){
      no_gene2<-c(no_gene2, assigned[i,2])
    }
  }
  if((i/100)%%1==0){
    print(i/325.01)
  }
}

colnames(GRange_compatible_data)<-c("catma_ID", colnames(GRange_compatible_data)[-1])


## Step 5: For multiple probes linked to a single gene, calculate mean measurements for each biorep
## Orders GRange data
ordered_GRange_data<-GRange_compatible_data[order(GRange_compatible_data[,2]),]

## Creates two loops - one containing single probe-gene pair numbers, other containing probes linked to multiple genes
single_gene_loop<-1
multi_gene_loop<-10
for(i in ((1:length(ordered_GRange_data[,2]))[-c(1,10,21363)])){
  if((!ordered_GRange_data[(i-1),2]==ordered_GRange_data[i,2])&(!ordered_GRange_data[(i+1),2]==ordered_GRange_data[i,2])){
    single_gene_loop<-c(single_gene_loop, i)
  }
  else{
    multi_gene_loop<-c(multi_gene_loop, i)
  }
}

## Manually checked: 1 is single-gene, 21363 is single-gene, 10 is multi-gene
## Must specify these manually
single_gene_loop<-c(single_gene_loop, 21363)

## Prepares the data for the loop
multi_gene_data<-ordered_GRange_data[multi_gene_loop,]
multivar<-c()
multi_data<-c()

## Calculate mean expression of gene for biorep from multiple probes
for(i in 2:length(multi_gene_data[,1])){
  
  ## First: If gene is followed by a repeat, note it, and keep going 
  if(multi_gene_data[i,2]==multi_gene_data[(i-1),2]){
    multivar<-c(multivar, i, (i-1))
  }
  ## Second: If next gene in list is a different gene, sum up the notes, and proceed
  else{
    multivar<-unique(multivar)
    multivar_data<-multi_gene_data[multivar,-(1:7)]
    multivar_processed_data<-c()
    ## This for-loop runs for each biorep, and calculates mean values
    for(o in 1:length(multivar_data[1,])){
      multivar_colvar<-mean(multivar_data[,o])
      multivar_processed_data<-cbind(multivar_processed_data, multivar_colvar)
    }
    ## Prepares output
    multivar_processed_data<-cbind(multi_gene_data[(i-1),1:7], multivar_processed_data)
    multi_data<-rbind(multi_data, multivar_processed_data)
    
    ## Resets multivar in preparation for next loop
    multivar<-c()
  }
}

## Same as before - however, for the last gene, it must be run again outside the for-loop
multivar<-unique(multivar)
multivar_data<-multi_gene_data[multivar,-(1:7)]
multivar_processed_data<-c()
for(o in 1:length(multivar_data[1,])){
  multivar_colvar<-mean(multivar_data[,o])
  multivar_processed_data<-cbind(multivar_processed_data, multivar_colvar)
}
multivar_processed_data<-cbind(multi_gene_data[(i-1),1:7], multivar_processed_data)
multi_data<-rbind(multi_data, multivar_processed_data)
multivar<-c()

## Collects single probe-gene data into an object...
single_gene_data<-ordered_GRange_data[single_gene_loop,]
## ... And reuses its colnames on multi_data
colnames(multi_data)<-colnames(single_gene_data)
## Binds single- and multi-gene data into a single object, to pass on to the GRange
GRange_compatible_data_corrected<-rbind(single_gene_data, multi_data)
## Re-orders the data.frame after rejoining single and multi-gene data
GRange_compatible_data_corrected<-GRange_compatible_data_corrected[order(GRange_compatible_data_corrected[,2]),]

## Step 6: Genomic Ranges formation
## Form a GRange using previously prepared data - only genes are submitted first.
atha_grange <- GRanges(
  seqnames = GRange_compatible_data_corrected$seqnames,
  ranges = IRanges(start = as.numeric(GRange_compatible_data_corrected$start), end = as.numeric(GRange_compatible_data_corrected$end), names = GRange_compatible_data_corrected$gene_id),
  strand = GRange_compatible_data_corrected$strand)

## Adds metadata columns to the GRange, containing TAIR gene ID (strictly speaking unnecessary, as it is also present in the main GRange body - however, it helps with debugging, and making sure
## that rows are aligned properly), as well as normalized expression data which has been prepared in previous steps
mcols(atha_grange)<-GRange_compatible_data_corrected[,-c(1,3,4,5,6,7)]

## Saves the Genomic Ranges object, so that previous time-intensive steps don't have to be re-run in the event of an environment restart
save(atha_grange, file = "atha_grange.RData")

## Create a second GRange object, without metadata, for methylation analysis
bare_atha_grange <- GRanges(
  seqnames = GRange_compatible_data_corrected$seqnames,
  ranges = IRanges(start = as.numeric(GRange_compatible_data_corrected$start), end = as.numeric(GRange_compatible_data_corrected$end), names = GRange_compatible_data_corrected$gene_id),
  strand = GRange_compatible_data_corrected$strand)

## Save the secondary GRange object
save(bare_atha_grange, file = "bare_atha_grange.RData")

## If ran on a device with an audio output, indicate this step is finished.
beep("ping")



#### 2. PCA - Principal component analysis
##
## 1. PCA - hili

## In case this was run without running part 1 first, load atha_grange 
load("atha_grange.RData")

## Load colour palette
setwd("../input")
col_palette<-as.data.frame(read_table("palette.txt", col_names = FALSE)[,3])
setwd("../output")

## Separate metadata from the GRange
atha_grange_meta<-mcols(atha_grange)
## Prepare the colour palette
## Palette: The 15-colour palette from http://mkweb.bcgsc.ca/colorblind
col_palette<-paste0("#",unlist(col_palette))
col_palette2<-col_palette[c(TRUE, FALSE)]

## Separates gene names from the data, and transposes it
dro_and_hili<-atha_grange_meta[,-1]
transpo_data<-t(as.data.frame.array(dro_and_hili))

## Separates hili from dro
thili_hili<-transpo_data[1:52,]
tdro_dro<-transpo_data[53:108,]


## Hili PCA
## Writing down hili timepoints for further use
tpoints_hili<-rep(c("0H","0.5H","1H","1.5H",'2H','2.5H','3H','3.5H','4H','4.5H','5H','5.5H','6H'),4)

## Calculates data for the PCA
data.PC_hili = prcomp(thili_hili,scale.=TRUE)

## Prepares point designations (for the plot) - first all, then specifically chosen points.
point_designations_hili<-paste0(c(rep("A", 13), rep("B", 13), rep("C", 13), rep("D", 13)), "", rep(c("0","0.5","1","1.5",'2','2.5','3','3.5','4','4.5','5','5.5','6'),4))
point_designations_2_hili<-c("A_0",rep("",15), "B_1.5", rep("",14), "C_2.5", "", "C_3.5", rep("",18)) 

## Forms a tibble, to pass on to ggplot
PC1_PC2_tibble_hili<-tibble(tpoints_hili, c(rep("A", 13), rep("B", 13), rep("C", 13), rep("D", 13)), point_designations_hili, point_designations_2_hili, data.frame(unlist(data.PC_hili$x[,1:2])))
colnames(PC1_PC2_tibble_hili)<-c("time", "Biorep", "ID", "ID_2", "PC1", "PC2")

## Makes time into a factor again, so that the plot is ordered
PC1_PC2_tibble_hili$time<-factor(PC1_PC2_tibble_hili$time, levels = unique(PC1_PC2_tibble_hili$time))

## Output plot: PCA hili
## ggplot-ing the PCA data
hili_PCA_plot<-ggplot(PC1_PC2_tibble_hili, aes(x = PC1, y = PC2, colour = time, group = time)) +
  geom_point(size=4) +
  geom_text(label=PC1_PC2_tibble_hili$ID_2, nudge_x = -1, nudge_y = 7, check_overlap = TRUE, colour = "black") +
  scale_color_manual(values=col_palette2[c(1:4, 6:9, 11:15)])+
  theme(axis.text.x = element_text(size = 16), 
        legend.key.height = unit(1, "cm"),
        axis.text.y = element_text(size = 16),
        text = element_text(size = 16))


## Drought PCA
## Unless specified otherwise, same as hili, only using drought data
tpoints_dro<-rep(c("0D","1D","2D","3D",'4D','5D','6D','7D','8D','9D','10D','11D','12D', '13D'),4)

data.PC_dro = prcomp(tdro_dro, scale.=TRUE)

## Different timepoints chosen to be labelled
point_designations_dro<-paste0(c(rep("a", 14), rep("b", 14), rep("c", 14), rep("d", 14)), "", rep(c("0","1","2","3",'4','5','6','7','8','9','10','11','12', '13'),4))
point_designations_2_dro<-c(rep("", 10), "A_10", rep("", 20), "C_3", rep("", 6),  "C_10", "C_11", "C_12", rep("",11), "D_10", "", "", "D_13")

PC1_PC2_tibble_dro<-tibble(tpoints_dro, c(rep("A", 14), rep("B", 14), rep("C", 14), rep("D", 14)), point_designations_dro, point_designations_2_dro,  data.frame(unlist(data.PC_dro$x[,1:2])))
colnames(PC1_PC2_tibble_dro)<-c("time", "Biorep", "ID", "ID_2", "PC1", "PC2")

PC1_PC2_tibble_dro$time<-factor(PC1_PC2_tibble_dro$time, levels = unique(PC1_PC2_tibble_dro$time))


## Output plot: PCA drought
## ggplot-ing the drought data
dro_pca_plot<-ggplot(PC1_PC2_tibble_dro, aes(x = PC1, y = PC2, colour = time)) +
  geom_point(size=4) +
  geom_text(label=PC1_PC2_tibble_dro$ID_2, nudge_x = -1, nudge_y = 7, check_overlap = FALSE, colour = "black") +
  scale_color_manual(values=c(col_palette2[c(1:4, 6:14)], "black"))+
  theme(axis.text.x = element_text(size = 16), 
        legend.key.height = unit(1, "cm"),
        axis.text.y = element_text(size = 16),
        text = element_text(size = 16))

## Save the two plots. WD should still be output.
ggsave(dro_pca_plot, filename = "Drought_PCA.png", width = 10, height = 10, device = "png", scale = 0.9)
ggsave(hili_PCA_plot, filename = "Hili_PCA.png", width = 10, height = 10, device = "png", scale = 0.9)

## Indicates part is completed.
beep("ping")


#### 3. Calculations
## In this part, variables will be calculated for each time point, as well as for the entire two experiments.
## Of these, 3 are calculated in a loop (var, sd and mean), while the remainder are derivative of these three.
## 

## As in 2, load GRanges object in case environment was cleared
load("atha_grange.RData")
## Create an object to hold the metadata values
atha_grange_meta<-mcols(atha_grange)


## Variance calculation.
## atha_variances object is split into 5 because of performance impact of concentrating large dataframes
atha_variances<-c()
atha_variances1<-c()
atha_variances2<-c()
atha_variances3<-c()
atha_variances4<-c()
atha_variances5<-c()
for(i in 1:length(atha_grange)){
  
  ## Extract data pertaining to a gene
  pergene<-unlist(atha_grange_meta[i,-1])
  genename<-unlist(atha_grange_meta[i,1])
  
  ## Calculate var within bioreps for hili.
  ## for-loop is broken up in order to manually exclude chosen bioreps (here 1.5B and 3.5C)
  hili<-c()
  for(z in 0:2){
    hili<-c(hili, var(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, var(pergene[c((1+3), (27+3), (40+3))]))
  for(z in 4:6){
    hili<-c(hili, var(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, var(pergene[c((1+7), (14+7), (40+7))]))
  for(z in 8:12){
    hili<-c(hili, var(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  
  ## Calculate var within bioreps for drought.
  ## for-loop is broken up in order to manually exclude chosen bioreps (here 3C, 11C, 12C, 13D)
  drought<-c()
  for(z in 0:2){
    drought<-c(drought, var(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, var(pergene[c((53+3), (67+3), (95+3))]))
  for(z in 4:10){
    drought<-c(drought, var(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, var(pergene[c((53+11), (67+11), (95+11))]))
  drought<-c(drought, var(pergene[c((53+12), (67+12), (95+12))]))
  drought<-c(drought, var(pergene[c((53+13), (67+13), (81+13))]))
  
  ## Calculate var for all high light (excluding selected bioreps)
  all_hili<-var(pergene[c(1:16, 18:33, 35:52)])
  ## Calculate var for all drought (excluding selected bioreps)
  all_drought<-var(pergene[c(53:83, 85:91, 94:107)])
  
  ## Prepare Output
  outputan<-data.frame(genename, t(c(hili, drought, all_drought, all_hili)))
  ## Append colnames
  colnames(outputan)<-c("genename", "Hili0HPS", "Hili0.5HPS", "Hili1HPS", "Hili1.5HPS", "Hili2HPS", "Hili2.5HPS", "Hili3HPS", "Hili3.5HPS",
                        "Hili4HPS", "Hili4.5HPS", "Hili5HPS", "Hili5.5HPS", "Hili6HPS", "Dro0DPS", "Dro1DPS", "Dro2DPS", "Dro3DPS", "Dro4DPS",
                        "Dro5DPS", "Dro6DPS", "Dro7DPS", "Dro8DPS", "Dro9DPS", "Dro10DPS", "Dro11DPS", "Dro12DPS", "Dro13DPS",
                        "DM", "HLM")
  ## Create output variables
  if(i<3846){
    atha_variances1<-rbind(atha_variances1, outputan)
  }
  else if(i<7693.6){
    atha_variances2<-rbind(atha_variances2, outputan)
  }
  else if(i<11540.4){
    atha_variances3<-rbind(atha_variances3, outputan)
  }
  else if(i<15387.2){
    atha_variances4<-rbind(atha_variances4, outputan)
  }
  else{
    atha_variances5<-rbind(atha_variances5, outputan)
  }
  
  if((i/100)%%1==0){
    print(i/192.34)
  }
}
## Bind outputs together into single object
atha_variances<-rbind(atha_variances1, atha_variances2, atha_variances3, atha_variances4, atha_variances5)
## Erase unnecessary objects
rm(atha_variances1, atha_variances2, atha_variances3, atha_variances4, atha_variances5)


## Standard Deviation calculations
## As per variability, only sd instead of var
atha_SDs<-c()
atha_SDs1<-c()
atha_SDs2<-c()
atha_SDs3<-c()
atha_SDs4<-c()
atha_SDs5<-c()
for(i in 1:length(atha_grange)){
  pergene<-unlist(atha_grange_meta[i,-1])
  genename<-unlist(atha_grange_meta[i,1])
  
  pse<-c()
  for(z in 0:12){
    pse<-c(pse, sd(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  
  hili<-c()
  for(z in 0:2){
    hili<-c(hili, sd(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, sd(pergene[c((1+3), (27+3), (40+3))]))
  for(z in 4:6){
    hili<-c(hili, sd(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, sd(pergene[c((1+7), (14+7), (40+7))]))
  for(z in 8:12){
    hili<-c(hili, sd(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  
  drought<-c()
  for(z in 0:2){
    drought<-c(drought, sd(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, sd(pergene[c((53+3), (67+3), (95+3))]))
  for(z in 4:10){
    drought<-c(drought, sd(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, sd(pergene[c((53+11), (67+11), (95+11))]))
  drought<-c(drought, sd(pergene[c((53+12), (67+12), (95+12))]))
  drought<-c(drought, sd(pergene[c((53+13), (67+13), (81+13))]))
  
  
  all_hili<-sd(pergene[c(1:16, 18:33, 35:52)])
  
  all_drought<-sd(pergene[c(53:83, 85:91, 94:107)])
  
  outputan<-data.frame(genename, t(c(hili, drought, all_drought, all_hili)))
  colnames(outputan)<-c("genename", "Hili0HPS", "Hili0.5HPS", "Hili1HPS", "Hili1.5HPS", "Hili2HPS", "Hili2.5HPS", "Hili3HPS", "Hili3.5HPS",
                        "Hili4HPS", "Hili4.5HPS", "Hili5HPS", "Hili5.5HPS", "Hili6HPS", "Dro0DPS", "Dro1DPS", "Dro2DPS", "Dro3DPS", "Dro4DPS",
                        "Dro5DPS", "Dro6DPS", "Dro7DPS", "Dro8DPS", "Dro9DPS", "Dro10DPS", "Dro11DPS", "Dro12DPS", "Dro13DPS",
                        "DM", "HLM")
  if(i<3846){
    atha_SDs1<-rbind(atha_SDs1, outputan)
  }
  else if(i<7693.6){
    atha_SDs2<-rbind(atha_SDs2, outputan)
  }
  else if(i<11540.4){
    atha_SDs3<-rbind(atha_SDs3, outputan)
  }
  else if(i<15387.2){
    atha_SDs4<-rbind(atha_SDs4, outputan)
  }
  else{
    atha_SDs5<-rbind(atha_SDs5, outputan)
  }
  
  if((i/100)%%1==0){
    print(i/192.34)
  }
}
atha_SDs<-rbind(atha_SDs1, atha_SDs2, atha_SDs3, atha_SDs4, atha_SDs5)
rm(atha_SDs1, atha_SDs2, atha_SDs3, atha_SDs4, atha_SDs5)

## Mean calculations
## As per variability, only mean instead of var
atha_means<-c()
atha_means1<-c()
atha_means2<-c()
atha_means3<-c()
atha_means4<-c()
atha_means5<-c()
for(i in 1:length(atha_grange)){
  pergene<-unlist(atha_grange_meta[i,-1])
  genename<-unlist(atha_grange_meta[i,1])
  
  pse<-c()
  for(z in 0:12){
    pse<-c(pse, mean(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  
  hili<-c()
  for(z in 0:2){
    hili<-c(hili, mean(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, mean(pergene[c((1+3), (27+3), (40+3))]))
  for(z in 4:6){
    hili<-c(hili, mean(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  hili<-c(hili, mean(pergene[c((1+7), (14+7), (40+7))]))
  for(z in 8:12){
    hili<-c(hili, mean(pergene[c((1+z), (14+z), (27+z), (40+z))]))
  }
  
  drought<-c()
  for(z in 0:2){
    drought<-c(drought, mean(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, mean(pergene[c((53+3), (67+3), (95+3))]))
  for(z in 4:10){
    drought<-c(drought, mean(pergene[c((53+z), (67+z), (81+z), (95+z))]))
  }
  drought<-c(drought, mean(pergene[c((53+11), (67+11), (95+11))]))
  drought<-c(drought, mean(pergene[c((53+12), (67+12), (95+12))]))
  drought<-c(drought, mean(pergene[c((53+13), (67+13), (81+13))]))
  
  
  all_hili<-mean(pergene[c(1:16, 18:33, 35:52)])
  
  all_drought<-mean(pergene[c(53:83, 85:91, 94:107)])
  
  outputan<-data.frame(genename, t(c(hili, drought, all_drought, all_hili)))
  colnames(outputan)<-c("genename", "Hili0HPS", "Hili0.5HPS", "Hili1HPS", "Hili1.5HPS", "Hili2HPS", "Hili2.5HPS", "Hili3HPS", "Hili3.5HPS",
                        "Hili4HPS", "Hili4.5HPS", "Hili5HPS", "Hili5.5HPS", "Hili6HPS", "Dro0DPS", "Dro1DPS", "Dro2DPS", "Dro3DPS", "Dro4DPS",
                        "Dro5DPS", "Dro6DPS", "Dro7DPS", "Dro8DPS", "Dro9DPS", "Dro10DPS", "Dro11DPS", "Dro12DPS", "Dro13DPS",
                        "DM", "HLM")
  if(i<3846){
    atha_means1<-rbind(atha_means1, outputan)
  }
  else if(i<7693.6){
    atha_means2<-rbind(atha_means2, outputan)
  }
  else if(i<11540.4){
    atha_means3<-rbind(atha_means3, outputan)
  }
  else if(i<15387.2){
    atha_means4<-rbind(atha_means4, outputan)
  }
  else{
    atha_means5<-rbind(atha_means5, outputan)
  }
  
  if((i/100)%%1==0){
    print(i/192.34)
  }
}
atha_means<-rbind(atha_means1, atha_means2, atha_means3, atha_means4, atha_means5)
rm(atha_means1, atha_means2, atha_means3, atha_means4, atha_means5)



## Calculate Coefficient of Variability
atha_CV<-cbind(atha_SDs[,1], atha_SDs[,-1]/atha_means[,-1])
colnames(atha_CV)<-c("genename", colnames(atha_CV)[-1])
## Calculate Fano Factor
atha_fano<-cbind(atha_SDs[,1], atha_variances[,-1]/atha_means[,-1])
colnames(atha_fano)<-c("genename", colnames(atha_fano)[-1])
## Calculate Coefficient of Variability Squared
atha_CV2<-cbind(atha_SDs[,1], (atha_CV[,-1])^2)
colnames(atha_CV2)<-c("genename", colnames(atha_CV2)[-1])


## Calculate Distance to Median


## [From this point on, it's imported from package scran, including comments]
DM <- function(mean, cv2, win.size=51) 
  # Computes the distance to median for the CV2 values across all genes, 
  # after fitting an abundance-dependent trend.
  # 
  # written by Jong Kyoung Kim
  # with modifications by Aaron Lun
  # created 12 March 2015 
{
  keep <- mean > 0 & !is.na(cv2) & cv2 > 0
  mean.expr <- log10(mean[keep])
  cv2.expr <- log10(cv2[keep])
  
  o <- order(mean.expr)
  if (win.size%%2L==0L) {
    win.size <- win.size+1L
  }
  med.trend <- runmed(cv2.expr[o], k=win.size)
  med.trend[o] <- med.trend
  
  dm.out <- cv2.expr - med.trend
  DM <- rep(NA_real_, length(keep))
  DM[keep] <- dm.out
  names(DM) <- names(mean)
  DM
}

### [/scran] From this point, it's original code again


## Calculate DMs
atha_DMs<-atha_CV2[,1]
for(i in 2:length(atha_CV2[1,])){
  col_DM<-DM(mean = atha_means[,i], cv2 = atha_CV2[,i])
  atha_DMs<-data.frame(atha_DMs, col_DM)
}
## Manually assign DMs colnames
colnames(atha_DMs)<-c("genename", "Hili0HPS", "Hili0.5HPS", "Hili1HPS", "Hili1.5HPS", "Hili2HPS", "Hili2.5HPS", "Hili3HPS", "Hili3.5HPS",
                        "Hili4HPS", "Hili4.5HPS", "Hili5HPS", "Hili5.5HPS", "Hili6HPS", "Dro0DPS", "Dro1DPS", "Dro2DPS", "Dro3DPS", "Dro4DPS",
                        "Dro5DPS", "Dro6DPS", "Dro7DPS", "Dro8DPS", "Dro9DPS", "Dro10DPS", "Dro11DPS", "Dro12DPS", "Dro13DPS",
                        "DM", "HLM")

## Set rownames to have genenames (stored in first column)
rownames(atha_means)<-atha_means[,1]
rownames(atha_variances)<-atha_variances[,1]
rownames(atha_SDs)<-atha_SDs[,1]
rownames(atha_CV)<-atha_CV[,1]
rownames(atha_fano)<-atha_fano[,1]
rownames(atha_CV2)<-atha_CV2[,1]
rownames(atha_DMs)<-atha_DMs[,1]

## Create tibbles from the above
ttibble2_means<-tibble((atha_means))
ttibble2_variances<-tibble((atha_variances))
ttibble2_SDs<-tibble((atha_SDs))
ttibble2_CV<-tibble((atha_CV))
ttibble2_fano<-tibble((atha_fano))
ttibble2_CV2<-tibble((atha_CV2))
ttibble2_DMs<-tibble((atha_DMs))

beep("ping")



#### 4. Plots - heatmaps
## In this part, heatmaps, venn diagram, and the distribution diagram are generated.
## 

####
## Patchwork
## For the following 4:
## Calculate p-value for wilcoxon rank sum test between them, and round to 6 digits
## First create a pivoted tibble
## Then draw a histogram with vertical line, text stating p-value, and value of the cutoff line


## CV
CV_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_CV[,29]), unlist(ttibble2_CV[,30]))[3]))
pivoted_CV<-pivot_longer(ttibble2_CV[,c(1,29, 30)], c(DM, HLM))
plot_CV<-ggplot(pivoted_CV, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  geom_vline(xintercept = 0.04, "dashed")+
  ggtitle('Coefficient of Variation') +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", CV_wilcox), x=0.25,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.04", x=0.25,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## CV squared
CV2_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_CV[,29]), unlist(ttibble2_CV[,30]))[3]))
pivoted_CV2<-pivot_longer(ttibble2_CV2[,c(1,29, 30)], c(DM, HLM))
plot_CV2<-ggplot(pivoted_CV2, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  ggtitle('CV2') +
  geom_vline(xintercept = 0.0016, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", CV2_wilcox), x=0.13,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.0016", x=0.13,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Fano factor
fano_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_fano[,29]), unlist(ttibble2_fano[,30]))[3]))
pivoted_fano<-pivot_longer(ttibble2_fano[,c(1,29, 30)], c(DM, HLM))
plot_fano<-ggplot(pivoted_fano, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  ggtitle('Fano Factor') +
  geom_vline(xintercept = 0.017, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", fano_wilcox), x=0.13,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.017", x=0.13,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Distance to Median
DM_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_DMs[,29]), unlist(ttibble2_DMs[,30]))[3]))
pivoted_DM<-pivot_longer(ttibble2_DMs[,c(1,29, 30)], c(DM, HLM))
plot_DM<-ggplot(pivoted_DM, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  ggtitle('Distance to Median') +
  geom_vline(xintercept = 0.43, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", DM_wilcox), x=0.55,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.043", x=0.55,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Join the 4 plots into one...
Dro_Hili_Plot <- (plot_CV | plot_CV2) / (plot_fano | plot_DM) + plot_layout(guides = 'collect')
## And export it
ggsave(Dro_Hili_Plot, filename = "Dro_Hili_mass_Plot.png", units = "mm", width = 150, height = 75, scale = 2.75)


## Make the venn diagram
venn.diagram(
  x = list(unlist((ttibble2_CV %>% filter(HLM > 0.04))[,1]), unlist((ttibble2_CV %>% filter(DM > 0.04))[,1])),
  category.names = c("HLM" , "DM"),
  filename = "venn_diagram_HLM_DM.png",
  output=TRUE,
  
  imagetype = "png", height = 1000, width = 1000,
  resolution = 750, compression = "lzw",
  
  lwd = 2, lty = 'blank', fill = c("#f99595", "#95abf9"),
  
  cex = .3, fontface = "plain", fontfamily = "sans",
  
  cat.cex = 0.3, cat.fontface = "plain", cat.default.pos = "outer", 
  cat.pos = c(-2, 2), cat.dist = c(0.04, 0.04), cat.fontfamily = "sans"
)



## Draw and save the mean plot
means_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_means[,29]), unlist(ttibble2_means[,30]))[3]))
pivoted_means<-pivot_longer(ttibble2_means[,c(1,29, 30)], c(DM, HLM))
plot_means<-ggplot(pivoted_means, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  geom_vline(xintercept = 8.4, "dashed")+
  geom_vline(xintercept = 12.05, "dashed") +
  ggtitle('Mean of expression') +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", CV_wilcox), x=0.25,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 8.4         h = 12.05", x=0.25,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14),
        text = element_text(size = 14))
ggsave(plot_means, filename = "plot_means.png", units = "mm", width = 150, height = 75, scale = 2)

## Calculate Quantile values for means (for later GO analysis)
quantile(unlist(ttibble2_means[,29]), c(0.1, 0.25, 0.5, 0.75, 0.9, 1))
quantile(unlist(ttibble2_means[,30]), c(0.1, 0.25, 0.5, 0.75, 0.9, 1))


## Draw and save the plot comparing CV and means for DM and HLM
DM_comparison<-tibble(ttibble2_means[,29], ttibble2_CV[,29], .name_repair = "unique")
HLM_comparison<-tibble(ttibble2_means[,30], ttibble2_CV[,30], .name_repair = "unique")
names(DM_comparison)<-c("mean", "CV")
names(HLM_comparison)<-c("mean", "CV")

DM_comparison<- DM_comparison %>% mutate(`CV value`=cut(CV,breaks=c(-1, 0.01, 0.02, 0.03, 0.04, 0.05,
                                                                    0.06, 0.07, 0.08, 0.09, 0.10, 0.20, max(CV,na.rm=T)),
                                                        labels=c("0-0.01","0.01-0.02","0.02-0.03","0.03-0.04",
                                                                 "0.04-0.05","0.05-0.06","0.06-0.07","0.07-0.08","0.08-0.09",
                                                                 "0.09-0.10", "0.10-0.20", ">0.20"))) %>%
  mutate(`CV value`=factor(as.character(`CV value`),levels=rev(levels(`CV value`))))

HLM_comparison<- HLM_comparison %>% mutate(`CV value`=cut(CV,breaks=c(-1, 0.01, 0.02, 0.03, 0.04, 0.05,
                                                                      0.06, 0.07, 0.08, 0.09, 0.10, 0.20, max(CV,na.rm=T)),
                                                          labels=c("0-0.01","0.01-0.02","0.02-0.03","0.03-0.04",
                                                                   "0.04-0.05","0.05-0.06","0.06-0.07","0.07-0.08","0.08-0.09",
                                                                   "0.09-0.10", "0.10-0.20", ">0.20"))) %>%
  mutate(`CV value`=factor(as.character(`CV value`),levels=rev(levels(`CV value`))))

DM_comparison_plot<-ggplot(DM_comparison, aes(mean, CV, colour = `CV value`)) +
  geom_point(size=1) +
  scale_colour_manual(values=rev(viridis(20))[c(2:8,14,16,18,20)]) +
  xlim(c(7.5, 16.5)) +
  ylim(c(0, 0.18)) +
  geom_hline(yintercept = 0.04, "dashed") +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        text = element_text(size = 15))
HLM_comparison_plot<-ggplot(HLM_comparison, aes(mean, CV, colour = `CV value`)) +
  geom_point(size=1) +
  scale_colour_manual(values=rev(viridis(20))[c(2:8,14,16,18,20)]) +
  xlim(c(7.5, 16.5)) +
  ylim(c(0, 0.18)) +
  geom_hline(yintercept = 0.04, "dashed") +
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16),
        text = element_text(size = 15))
complot<-(HLM_comparison_plot | DM_comparison_plot)  + plot_layout(guides = 'collect')
ggsave(complot, filename = "ComparisonPlot.png", device = "png", units = "mm", width = 150, height = 70, scale = 2.75)


##
## Heatmaps
##

######
##CV

## Select genes matching criteria: Only high light bioreps, CV Higher than 0.04 for high light...
ttible2_0.04_Hili_alone<-ttibble2_CV %>% filter(HLM > 0.04) %>% arrange((HLM)) %>% dplyr:::select(everything()[c(1:14, 29, 30)])
## Here higher than 0.04 for both high light and drought, arranged by high light
ttible2_0.04_Hili<-ttibble2_CV %>% filter(HLM > 0.04) %>% filter(DM > 0.04) %>% arrange((HLM)) %>% dplyr::select(everything()[c(1:14, 29, 30)])

## Higher than 0.04 for drought, only drought bioreps, arranged by drought
ttible2_0.04_Dro_alone<-ttibble2_CV %>% filter(DM > 0.04) %>% arrange((DM)) %>% dplyr::select(everything()[c(1, 15:30)])
## Higher than 0.04 for drought and high light, only drought bioreps, arranged by drought
ttible2_0.04_Dro<-ttibble2_CV %>% filter(HLM > 0.04) %>% filter(DM > 0.04) %>% arrange((DM)) %>% dplyr::select(everything()[c(1, 15:30)])

## Higher than 0.04 for drought and high light, BOTH bioreps, arranged by drought
ttible2_0.04_Both<-ttibble2_CV %>% filter(HLM > 0.04) %>% filter(DM > 0.04) %>% arrange((DM))

## All drought and all highlight, unarranged
ttible2_totaldro<-ttibble2_CV %>% dplyr::select(everything()[c(1, 15:30)]) %>% arrange((DM))
ttible2_totalhili<-ttibble2_CV %>% dplyr::select(everything()[c(1:14, 29, 30)]) %>% arrange((HLM))

## Both, arranged
ttible2_arranged<-ttibble2_CV %>% arrange((DM))

## For-loop making the plots. For each of the tibbles...
## "Thin" plot, either Dro or Hili, but not both
for(i in c("ttible2_0.04_Hili_alone", 'ttible2_0.04_Dro_alone', 'ttible2_totaldro', 'ttible2_totalhili')){
  ## Generate a dendrogram plot
  dendro_plot<-ggdendrogram(data = as.dendrogram(hclust(d = dist(x = t(as.data.frame(get(i)[,-1]))))), labels = FALSE)
  ## Use dendrogram data to order columns (so that when plot is made, columns will be ordered like in the dendogram)
  dendro_ordered<-order.dendrogram(as.dendrogram(hclust(d = dist(x = t(as.data.frame(get(i)[,-1]))))))
  dendro_ordered<-c(0, dendro_ordered)
  dendro_ordered<-dendro_ordered+1
  ordered_i<-get(i)[,dendro_ordered]
  ## "Melt" selected tibble (read: change from wide to long)
  molten_i<-melt(ordered_i)
  
  ## Assign discrete values for the heatmap
  molten_tibble<-tibble(molten_i) %>%
    mutate(variable=factor(variable,levels=(sort(unique(variable))))) %>%
    mutate(value_2=cut(value,breaks=c(-1, 0.01, 0.02, 0.03, 0.04, 0.05,
                                      0.06, 0.07, 0.08, 0.09, 0.10, 0.20, max(value,na.rm=T)),
                       labels=c("0-0.01","0.01-0.02","0.02-0.03","0.03-0.04",
                                "0.04-0.05","0.05-0.06","0.06-0.07","0.07-0.08","0.08-0.09",
                                "0.09-0.10", "0.10-0.20", ">0.20"))) %>%
    mutate(value_2=factor(as.character(value_2),levels=rev(levels(value_2))))
  ## Coerce column names to reflect their contents
  colnames(molten_tibble)<-c("genename", "variable", "value", "CV")
  ## Ensure rows appear in the proper order
  molten_tibble<-data.frame(factor(molten_tibble$genename, levels = unique(molten_tibble$genename)), molten_tibble[,-1])
  ## Ensure columns appear in the proper order
  molten_tibble[,2]<-factor(x = molten_tibble$variable,
                            levels = colnames(ordered_i), 
                            ordered = TRUE)
  ## Coerce column names to reflect their contents again.
  colnames(molten_tibble)<-c("genename", "variable", "value", "CV")
  ## Generate the plot.
  assign(paste0(i,"_test_plot"), ggplot(molten_tibble, aes(variable, genename, fill= CV)) + 
           geom_tile() +
           guides(colour = guide_colourbar(barheight = unit(30, "cm"))) +
           theme(axis.text.x = element_text(size = 25, angle = 270), 
                 legend.key.height = unit(2.5, "cm"),
                 axis.text.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank(),
                 text = element_text(size = 19),
                 axis.title.x = element_text(size = 20),
                 legend.text = element_text(size = 19),
                 legend.title = element_text(size = 19)) +
           xlab("Measurement") +
           scale_fill_manual(values=rev(viridis(20))[c(1:8,14,16,18,20)]))
  
  ## Start the graphics device
  ## Here, the plot is "thin", meaning it is fit for either drought or highlight, but not both, as reflected in the tibbles chosen.
  png(filename = paste0(i, "_CV_thin.png"), width = 1000, height = 1080, units = "px")
  ## Generate new grid page
  grid.newpage()
  ## Align dendrogram...
  print(dendro_plot, vp = viewport(x = 0.430, y = 0.9, width = 0.89, height = 0.15))
  ## And the heatmap.
  ## These x and y values were entered by hand - changing number of columns can result in dendrogram not aligned to heatmap
  print(get(paste0(i,"_test_plot")), vp = viewport(x = 0.5, y = 0.45, width = 1.0, height = 0.87))
  ## Shut off graphics device, save plot
  dev.off()
}

## As before, unless stated otherwise
## "Thick" plot, both Dro and Hili
for(i in c('ttible2_0.04_Both', 'ttible2_arranged')){
  dendro_plot<-ggdendrogram(data = as.dendrogram(hclust(d = dist(x = t(as.data.frame(get(i)[,-1]))))), labels = FALSE)
  dendro_ordered<-order.dendrogram(as.dendrogram(hclust(d = dist(x = t(as.data.frame(get(i)[,-1]))))))
  dendro_ordered<-c(0, dendro_ordered)
  dendro_ordered<-dendro_ordered+1
  ordered_i<-get(i)[,dendro_ordered]
  molten_i<-melt(ordered_i)
  
  molten_tibble<-tibble(molten_i) %>%
    mutate(variable=factor(variable,levels=(sort(unique(variable))))) %>%
    mutate(value_2=cut(value,breaks=c(-1, 0.01, 0.02, 0.03, 0.04, 0.05,
                                      0.06, 0.07, 0.08, 0.09, 0.10, 0.20, max(value,na.rm=T)),
                       labels=c("0-0.01","0.01-0.02","0.02-0.03","0.03-0.04",
                                "0.04-0.05","0.05-0.06","0.06-0.07","0.07-0.08","0.08-0.09",
                                "0.09-0.10", "0.10-0.20", ">0.20"))) %>%
    mutate(value_2=factor(as.character(value_2),levels=rev(levels(value_2))))
  colnames(molten_tibble)<-c("genename", "variable", "value", "CV")
  molten_tibble<-data.frame(factor(molten_tibble$genename, levels = unique(molten_tibble$genename)), molten_tibble[,-1])
  
  molten_tibble[,2]<-factor(x = molten_tibble$variable,
                            levels = colnames(ordered_i), 
                            ordered = TRUE)
  colnames(molten_tibble)<-c("genename", "variable", "value", "CV")
  assign(paste0(i,"_test_plot"), ggplot(molten_tibble, aes(variable, genename, fill= CV)) + 
           geom_tile() +
           guides(colour = guide_colourbar(barheight = unit(30, "cm"))) +
           theme(axis.text.x = element_text(size = 25, angle = 270), 
                 legend.key.height = unit(2.5, "cm"),
                 axis.text.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank(),
                 text = element_text(size = 19),
                 axis.title.x = element_text(size = 20),
                 legend.text = element_text(size = 19),
                 legend.title = element_text(size = 19)) +
           xlab("Measurement") +
           scale_fill_manual(values=rev(viridis(20))[c(1:8,14,16,18,20)]))
  
  ## Here, the plot is "wide", meaning it is sized for all timepoints, and not only drought or highlight
  png(filename = paste0(i, "_CV_wide.png"), width = 2000, height = 1080, units = "px")
  grid.newpage()
  print(dendro_plot, vp = viewport(x = 0.465, y = 0.9, width = 0.99, height = 0.15))
  print(get(paste0(i,"_test_plot")), vp = viewport(x = 0.5, y = 0.45, width = 1.00, height = 0.87))
  dev.off()
}


## Here make a plot comparing CV
## First: Unfiltered
## Extract values into long format
piv1<-(pivot_longer(ttible2_totaldro, everything()[-1]))
piv2<-(pivot_longer(ttible2_totalhili, everything()[-1]))
piv1[,2]<-"Drought"
piv2[,2]<-"Hili"
## Join them into single variable
piv_joined<-rbind(piv1, piv2)
## Calculate p-value, drought mean, and difference between means
piv_wilcox<-wilcox.test(unlist(piv1[,3]), unlist(piv2[,3]))[3]
piv_mean_diff<-signif(mean(unlist(piv1[,3])) - mean(unlist(piv2[,3])))
piv_mean_dro<-signif(mean(unlist(piv1[,3])))
## Make plot
plot_piv_1<-ggplot(piv_joined, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100)+
  ggtitle('Comparison between CV of all genes for all time-points') +
  geom_vline(xintercept = 0.04, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", piv_wilcox), x=0.22,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.04", x=0.22,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Difference between means equals ", piv_mean_diff), x=0.22,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Drought mean equals ", piv_mean_dro), x=0.22,  y=0.80, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Second: Only drought/only hili filtered
## Extract values into long format
piv1<-(pivot_longer(ttible2_0.04_Dro_alone, everything()[-1]))
piv2<-(pivot_longer(ttible2_0.04_Hili_alone, everything()[-1]))
piv1[,2]<-"Drought"
piv2[,2]<-"Hili"
## Join them into single variable
piv_joined<-rbind(piv1, piv2)
## Calculate p-value, drought mean, and difference between means
piv_wilcox<-wilcox.test(unlist(piv1[,3]), unlist(piv2[,3]))[3]
piv_mean_diff<-signif(mean(unlist(piv1[,3])) - mean(unlist(piv2[,3])))
piv_mean_dro<-signif(mean(unlist(piv1[,3])))
## Make plot
plot_piv_2<-ggplot(piv_joined, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  ggtitle('Comparison between CV of filtered genes (only HLM or only DM)') +
  geom_vline(xintercept = 0.04, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", piv_wilcox), x=0.22,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.04", x=0.22,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Difference between means equals ", piv_mean_diff), x=0.22,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Drought mean equals ", piv_mean_dro), x=0.22,  y=0.80, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Third: Both filtered
## Extract values into long format
piv1<-(pivot_longer(ttible2_0.04_Dro, everything()[-1]))
piv2<-(pivot_longer(ttible2_0.04_Hili, everything()[-1]))
piv1[,2]<-"Drought"
piv2[,2]<-"Hili"
## Join them into single variable
piv_joined<-rbind(piv1, piv2)
## Calculate p-value, drought mean, and difference between means
piv_wilcox<-wilcox.test(unlist(piv1[,3]), unlist(piv2[,3]))[3]
piv_mean_diff<-signif(mean(unlist(piv1[,3])) - mean(unlist(piv2[,3])))
piv_mean_dro<-signif(mean(unlist(piv1[,3])))
## Make plot
plot_piv_3<-ggplot(piv_joined, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100)+
  ggtitle('Comparison between CV of filtered genes (both HLM and DM)') +
  geom_vline(xintercept = 0.04, "dashed") +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", piv_wilcox), x=0.22,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.04", x=0.22,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Difference between means equals ", piv_mean_diff), x=0.22,  y=0.85, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob(paste0("Drought mean equals ", piv_mean_dro), x=0.22,  y=0.80, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))

## Bind the 3 plots together, and export
full_piv_lot <- (plot_piv_1 | plot_piv_2) / plot_piv_3 + plot_layout(guides = 'collect') 
ggsave(full_piv_lot, filename = "full_pivot_plot.png", units = "mm", width = 150, height = 100, scale = 2.8)


setwd("../output")

## Make a noise to indicate finished part
beep()

#### 5. Gene Ontology
## In this part, GO is applied.
## As two tools are not built into R, manual action is required.
## 
## 

## Read the TAIR gene-GO information
## Note: Code assumes previous part was run immediately before
setwd("../input")
GA<-read.csv("gene_association.csv.gz", header = FALSE)
setwd("../output")

## Select Biological Process annotation
GA_DF<-as.data.frame(tibble(GA) %>% filter(V9 == "P") %>% dplyr::select(V10, V5))
colnames(GA_DF)<-c("gene", "go_id")

## Perform enrichment on drought, using TAIR data
enrichment_drought<-go_enrich(atha_CV[,c(1,29)], test='wilcoxon', n_randsets=100, annotations = GA_DF)
## Perform enrichment on high light
enrichment_high.light<-go_enrich(atha_CV[,c(1,30)], test='wilcoxon', n_randsets=100, annotations = GA_DF)

## DEBUG: This calculates the number of analysed genes without biological process annotation in TAIR files
## 
# rejects<-c()
# for(i in 1:length(atha_CV[,1])){
#   if(length(grep(atha_CV[i,1], GA_DF[,1]))==0){
#     rejects<-c(rejects, i)
#   }
# }

## Saving results of enrichment separately
enrichment_tibble_dro<-tibble(enrichment_drought$results)
enrichment_tibble_hili<-tibble(enrichment_high.light$results)


## Data output:
## Write 4 files, to be entered into Gene Ontology
## Both Hili and Dro pass threshold
write(unlist(as.data.frame((ttibble2_CV %>% filter(HLM > 0.04) %>% filter(DM > 0.04))[,1])), "input_GO.txt")
## Background
write(unlist(as.data.frame((ttibble2_CV)[,1])), "background_GO.txt")
## Just drought passing threshold
write(unlist(as.data.frame((ttibble2_CV %>% filter(DM > 0.04))[,1])), "just_dro_GO.txt")
## Just High Light passing threshold
write(unlist(as.data.frame((ttibble2_CV %>% filter(HLM > 0.04))[,1])), "just_hili_GO.txt")

## INSTRUCTIONS:
## Run 6 analysis methods, all Biological Process
## For all 6, select Araport/TAIR IDs, and run with default settings besides.
## Download output, and coerce it into .csv (trim unnecessary info from Panther result text file)
##
## DAVID with input_GO and background_GO background - DAVID_background.csv
## Panther with input_GO and background_GO background - Panther_background.csv
## DAVID with just_dro_GO and background_GO background - Drought_DAVID_background.csv
## DAVID with just_hili_GO and background_GO background - Hili_DAVID_background.csv
##


## Read the files
setwd("../input")
DAVID_w_background<-read.csv("DAVID_background.csv")
PANTHER_w_background<-read.csv("Panther_background.csv")
setwd("../output")

## Select only relevant q-values
## Use substr and str_sub to produce GO term IDs following consistent format
DAVID_IDs_b_F<-substr(as.data.frame(tibble(DAVID_w_background) %>% filter(FDR < 0.05))[,2], 1, 10)
PANTHER_IDs_b_F<-str_sub(as.data.frame(tibble(PANTHER_w_background) %>% filter(input_GO.txt..FDR. < 0.05))[,1], -11,-2)

## Prepare data for merging into long format.
## Format: GO term - tool used - whether it is enriched or depleted
## DAVID only does enrichment...
DAVID_IDS_b_sign_F<-cbind(DAVID_IDs_b_F, "DAVID w/ background", "+")
## Whereas Panther does both enrichment and depletion
PANTHER_IDS_b_sign_F<-cbind(PANTHER_IDs_b_F, "PANTHER w/ background", as.data.frame(tibble(PANTHER_w_background) %>% filter(input_GO.txt..FDR. < 0.05))[,5])
colnames(DAVID_IDS_b_sign_F)<-c("GO terms", "Tool", "Sign")
colnames(PANTHER_IDS_b_sign_F)<-c("GO terms", "Tool", "Sign")


## Work with GOfuncR output. Depleted Dro
negative_enrichment_dro<-unlist(as.data.frame(enrichment_tibble_dro[,4]))
names(negative_enrichment_dro)<-c()
## Enriched Dro
positive_enrichment_dro<-unlist(as.data.frame(enrichment_tibble_dro[,5]))
names(positive_enrichment_dro)<-c()

## Depleted Hili
negative_enrichment_hili<-unlist(as.data.frame(enrichment_tibble_hili[,4]))
names(negative_enrichment_hili)<-c()
## Enriched Hili
positive_enrichment_hili<-unlist(as.data.frame(enrichment_tibble_hili[,5]))
names(positive_enrichment_hili)<-c()

## GOfuncR does not output FDR-adjusted p value. Therefore, calculate it using p.adjust for all the results below.
negative_enrichment_F_dro<-p.adjust(negative_enrichment_dro, "fdr")
positive_enrichment_F_dro<-p.adjust(positive_enrichment_dro, "fdr")
negative_enrichment_F_hili<-p.adjust(negative_enrichment_hili, "fdr")
positive_enrichment_F_hili<-p.adjust(positive_enrichment_hili, "fdr")

## Select only relevant q-values (q<0.05)
negative_IDs_F_dro<-enrichment_tibble_dro[negative_enrichment_F_dro<0.05,2]
positive_IDs_F_dro<-enrichment_tibble_dro[positive_enrichment_F_dro<0.05,2]
negative_IDs_F_hili<-enrichment_tibble_hili[negative_enrichment_F_hili<0.05,2]
positive_IDs_F_hili<-enrichment_tibble_hili[positive_enrichment_F_hili<0.05,2]

## Prepare them to match the format
negative_IDs_F_dro<-cbind(negative_IDs_F_dro, "Wilcoxon DM", "-")
positive_IDs_F_dro<-cbind(positive_IDs_F_dro, "Wilcoxon DM", "+")
negative_IDs_F_hili<-cbind(negative_IDs_F_hili, "Wilcoxon HLM", "-")
positive_IDs_F_hili<-cbind(positive_IDs_F_hili, "Wilcoxon HLM", "+")

## Coerce colnames to match
colnames(negative_IDs_F_dro)<-c("GO terms", "Tool", "Sign")
colnames(positive_IDs_F_dro)<-c("GO terms", "Tool", "Sign")
colnames(negative_IDs_F_hili)<-c("GO terms", "Tool", "Sign")
colnames(positive_IDs_F_hili)<-c("GO terms", "Tool", "Sign")

## Merge data together
background_data<-rbind(DAVID_IDS_b_sign_F, PANTHER_IDS_b_sign_F, positive_IDs_F_dro, negative_IDs_F_dro, positive_IDs_F_hili, negative_IDs_F_hili)

## Change + and - to read "Enriched" and "Depleted"
background_data[background_data[,3]=="+", 3] <- "Enriched"
background_data[background_data[,3]=="-", 3] <- "Depleted"

## Separate data into objects, based on whether or not they appear in selected datasets
## Merged object has results sorted by GO term, rather than tool
loop1<-c()
loop3<-c()
loop5<-c()
loop6<-c()
GO_table<-sort(table(background_data[,1]), decreasing = TRUE)
for(i in 1:length(GO_table)){
  IDs_loop<-background_data[background_data[,1]==names(GO_table[i]),]
  
  if(sum(IDs_loop[,2]=="DAVID w/ background")!=0){
    loop1<-rbind(loop1, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="PANTHER w/ background")!=0){
    loop3<-rbind(loop3, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="Wilcoxon DM")!=0){
    loop5<-rbind(loop5, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="Wilcoxon HLM")!=0){
    loop6<-rbind(loop6, IDs_loop)
  }
}

## Merge loops for plot:
background_data_plot<-rbind(loop1, loop3, loop5, loop6)
background_data_plot[,1]<-factor(background_data_plot[,1], rev(unique(background_data_plot[,1])))

## Generate plot
GO_w_background_plot_1<-ggplot(background_data_plot, aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 10, angle = 270), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_blank()) +
  scale_fill_manual(values=rev(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")

## Save the plot
ggsave(GO_w_background_plot_1, filename = "GO_w_background_plot_1.png", width = 7, height = 10, device = "png")


## Generate object with terms that appear in DAVID and Panther output, and either GOfuncR Drought or High light output
termplot<-c()
GO_table<-sort(table(background_data[,1]), decreasing = TRUE)
for(i in 1:length(GO_table)){
  IDs_loop<-background_data[background_data[,1]==names(GO_table[i]),]
  
  if(sum(IDs_loop[,2]=="DAVID w/ background")!=0){
    if(sum(IDs_loop[,2]=="PANTHER w/ background")!=0){
      if(sum(IDs_loop[,2]=="Wilcoxon DM" | IDs_loop[,2]=="Wilcoxon HLM")!=0){
        termplot<-rbind(termplot, IDs_loop)
      }
    }
  }
}


## GO-to-term translate plug
## As before
goterms <- Term(GOTERM)
goterm_names<-names(goterms)

term_GO<-c()
for(i in termplot[,1]){
  term_GO<-c(term_GO, goterms[i==goterm_names])
}
##

## Update term object with GO term names
termplot[,1]<-term_GO

## Generate the plot
GO_termplot<-ggplot(termplot, aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 12, angle = 270), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values=(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")
## Save the plot
ggsave(GO_termplot, filename = "GO_termplot.png", width = 7, height = 10, device = "png")


## 
setwd("../input")

## As before - read data for each tool, and process it.
Drought_DAVID_c<-read.csv("Drought_DAVID_background.csv")
Hili_DAVID_c<-read.csv("Hili_DAVID_background.csv")
setwd("../output")

Drought_DAVID<-substr(as.data.frame(tibble(Drought_DAVID_c) %>% filter(FDR < 0.05))[,2], 1, 10)
Hili_DAVID<-substr(as.data.frame(tibble(Hili_DAVID_c) %>% filter(FDR < 0.05))[,2], 1, 10)

Drought_DAVID<-cbind(Drought_DAVID, "DM DAVID", "+")
Hili_DAVID<-cbind(Hili_DAVID, "HLM DAVID", "+")

colnames(Drought_DAVID)<-c("GO terms", "Tool", "Sign")
colnames(Hili_DAVID)<-c("GO terms", "Tool", "Sign")
setwd("../output")

hili_dro_comp2<-rbind(Drought_DAVID, Hili_DAVID)
hili_dro_comp2<-as.data.frame(hili_dro_comp2)


## GO-to-term translate plug
## As before
## If terms are absent from annotation, mark them as obsolete
library(GO.db)
goterms <- Term(GOTERM)
goterm_names<-names(goterms)

term_GO<-c()
obsolete_terms_1<-c()
for(i in hili_dro_comp2[,1]){
  if(length(goterms[i==goterm_names])!=0){
    term_GO<-c(term_GO, goterms[i==goterm_names])
  }
  else{
    term_GO<-c(term_GO, "OBSOLETE TERM")
    obsolete_terms_1<-c(obsolete_terms_1, i)
  }
}
##

## Apply the GO terms
hili_dro_comp2[,1]<-term_GO
## Only select terms which *do not* appear in both dro and hili
## Also erase obsolete terms
hili_dro_names<-names(table(term_GO)[(table(term_GO))<2])

nomatch<-c()
for(i in hili_dro_names){
  nomatch<-rbind(nomatch, hili_dro_comp2[i==hili_dro_comp2[,1],])
}

## As before, replace "+" and "-" with "enriched" and "depleted"
nomatch[nomatch[,3]=="+", 3] <- "Enriched"
nomatch[nomatch[,3]=="-", 3] <- "Depleted"


### Unused in the thesis itself
# ## Generate plot of differences between high-light and drought
# GO_nomatch_plot<-ggplot(nomatch, aes(Tool, `GO terms`, fill = Sign)) + 
#   geom_tile(color = "gray") +
#   theme(axis.text.x = element_text(size = 12, angle = 270), 
#         legend.key.height = unit(2.5, "cm"),
#         axis.title.x = element_text(size = 12), 
#         axis.text.y = element_text(size = 12)) +
#   scale_fill_manual(values=(c("#009E73", "#D55E00"))) +
#   labs(fill = "Value")
# ## Save the plot
# ggsave(GO_nomatch_plot, filename = "GO_nomatch_plot.png", width = 6, height = 10, device = "png")
###

## Generate the plot of *relevant* differences between high-light and drought
relevant<-ggplot(nomatch[c(2,3,5,10, 26),], aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 12, angle = 270), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values=(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")

## Save the plot
ggsave(relevant, filename = "GO_relevant.png", width = 6, height = 3, device = "png")


##
## Reverse GO:
## Analyse downregulated terms
##
## Data output:
## Write three file, to be entered into Gene Ontology
## Use the background file generated for "regular" GO - it's the same.
setwd("../output")
write(unlist(as.data.frame((ttibble2_CV %>% filter(HLM < 0.018) %>% filter(DM < 0.018))[,1])), "input_GO_reverse.txt")

## Reverse GO is not conducted on a per-experiment basis by default
# write(unlist(as.data.frame((ttibble2_CV %>% filter(DM < 0.018))[,1])), "just_dro_GO_reverse.txt")
# write(unlist(as.data.frame((ttibble2_CV %>% filter(HLM < 0.018))[,1])), "just_hili_GO_reverse.txt")

## Make separate directory for reverse GO input
setwd("../input")
dir.create("GO_reverse")


## Generate a plot informing the choice of cutoff value
## Histogram of all all_drought and all_highlight CV values
CV_wilcox<-signif(unlist(wilcox.test(unlist(ttibble2_CV[,29]), unlist(ttibble2_CV[,30]))[3]))
pivoted_CV<-pivot_longer(ttibble2_CV[,c(1,29, 30)], c(DM, HLM))
plot_CV_REV<-ggplot(pivoted_CV, aes(x=value, color = name, fill = name)) + 
  geom_histogram(position="identity", alpha=0.5, bins = 100) +
  geom_vline(xintercept = 0.018, "dashed")+
  ggtitle('Coefficient of Variation') +
  annotation_custom(grobTree(textGrob(paste0("P-value equals ", CV_wilcox), x=0.15,  y=0.95, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  annotation_custom(grobTree(textGrob("h = 0.018", x=0.15,  y=0.9, hjust=0, gp=gpar(col="black", fontsize=13)))) +
  theme(legend.text=element_text(size=14)) +
  guides(fill=guide_legend(title="sample"), colour=guide_legend(title="sample"))
  



## Import DAVID and Panther results
setwd("GO_reverse")
DAVID_w_background_rev<-read.csv("DAVID_background.csv")
PANTHER_w_background_rev<-read.csv("Panther_background.csv")

## Select only relevant IDs by q-values
DAVID_IDs_b_F_rev<-substr(as.data.frame(tibble(DAVID_w_background_rev) %>% filter(FDR < 0.05))[,2], 1, 10)
PANTHER_IDs_b_F_rev<-str_sub(as.data.frame(tibble(PANTHER_w_background_rev) %>% filter(input_GO_reverse.txt..FDR. < 0.05))[,1], -11,-2)

## Prepare data for merging into long format.
## Format: GO term - tool used - whether it is enriched or depleted
## DAVID only does enrichment...
DAVID_IDS_b_sign_F_rev<-cbind(DAVID_IDs_b_F_rev, "DAVID w/ background", "+")
## Whereas Panther does both enrichment and depletion
PANTHER_IDS_b_sign_F_rev<-cbind(PANTHER_IDs_b_F_rev, "PANTHER w/ background", as.data.frame(tibble(PANTHER_w_background_rev) %>% filter(input_GO_reverse.txt..FDR. < 0.05))[,5])


## Use GOfuncR data from before - just need to reverse it
## i.e. - negatives are positives, and positives are negatives.
negative_IDs_F_dro_rev<-enrichment_tibble_dro[negative_enrichment_F_dro<0.05,2]
positive_IDs_F_dro_rev<-enrichment_tibble_dro[positive_enrichment_F_dro<0.05,2]
negative_IDs_F_hili_rev<-enrichment_tibble_hili[negative_enrichment_F_hili<0.05,2]
positive_IDs_F_hili_rev<-enrichment_tibble_hili[positive_enrichment_F_hili<0.05,2]

## Prepare them to match the format
negative_IDs_F_dro_rev<-cbind(negative_IDs_F_dro_rev, "Wilcoxon DM", "+")
positive_IDs_F_dro_rev<-cbind(positive_IDs_F_dro_rev, "Wilcoxon DM", "-")
negative_IDs_F_hili_rev<-cbind(negative_IDs_F_hili_rev, "Wilcoxon HLM", "+")
positive_IDs_F_hili_rev<-cbind(positive_IDs_F_hili_rev, "Wilcoxon HLM", "-")

reverse_colnames<-c("GO terms", "Tool", "Sign")
colnames(DAVID_IDS_b_sign_F_rev)<-reverse_colnames
colnames(PANTHER_IDS_b_sign_F_rev)<-reverse_colnames
colnames(negative_IDs_F_dro_rev)<-reverse_colnames
colnames(positive_IDs_F_dro_rev)<-reverse_colnames
colnames(negative_IDs_F_hili_rev)<-reverse_colnames
colnames(positive_IDs_F_hili_rev)<-reverse_colnames

background_data_reversed<-rbind(DAVID_IDS_b_sign_F_rev, PANTHER_IDS_b_sign_F_rev, negative_IDs_F_dro_rev, positive_IDs_F_dro_rev, negative_IDs_F_hili_rev, positive_IDs_F_hili_rev)

## As before, change + and - to read "Enriched" and "Depleted"
background_data_reversed[background_data_reversed[,3]=="+", 3] <- "Enriched"
background_data_reversed[background_data_reversed[,3]=="-", 3] <- "Depleted"
###

## Sorting loop, as before
loop1<-c()
loop3<-c()
loop4<-c()
loop5<-c()
loop6<-c()
GO_table<-sort(table(background_data[,1]), decreasing = TRUE)
for(i in 1:length(GO_table)){
  IDs_loop<-background_data_reversed[background_data_reversed[,1]==names(GO_table[i]),]
  
  if(sum(IDs_loop[,2]=="DAVID w/ background")!=0){
    loop1<-rbind(loop1, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="PANTHER w/ background")!=0){
    loop3<-rbind(loop3, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="Wilcoxon DM")!=0){
    loop5<-rbind(loop5, IDs_loop)
  }
  else if(sum(IDs_loop[,2]=="Wilcoxon HLM")!=0){
    loop6<-rbind(loop6, IDs_loop)
  }
}

## Merge loops for plot:
background_data_plot_rev<-rbind(loop1, loop3, loop5, loop6)
background_data_plot_rev[,1]<-factor(background_data_plot_rev[,1], rev(unique(background_data_plot_rev[,1])))

## Generate plot
GO_w_background_plot_1_rev<-ggplot(background_data_plot_rev, aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 10, angle = 270), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 10),
        axis.text.y = element_blank()) +
  scale_fill_manual(values=rev(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")

setwd("../../output")
dir.create("GO_reverse")
setwd("GO_reverse")

## Save the plot
ggsave(GO_w_background_plot_1_rev, filename = "REVERSE_GO_w_background_plot_1.png", width = 7, height = 10, device = "png")

## Save the cutoff value plot for reverse GO
ggsave(plot_CV_REV, filename = "plot_CV_REV.png", units = "mm", width = 150, height = 60, scale = 2.75)


## Generate data for plot where a term must appear in output of at least 3 tools
names_background_data_plot_3_rev<-names(table(background_data_plot_rev[,1])[table(background_data_plot_rev[,1]) > 2])
rep_background_3_rev<-c()
for(i in 1:length(names_background_data_plot_3_rev)){
  loop3<-background_data_plot_rev[background_data_plot_rev[,1]==names_background_data_plot_3_rev[i],]
  rep_background_3_rev<-rbind(rep_background_3_rev, loop3)
}

## generate plot:
match3_plot_rev<-ggplot(rep_background_3_rev, aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 10, angle = 0), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 10)) +
  scale_fill_manual(values=rev(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")

## GO-to-term translate plug
## Translates GO-IDs to GO terms
goterms <- Term(GOTERM)
goterm_names<-names(goterms)

term_GO<-c()
for(i in rep_background_3_rev[,1]){
  term_GO<-c(term_GO, goterms[i==goterm_names])
}
##

## Applies GO term names, and coerces the result to follow a selected order (needed for plot)  
rep_background_3_rev[,1]<-term_GO
rep_background_3_rev[,1]<-factor(rep_background_3_rev[,1], levels = (unique(rep_background_3_rev[,1])))

## Generate the plot
GO_w_background_plot_rev<-ggplot(rep_background_3_rev, aes(Tool, `GO terms`, fill = Sign)) + 
  geom_tile(color = "gray") +
  theme(axis.text.x = element_text(size = 12, angle = 270), 
        legend.key.height = unit(2.5, "cm"),
        axis.title.x = element_text(size = 12), 
        axis.text.y = element_text(size = 12)) +
  scale_fill_manual(values=rev(c("#009E73", "#D55E00"))) +
  labs(fill = "Value")

## Save the plot
ggsave(GO_w_background_plot_rev, filename = "GO_w_background_plot_rev.png", width = 7, height = 10, device = "png")

setwd("../")


## Saves atha_CV and atha_means, for epigenetics analysis
save(atha_CV, file = "atha_CV.R")
save(atha_means, file = "atha_means.R")


## Plot to table converter, for appendix
for(z in c("nomatch", "background_data_plot", "background_data_plot_rev")){
  
  ## Select the plot
  sel_plot<-get(z)
  columns<-unique(sel_plot[,2])
  terms<-unique(sel_plot[,1])
  
  termloop<-c()
  colloop<-c()
  
  ## For each term, for each column (measurement), see if it is enriched or depleted.
  ## If absent, leave empty space
  for(o in terms){
    for(i in columns){
      if(any(sel_plot[sel_plot[,2]==i,1]==o)){
        colloop<-c(colloop, sel_plot[sel_plot[,2]==i,3][sel_plot[sel_plot[,2]==i,1]==o])
      }
      else{
        colloop<-c(colloop, "")
      }
    }
    termloop<-rbind(termloop, colloop)
    colloop<-c()
  }
  
  ## Correct column names
  colnames(termloop)<-columns
  
  
  ## GO-to-term translate plug
  goterms <- Term(GOTERM)
  goterm_names<-names(goterms)
  
  ## Translate GO terms
  ## If GO terms are absent from annotation because of GO version mismatch between TAIR annotation and GO.db, make a note.
  term_GO<-c()
  for(i in terms){
    if(length(goterms[i==goterm_names])!=0){
      term_GO<-c(term_GO, paste0(goterms[i==goterm_names], " - ", i))
    }
    else{
      term_GO<-c(term_GO, paste0("!!! absent from annotation - ", i))
    }
  }
  ##
  ## Apply updated term-names
  rownames(termloop)<-term_GO
  
  ## Save the files. Edit names manually.
  write.csv(termloop, paste0("Supplementary_", z, "_table.csv"))
  
}













