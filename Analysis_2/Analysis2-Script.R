##
## Run file - Analysis 2
## Requires :
## - DLM and HLM from analysis 1
## - Atha and Prom Rdata files from command-line tools, and separate DMRcaller scripts
## - CX report files from command-line tools
## Requires 1 folder underneath wd: results, split into epi1, epi2, and epi3
##

## Load the libraries
library(viridis)
library(ggplot2)
library(DMRcaller)
library(patchwork)
library(tidyverse)
library(ggExtra)
library(ggpmisc)
library(scales)
library(grid)
library(patchwork)
library(ggforce)

## Step 1: Preparation.
## Load the R files necessary
#####

## Import the data
## Import data for the first bio-replicate
setwd("results/epi1/")
load("CG-Atha.RData")
load("CHG-Atha.RData")
load("CHH-Atha.RData")
load("CG-Prom.RData")
load("CHG-Prom.RData")
load("CHH-Prom.RData")

## The variables are re-used for importing, so assign new ones
out_grange_CG_A_1<-out_grange_CG_A
out_grange_CHG_A_1<-out_grange_CHG_A
out_grange_CHH_A_1<-out_grange_CHH_A

out_grange_CG_P_1<-out_grange_CG_P
out_grange_CHG_P_1<-out_grange_CHG_P
out_grange_CHH_P_1<-out_grange_CHH_P

## Remove old variabiles, just in case.
rm(out_grange_CG_A, out_grange_CHG_A, out_grange_CHH_A, out_grange_CG_P, out_grange_CHG_P, out_grange_CHH_P)

## Repeat for Biorep 2
setwd("../epi2")
load("CG-Atha.RData")
load("CHG-Atha.RData")
load("CHH-Atha.RData")
load("CG-Prom.RData")
load("CHG-Prom.RData")
load("CHH-Prom.RData")

out_grange_CG_A_2<-out_grange_CG_A
out_grange_CHG_A_2<-out_grange_CHG_A
out_grange_CHH_A_2<-out_grange_CHH_A

out_grange_CG_P_2<-out_grange_CG_P
out_grange_CHG_P_2<-out_grange_CHG_P
out_grange_CHH_P_2<-out_grange_CHH_P

rm(out_grange_CG_A, out_grange_CHG_A, out_grange_CHH_A, out_grange_CG_P, out_grange_CHG_P, out_grange_CHH_P)

## Repeat for Biorep 3
setwd("../epi3")
load("CG-Atha.RData")
load("CHG-Atha.RData")
load("CHH-Atha.RData")
load("CG-Prom.RData")
load("CHG-Prom.RData")
load("CHH-Prom.RData")

out_grange_CG_A_3<-out_grange_CG_A
out_grange_CHG_A_3<-out_grange_CHG_A
out_grange_CHH_A_3<-out_grange_CHH_A

out_grange_CG_P_3<-out_grange_CG_P
out_grange_CHG_P_3<-out_grange_CHG_P
out_grange_CHH_P_3<-out_grange_CHH_P

rm(out_grange_CG_A, out_grange_CHG_A, out_grange_CHH_A, out_grange_CG_P, out_grange_CHG_P, out_grange_CHH_P)

## Reset the WD
setwd("../..")

## Import: Coefficient of Variability, Mean
## As calculated before.
load("atha_CV.R")
load("atha_means.R")
ttibble2_CV<-tibble(atha_CV)
ttibble2_means<-tibble(atha_means)

colnames(ttibble2_CV)<-c(colnames(ttibble2_CV)[-c(29:30)], "DM", "HLM")
colnames(ttibble2_means)<-c(colnames(ttibble2_means)[-c(29:30)], "DM", "HLM")

#####
## Step 2: Calculations
## Process the data, and generate files for later plotting
#####

## Extract Proportions (so, proportion CG, CHG, CHH)
all_bound<-c()
for(i in c("out_grange_CG_A_1", "out_grange_CG_A_2", "out_grange_CG_A_3", "out_grange_CHG_A_1", "out_grange_CHG_A_2", "out_grange_CHG_A_3", "out_grange_CHH_A_1", "out_grange_CHH_A_2", "out_grange_CHH_A_3",
           "out_grange_CG_P_1", "out_grange_CG_P_2", "out_grange_CG_P_3", "out_grange_CHG_P_1", "out_grange_CHG_P_2", "out_grange_CHG_P_3", "out_grange_CHH_P_1", "out_grange_CHH_P_2", "out_grange_CHH_P_3")){
  workrange<-get(i)
  all_bound<-cbind(all_bound, unlist(mcols(workrange)[3]))
}
all_bound<-data.frame(names(out_grange_CG_P_1), all_bound)

## Change column name scheme
colnames(all_bound)<-c("genename", "CG_gbM_1", "CG_gbM_2", "CG_gbM_3", "CHG_gbM_1", "CHG_gbM_2", "CHG_gbM_3", "CHH_gbM_1", "CHH_gbM_2", "CHH_gbM_3",
                       "CG_promoter_1", "CG_promoter_2", "CG_promoter_3", "CHG_promoter_1", "CHG_promoter_2", "CHG_promoter_3", "CHH_promoter_1", "CHH_promoter_2", "CHH_promoter_3")
rownames(all_bound)<-names(out_grange_CG_P_1)

## Remove NAs from all_bound
vector_NAs<-c()
for(i in 1:length(all_bound[1,])){
  number_NAs<-(1:19239)[(is.na(all_bound[,i]))]
  vector_NAs<-c(vector_NAs, number_NAs)
}
NA_containing_entries<-all_bound[unique(vector_NAs),]
all_bound<-all_bound[-unique(vector_NAs),]

## Remove NAs from ttibble2_CV and ttibble2_means
vector_NAs_2<-c()
for(i in rownames(NA_containing_entries)){
  number_NAs_2<-(1:19239)[unlist(ttibble2_CV[,1])==i]
  vector_NAs_2<-c(vector_NAs_2, number_NAs_2)
}

ttibble2_CV<-ttibble2_CV[-unique(vector_NAs_2),]
ttibble2_means<-ttibble2_means[-unique(vector_NAs_2),]

## Calculate mean value of each proportion
bound_means<-c()
for(i in 1:length(all_bound[,1])){
  loop_means<-c(mean(unlist(all_bound[i,2:4])), mean(unlist(all_bound[i,5:7])), mean(unlist(all_bound[i,8:10])), mean(unlist(all_bound[i,11:13])), mean(unlist(all_bound[i,14:16])), mean(unlist(all_bound[i,17:19])))
  bound_means_l<-t(loop_means)
  bound_means<-rbind(bound_means, bound_means_l)
}

bound_means<-data.frame(all_bound[,1], bound_means)
colnames(bound_means)<-c("genenames", "CG_gbM", "CHG_gbM", "CHH_gbM", "CG_promoter", "CHG_promoter", "CHH_promoter")

## Append DM and HLM to the proportion values
joined_all<-c()
for(i in unlist(ttibble2_CV[,1])){
  b_loop<-all_bound[all_bound[,1]==i,]
  joined_all<-rbind(joined_all, b_loop)
}
joined_all<-data.frame(joined_all, ttibble2_CV[,29:30], ttibble2_means[,29:30])
colnames(joined_all)<-c(colnames(joined_all)[-(20:23)], "DM_CV", "HLM_CV", "DM_mean", "HLM_mean")

## Append DM and HLM to the mean values
joined_means<-c()
for(i in unlist(ttibble2_CV[,1])){
  b_loop<-bound_means[bound_means[,1]==i,]
  joined_means<-rbind(joined_means, b_loop)
}
joined_means<-data.frame(joined_means, ttibble2_CV[,29:30], ttibble2_means[,29:30])
colnames(joined_means)<-c(colnames(joined_means)[-(8:11)], "DM_CV", "HLM_CV", "DM_mean", "HLM_mean")

joined_means_vector<-rep("gbM", (dim(joined_means)[1]))
joined_means_vector[joined_means[,3] > 0.05 | joined_means[,4] > 0.05]<-"Transposable-like"

joined_means<-cbind(joined_means, joined_means_vector)

pivoted_means_1<-pivot_longer(joined_means, everything()[-c(1,8:12)])



## Create function to assign methylation types - gbM and transposable-like, and promoter
## Because entries are all in a single table, gene methylation data (gbm and transposable-like) lies next to promoters
## But because of this, plots separate them.
type_assigner<-function(x){
  if(unlist(strsplit(unlist(x[7]), "_"))[2]=="promoter"){
    return("promoter")
  }
  else{
    return(unlist(x[6]))
  }
}

pivoted_means_1[,6]<-apply(pivoted_means_1, 1, type_assigner)

pivoted_means_1<-cbind(pivoted_means_1, c(rep(c("CG", "CHG", "CHH"), 2)))
colnames(pivoted_means_1)<-c("genename", "DM_CV", "HLM_CV", "DM_mean", "HLM_mean", "Measurement_Type", "Measurement", "Mean_Epi", "Methylation_Type")
pivoted_means_1[,6]<-factor(pivoted_means_1[,6], levels = c("gbM", "Transposable-like", "promoter"))


#####
## Step 3: Make plots!
## Process the data, and generate files for later plotting
#####

## 1. Mean vs Mean


pivoted_means_0<-pivoted_means_1 %>% mutate(state = ifelse(Mean_Epi > 0.1,"Methylated","Non-Methylated"))

## Create a function to categorise inputs for drought
assign_state_mean_mean_dro<-function(ent){
  ## Methylation thresholds:
  ## For CG, 0.1
  ## For CHG and CHH, 0.05
  if(ent[9]=="CG"){
    if(ent[8]>0.1){
      if(ent[4]<(8.4)){
        return("Methylated, Low expression")
      }
      if(ent[4]>(12.05)){
        return("Methylated, High expression")
      }
      else{
        return("Methylated, Medium expression")
      }
    }
    else{
      if(ent[4]<8.4){
        return("Non-methylated, Low expression")
      }
      if(ent[4]>12.05){
        return("Non-methylated, High expression")
      }
      else{
        return("Non-methylated, Medium expression")
      }
    }
  }
  else{
    if(ent[8]>0.05){
      if(ent[4]<8.4){
        return("Methylated, Low expression")
      }
      if(ent[4]>12.05){
        return("Methylated, High expression")
      }
      else{
        return("Methylated, Medium expression")
      }
    }
    else{
      if(ent[4]<8.4){
        return("Non-methylated, Low expression")
      }
      if(ent[4]>12.05){
        return("Non-methylated, High expression")
      }
      else{
        return("Non-methylated, Medium expression")
      }
    }
  }
}
## And a separate one to do the same for high light
assign_state_mean_mean_hili<-function(ent){
  ## For CG, 0.1
  ## For CHG and CHH, 0.5
  if(ent[9]=="CG"){
    if(ent[8]>0.1){
      if(ent[5]<8.4){
        return("Methylated, Low expression")
      }
      if(ent[5]>12.05){
        return("Methylated, High expression")
      }
      else{
        return("Methylated, Medium expression")
      }
    }
    else{
      if(ent[5]<8.4){
        return("Non-methylated, Low expression")
      }
      if(ent[5]>12.05){
        return("Non-methylated, High expression")
      }
      else{
        return("Non-methylated, Medium expression")
      }
    }
  }
  else{
    if(ent[8]>0.05){
      if(ent[5]<8.4){
        return("Methylated, Low expression")
      }
      if(ent[5]>12.05){
        return("Methylated, High expression")
      }
      else{
        return("Methylated, Medium expression")
      }
    }
    else{
      if(ent[5]<8.4){
        return("Non-methylated, Low expression")
      }
      if(ent[5]>12.05){
        return("Non-methylated, High expression")
      }
      else{
        return("Non-methylated, Medium expression")
      }
    }
  }
}

## Apply to drought
## Normally this would be done using apply, but for some reason it does not work here
## So a for-loop it is instead
assign_loop_d<-c()
for(i in 1:length(pivoted_means_1[,1])){
  assign_loop_d<-c(assign_loop_d, assign_state_mean_mean_dro(pivoted_means_1[i,]))
}
pivoted_means_1_dro<-cbind(pivoted_means_1, assign_loop_d)
colnames(pivoted_means_1_dro)<-c(colnames(pivoted_means_1_dro)[1:9], "state")

## Apply to high light
assign_loop_h<-c()
for(i in 1:length(pivoted_means_1[,1])){
  assign_loop_h<-c(assign_loop_h, assign_state_mean_mean_hili(pivoted_means_1[i,]))
}
pivoted_means_1_hili<-cbind(pivoted_means_1, assign_loop_h)
colnames(pivoted_means_1_hili)<-colnames(pivoted_means_1_dro)

## Order levels properly, so the plot is easier to read
pivoted_means_1_dro$state <- factor(pivoted_means_1_dro$state, levels = c("Methylated, Low expression", "Non-methylated, Low expression", "Methylated, Medium expression",
                                                                          "Non-methylated, Medium expression", "Methylated, High expression", "Non-methylated, High expression"))
pivoted_means_1_hili$state <- factor(pivoted_means_1_hili$state, levels = c("Methylated, Low expression", "Non-methylated, Low expression", "Methylated, Medium expression",
                                                                          "Non-methylated, Medium expression", "Methylated, High expression", "Non-methylated, High expression"))

## Create the plots. First the 3 Drought ones...
D_MVM <- ggplot(pivoted_means_1_dro, aes(DM_mean, Mean_Epi, color = state)) + 
  geom_point() +
  facet_grid(cols = vars(Measurement_Type), rows = vars(Methylation_Type)) +
  xlab("Mean of expression (DM)") +
  ylab("Mean of methylation proportions in the context") +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values=(c("#FF3636", "#721D1D","#BDBDBD", "#737373","#6FAFFF", "#2A2C62"))) +
  xlim(c(7.5,16.5)) +
  theme(strip.text.x = element_text(size = 13, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"), legend.position="top",
        legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Then the 3 high light ones...
H_MVM <- ggplot(pivoted_means_1_hili, aes(HLM_mean, Mean_Epi, color = state)) + 
  geom_point() +
  facet_grid(cols = vars(Measurement_Type), rows = vars(Methylation_Type)) +
  xlab("Mean of expression (HLM)") +
  ylab("Mean of methylation proportions in the context") +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values=(c("#FF3636", "#721D1D","#BDBDBD", "#737373","#6FAFFF", "#2A2C62"))) +
  xlim(c(7.5,16.5)) +
  theme(strip.text.x = element_text(size = 13, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"), legend.position="top",
        legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Then join them together...
p_MVM<-(D_MVM | H_MVM) + plot_layout(guides = 'collect') & theme(legend.position = 'top')



## 2. CV vs Mean
## Could be done using mutate, but forloop gives greater control

## Generate the "state" layer for hili
## Split by methylation (0.1 for CG, 0.05 for CHG and CHH) and expression (0.04 for all)
hili_state<-c()
for(i in 1:length(pivoted_means_1[,1])){
  ## For CG, 0.1
  ## For CHG and CHH, 0.5
  if(pivoted_means_1[i,9]=="CG"){
    if(pivoted_means_1[i,8]>0.1){
      if(pivoted_means_1[i,3]>0.04){
        loopvar<-"Methylated, High variability"
      }
      else{
        loopvar<-"Methylated, Low variability"
      }
    }
    else{
      if(pivoted_means_1[i,3]>0.04){
        loopvar<-"Non-methylated, High variability"
      }
      else{
        loopvar<-"Non-methylated, Low variability"
      }
    }
  }
  else{
    if(pivoted_means_1[i,8]>0.05){
      if(pivoted_means_1[i,3]>0.04){
        loopvar<-"Methylated, High variability"
      }
      else{
        loopvar<-"Methylated, Low variability"
      }
    }
    else{
      if(pivoted_means_1[i,3]>0.04){
        loopvar<-"Non-methylated, High variability"
      }
      else{
        loopvar<-"Non-methylated, Low variability"
      }
    }
  }
  hili_state<-c(hili_state, loopvar)
}

## Apply the state layer, and generate the hili part of the plot
pivoted_means_2_h<-data.frame(pivoted_means_1, hili_state)
colnames(pivoted_means_2_h)<-c(colnames(pivoted_means_1), "state")
hili_CV_vs_mean_plot <- ggplot(pivoted_means_2_h, aes(HLM_CV, Mean_Epi, color = state)) + 
  geom_point() +
  facet_grid(cols = vars(Measurement_Type), rows = vars(Methylation_Type)) +
  xlab("CV of expression (HLM)") +
  geom_vline(xintercept = 0.04, "dashed", alpha = 0.5) +
  ylab("Mean of methylation proportions in the context") +
  theme(axis.text.x = element_text(size = 13.5), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535", "#83239E"))) +
  xlim(c(0,0.18)) +
  theme(strip.text.x = element_text(size = 13, face = "bold"),
        strip.text.y = element_text(size = 13, face = "bold"), legend.position="top",
        legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Generate the "state" layer for drought
dro_state<-c()
for(i in 1:length(pivoted_means_1[,1])){
  ## For CG, 0.1
  ## For CHG and CHH, 0.5
  if(pivoted_means_1[i,9]=="CG"){
    if(pivoted_means_1[i,8]>0.1){
      if(pivoted_means_1[i,2]>0.04){
        loopvar<-"Methylated, High variability"
      }
      else{
        loopvar<-"Methylated, Low variability"
      }
    }
    else{
      if(pivoted_means_1[i,2]>0.04){
        loopvar<-"Non-methylated, High variability"
      }
      else{
        loopvar<-"Non-methylated, Low variability"
      }
    }
  }
  else{
    if(pivoted_means_1[i,8]>0.05){
      if(pivoted_means_1[i,2]>0.04){
        loopvar<-"Methylated, High variability"
      }
      else{
        loopvar<-"Methylated, Low variability"
      }
    }
    else{
      if(pivoted_means_1[i,2]>0.04){
        loopvar<-"Non-methylated, High variability"
      }
      else{
        loopvar<-"Non-methylated, Low variability"
      }
    }
  }
  dro_state<-c(dro_state, loopvar)
}

## Apply the state layer, and generate the plot
pivoted_means_2_d<-data.frame(pivoted_means_1, dro_state)
colnames(pivoted_means_2_d)<-c(colnames(pivoted_means_1), "state")
dro_CV_vs_mean_plot <- ggplot(pivoted_means_2_d, aes(DM_CV, Mean_Epi, color = state)) + 
  geom_point() +
  facet_grid(cols = vars(Measurement_Type), rows = vars(Methylation_Type)) +
  xlab("CV of expression (DM)") +
  geom_vline(xintercept = 0.04, "dashed", alpha = 0.5) +
  ylab("Mean of methylation proportions in the context") +
  theme(axis.text.x = element_text(size = 13), 
        axis.text.y = element_text(size = 15),
        text = element_text(size = 15)) +
  scale_color_manual(values=rev(c("#9C9CA8", "#3571CA", "#CA3535", "#83239E"))) +
  xlim(c(0,0.18)) +
  theme(strip.text.x = element_text(size = 13, face = "bold"),
    strip.text.y = element_text(size = 13, face = "bold"), legend.position="top",
    legend.text = element_text(size=15)) +
  guides(color = guide_legend(override.aes = list(size = 4)))

## Merge the plots
joined_CV_mean<-(dro_CV_vs_mean_plot | hili_CV_vs_mean_plot) + plot_layout(guides = 'collect') & theme(legend.position = 'top')


## Saves both plots
dir.create("plots")
setwd("plots")
ggsave(joined_CV_mean, filename = "joined_CV_mean.png", width = 10, height = 10, device = "png", scale = 1.5)
ggsave(p_MVM, filename = "joined_mean_mean.png", width = 10, height = 10, device = "png", scale = 1.5)
setwd("../")




## Boxplots
## Generate boxplots of distributions

setwd("plots")

methylation<-unique(pivoted_means_1_dro$Methylation_Type)
measurement<-c("gbM", "Transposable-like", "promoter")

## Join states from dro and hili
pivoted_means_1_mean_mean<-cbind(pivoted_means_1_dro, pivoted_means_1_hili[,10])
colnames(pivoted_means_1_mean_mean)<-c(colnames(pivoted_means_1_mean_mean)[1:9], "state_d", "state_h")

## Begin boxplot generation
png("boxplot_mean_nooutliers_2.jpg", width = 2000, height = 2000, res = 250)
par(mfrow = c(3, 3), mar = c(2.5, 2.5, 2, 2.5))
ylims<-c(1.3, 1.2, 2.5, 1.3, 1.2, 2.0, 1.2, 1.2, 1.5)
textloop<-c()
## Ylims have to be specified manually, decided by iteration, to leave space for manually-added annotations
iteration<-0

for(i in methylation){
  for(z in measurement){
    iteration<-iteration+1
    work_tibble2<-tibble(pivoted_means_1_mean_mean) %>% filter(Measurement_Type == z & Methylation_Type == i)
    ## Split data into genes that are hi, low and mid in both samples
    both_hi <- (work_tibble2 %>% filter(state_d == "Methylated, High expression" | state_d == "Non-methylated, High expression") %>% filter(state_h == "Methylated, High expression" | state_h == "Non-methylated, High expression"))[,8]
    both_low <- (work_tibble2 %>% filter(state_d == "Methylated, Low expression" | state_d == "Non-methylated, Low expression") %>% filter(state_h == "Methylated, Low expression" | state_h == "Non-methylated, Low expression"))[,8]
    both_mid <- (work_tibble2 %>% filter(state_d == "Methylated, Medium expression" | state_d == "Non-methylated, Medium expression") %>% filter(state_h == "Methylated, Medium expression" | state_h == "Non-methylated, Medium expression"))[,8]
    ## Make boxplot
    boxplot_data<-data.frame(c(unlist(both_low), rep(NA, (length(unlist(both_mid))-length(unlist(both_low))))),
                             both_mid, 
                             c(unlist(both_hi), rep(NA, (length(unlist(both_mid))-length(unlist(both_hi))))))
    colnames(boxplot_data)<-c("low", "medium", "high")
    ## Write the p-values
    textloop<-rbind(textloop,unlist(c(wilcox.test(boxplot_data[,1], boxplot_data[,2])[3], (wilcox.test(boxplot_data[,1], boxplot_data[,3])[3]), (wilcox.test(boxplot_data[,2], boxplot_data[,3])[3]))))
    boxplot(boxplot_data, main = paste0(i, " ",z), outline=FALSE, ylim = c(boxplot.stats(boxplot_data)$stats[1]*0.80, boxplot.stats(boxplot_data)$stats[5]*ylims[iteration]))
    
  }
}
colnames(textloop)<-c("low-to-mid", "low-to-high", "mid-to-high")
## Save the p-values
## order goes: left-to-right, up-to-down
write.table(textloop, "meansplot_p.txt")
dev.off()



## Set function to assign variability thresholds
assign_func_box<-function(x){
  if(x[2]<=0.04 & x[3]<=0.04){
    return("Below_threshold")
  }
  if(x[2]>0.04){
    if(x[3]>0.04){
      return("Variable_both")
    }
    else{
      return("Variable_drought")
    }
  }
  else{
    return("Variable_hili")
  }
}

## Generate the state layer
state_layer<-(apply(pivoted_means_1, 1, assign_func_box))


pivoted_means_1_w_state<-cbind(pivoted_means_1, state_layer)
colnames(pivoted_means_1_w_state)<-c(colnames(pivoted_means_1), "state")
textloop2<-c()
meandiff_CV<-c()
png("boxplot_CV_nooutliers_2.jpg", width = 2000, height = 2000, res = 250)
par(mfrow = c(3, 3), mar = c(2.5, 2.5, 2, 2.5))
for(i in methylation){ # methylation
  for(z in measurement){ #measurement
    work_tibble2<-tibble(pivoted_means_1_w_state) %>% filter(Measurement_Type == z & Methylation_Type == i)
    ## Here, instead gene must be variable or non-variable in both
    both_data<-(work_tibble2 %>% filter(state == "Variable_both"))[,8]
    neither_data<-(work_tibble2 %>% filter(state == "Below_threshold"))[,8]
    
    boxplot_data<-data.frame(neither_data, 
                             c(unlist(both_data), rep(NA, (length(unlist(neither_data))-length(unlist(both_data))))))
    colnames(boxplot_data)<-c("Neither variable", "Both variable")
    textloop2<-c(textloop2,unlist(c(wilcox.test(boxplot_data[,1], boxplot_data[,2])[3])))
    meandiff_CV<-c(meandiff_CV, abs(mean(unlist(both_data))-mean(unlist(neither_data))))
        
    boxplot(boxplot_data, main = paste0(i, " ", z), outline=FALSE, ylim = c(boxplot.stats(boxplot_data)$stats[1]*0.80, boxplot.stats(boxplot_data)$stats[5]*1.20))
  }
}
##
write(textloop2, "CVplot_p.txt")
write(meandiff_CV, "meandiff_CV.txt")
dev.off()

setwd("../")


#####
## Step 4: Gene Ontology
## Perform the GO analysis of the 4 categories in CG methylation
## Because of significantly greater amount of categories, no wilcox analysis is performed.
## Instead, this just generates files to be used with PANTHER
#####

## Isolate CG data
CG_only_CG_met<-pivoted_means_1_w_state %>% filter(Measurement_Type == "gbM") %>% filter(Measurement == "CG_gbM")

## Write state assignment function
CG_only_CG_met_assign<-function(x){
  if(x[2]>0.04 & x[3]>0.04){
    if(x[8]>0.1){
      return("Variable, Methylated")
    }
    else{
      return("Variable, Non-methylated")
    }
  }
  if(x[2]<=0.04 & x[3]<=0.04){
    if(x[8]>0.1){
      return("Non-variable, Methylated")
    }
    else{
      return("Non-variable, Non-methylated")
    }
  }
  return("HLM-DM mismatch")
}

## Assign the state
CG_only_CG_met<-cbind(CG_only_CG_met, apply(CG_only_CG_met, 1, CG_only_CG_met_assign))
colnames(CG_only_CG_met)<-c(colnames(CG_only_CG_met)[-11], "state2")

## Write the GO files for use with PANTHER
## Per-state genes + background
dir.create("GeneOntology")
setwd("GeneOntology")
write((CG_only_CG_met %>% filter(state2 == "Variable, Methylated"))[,1], "var_met.txt")
write((CG_only_CG_met %>% filter(state2 == "Variable, Non-methylated"))[,1], "var_nonmet.txt")
write((CG_only_CG_met %>% filter(state2 == "Non-variable, Methylated"))[,1], "nonvar_met.txt")
write((CG_only_CG_met %>% filter(state2 == "Non-variable, Non-methylated"))[,1], "nonvar_nonmet.txt")
write(CG_only_CG_met[,1], "background.txt")
setwd("../")

#####
## Step 5: Base analysis.
## Generate coverage plots, profiles, and correlation
#####

## Load data
setwd("epi1")
epi1<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../epi2")
epi2<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../epi3")
epi3<-readBismark("aligned.deduplicated.CX_report.txt.CX_report.txt")
setwd("../")

## Correct data for chloroplast methylation

chloroplast_data<-epi1[seqnames(epi1)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
epi1_corrected<-epi1
epi1_corrected$readsM<-round((epi1$readsM)-(epi1$readsN)*(1-rate))
epi1_corrected$readsM[epi1_corrected$readsM<0]<-0
epi1_corrected$readsN<-round(epi1$readsN*rate)

chloroplast_data<-epi2[seqnames(epi2)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
epi2_corrected<-epi2
epi2_corrected$readsM<-round((epi2$readsM)-(epi2$readsN)*(1-rate))
epi2_corrected$readsM[epi2_corrected$readsM<0]<-0
epi2_corrected$readsN<-round(epi2$readsN*rate)


chloroplast_data<-epi3[seqnames(epi3)=="chloroplast"]
rate<-(1-(sum(mcols(chloroplast_data)$readsM)/sum(mcols(chloroplast_data)$readsN)))
epi3_corrected<-epi3
epi3_corrected$readsM<-round((epi3$readsM)-(epi3$readsN)*(1-rate))
epi3_corrected$readsM[epi3_corrected$readsM<0]<-0
epi3_corrected$readsN<-round(epi3$readsN*rate)

## Remove the space-consuming cx report files
rm(epi1, epi2, epi3)

## Make chromosome 1 object
regions_loop<-GRanges(seqnames = Rle(1), ranges = IRanges(1,30427670))

###########################################################
#####
#### CG

## Extract data from 3 bio-replicates. Low resolution profile, coverage, and spatial correlation 
profile_1_corr_CG <- computeMethylationProfile(epi1_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_1_corr_CG <-computeMethylationDataCoverage(epi1_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_1_corr_CG <-computeMethylationDataSpatialCorrelation(epi1_corrected,
                                                             context="CG",
                                                             distances=c(1,10,100,1000,10000))

profile_2_corr_CG <- computeMethylationProfile(epi2_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_2_corr_CG <-computeMethylationDataCoverage(epi2_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_2_corr_CG <-computeMethylationDataSpatialCorrelation(epi2_corrected,
                                                             context="CG",
                                                             distances=c(1,10,100,1000,10000))


profile_3_corr_CG <- computeMethylationProfile(epi3_corrected,
                                               regions_loop,
                                               windowSize = 500000,
                                               context = "CG")
coverage_3_corr_CG <-computeMethylationDataCoverage(epi3_corrected, context = "CG", breaks = c(1,5,10,15, 20, 25))
spacorr_3_corr_CG <-computeMethylationDataSpatialCorrelation(epi3_corrected,
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

profile_1_corr_CHG <- computeMethylationProfile(epi1_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_1_corr_CHG <-computeMethylationDataSpatialCorrelation(epi1_corrected,
                                                              context="CHG",
                                                              distances=c(1,10,100,1000,10000))

profile_2_corr_CHG <- computeMethylationProfile(epi2_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_2_corr_CHG <-computeMethylationDataSpatialCorrelation(epi2_corrected,
                                                              context="CHG",
                                                              distances=c(1,10,100,1000,10000))


profile_3_corr_CHG <- computeMethylationProfile(epi3_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHG")
spacorr_3_corr_CHG <-computeMethylationDataSpatialCorrelation(epi3_corrected,
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
profile_1_corr_CHH <- computeMethylationProfile(epi1_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_1_corr_CHH <-computeMethylationDataSpatialCorrelation(epi1_corrected,
                                                              context="CHH",
                                                              distances=c(1,10,100,1000,10000))

profile_2_corr_CHH <- computeMethylationProfile(epi2_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_2_corr_CHH <-computeMethylationDataSpatialCorrelation(epi2_corrected,
                                                              context="CHH",
                                                              distances=c(1,10,100,1000,10000))


profile_3_corr_CHH <- computeMethylationProfile(epi3_corrected,
                                                regions_loop,
                                                windowSize = 500000,
                                                context = "CHH")
spacorr_3_corr_CHH <-computeMethylationDataSpatialCorrelation(epi3_corrected,
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
####
#####
###########################################################

## Remove unnecessary files
rm(epi2_corrected, epi3_corrected, epi1_corrected)
gc()
dir.create("met_profiles")
setwd("met_profiles")

## Merge the plots into a single figure
allplot<- ((profileplot_CG | profileplot_CHG) / (profileplot_CHH | (correlationplot_CG | correlationplot_CHG) / (correlationplot_CHH | coverageplot))) + plot_layout(guides = 'collect')

ggsave(allplot, filename = "base_analysis_plots.png", width = 10, height = 10, device = "png", scale = 1.5)

setwd("../")











