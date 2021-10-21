## Sankey plot
## Run after analysis one, IN THE SAME ENVIRONMENT
## It's separate from the file, since it requires manual exporting, and might not function on some workstations

library(plotly) ##requires ggplot2


## Acquire node size
## Object names explain contents
consensus_probeset<-length(bound_experiments[,1])
normalized_probeset<-length(normalized_data[,1])
probes_matched_with_genes<-sum(sort(table(assigned[,2]), decreasing = TRUE)[-c(1:2)])
probes_matched_without_genes<-sum(assigned[,2]=="No linked genes")
probes_matched_to_multiple_genes<-sum(assigned[,2]=="Multiple linked genes")
locations_assigned_TxDb<-length(GRange_compatible_data_tst[,1])
locations_assigned_noloc<-length(no_gene)
locations_assigned_fixed<-length(replacements)
locations_assigned_missing<-length(obsolete_removed)


single_probe_genes_probes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==1))
two_probe_genes_probes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==2)*2)
three_probe_genes_probes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==3)*3)
four_probe_genes_probes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==4)*4)
five_probe_genes_probes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==5)*5)


single_probe_genes_genes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==1))
two_probe_genes_genes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==2))
three_probe_genes_genes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==3))
four_probe_genes_genes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==4))
five_probe_genes_genes<-(sum(sort(table(ordered_GRange_data[,2]), decreasing = TRUE)==5))


## Construct the sankey plot itself
## To edit, open the object in the window, and move nodes around until desirable
Sankey_plot <- plot_ly(
  type = "sankey",
  orientation = "h",
  node = list(
    label = c(
              ),
    color = c("darkblue", "darkblue", "darkblue", "darkred", "darkred",
              "darkblue", "darkorange1", "darkblue", "darkorange1",
              "darkblue", "darkblue", "darkblue","darkblue", "darkblue",
              "darkblue", "darkblue", "darkblue","darkblue", "darkblue",
              "darkblue", "darkblue", "darkblue","darkblue", "darkblue"),
    pad = 15,
    thickness = 20,
    line = list(
      color = "grey",
      width = 0.1
    )
  ),
  link = list(
    source = c(0,1,1,1,2,
               2,6,6,7,
               5,5,5,5,5,
               10,11,12,13,14,
               15,16,17,18,19),
    target = c(1,2,3,4,5,
               6,7,8,5,
               10,11,12,13,14,
               15,16,17,18,19,
               20,20,20,20,20),
    value =  c(
      normalized_probeset, probes_matched_with_genes, probes_matched_without_genes, probes_matched_to_multiple_genes, locations_assigned_TxDb,
      locations_assigned_noloc, locations_assigned_fixed, locations_assigned_missing, locations_assigned_fixed,
      single_probe_genes_probes, two_probe_genes_probes, three_probe_genes_probes, four_probe_genes_probes, five_probe_genes_probes,
      single_probe_genes_genes, two_probe_genes_genes, three_probe_genes_genes, four_probe_genes_genes, five_probe_genes_genes,
      single_probe_genes_genes, two_probe_genes_genes, three_probe_genes_genes, four_probe_genes_genes, five_probe_genes_genes
    )
  )
) %>% 
  layout(
    font = list(
      size = 15
    )
  )


