## Some of the analyses are time consuming.
## Due to this reason, we stored the results by running the script
## The link of the stored data is provided in the manuscript 

## Required libraries for this study
library(treeio)
library(ggtree)
library(geiger)
library(ape)
library(phytools)
library(stringr)
library(ggrepel)
library(dplyr)
library(parallel)
library(digest)
library(magrittr)
library(tidyverse)
library(caper)

## Setting working directory s
setwd("/Users/admin/Desktop/nwork/TMRR/")
set.seed(23456)

## Functions required to run this script
source("functions_Dunn.R")  ## Original functions provided by Dunn et al.(2018)
source("functions_TM_new.R") ## Functions written for this study 

## Set system computational parameters
cores <- detectCores()

## Specifying the number of loops 
run<-100

## To perform analyses, we loaded our newly processed dataset, which were previously saved by running "Model_fitting.R" script
## That also includes the results of Dunn et al. (2018)
load("Dunn_tree_index_mapping_new.rda")

## To perform analyses, we filtered gene trees with negative branchlengths from the dataset of Dunn et al. (2018)
## Identifying calibrated trees with negative edge lengths
trees_negative <- negative_edgelength(gene_trees_calibrated)
## Removing trees with negative branch length to avoid problems in diagnostic tests 
gene_trees_pic1<-gene_trees_pic[-trees_negative]  

## Adding node heights (required for diagnostic tests) to the trees
gene_trees_pic1<-mclapply(gene_trees_pic1,tree_height,mc.cores=cores) 

## Identifying trees with at least one duplication, and one speciation events in Dunn et al. provided gene treeset 
Dunn_trees_with_duplication<-mclapply(gene_trees_pic1,tree_with_duplication, mc.cores = cores)
Dunn_trees_with_duplication<-Dunn_trees_with_duplication[!is.na(Dunn_trees_with_duplication)]    
Dunn_trees_with_duplication<-mclapply(Dunn_trees_with_duplication,tree_height,mc.cores=cores) ## Added node height to the trees

## Identifying trees with strong phylogenetic signals (following the criteria of Dunn et al. 2018))
## We make a vector to collect calibrated tree index with Blomberg's K > 0.551 
tree_index <-as.numeric(row.names(tree_summary)[(!(is.na(tree_summary$K)) & (tree_summary$K > 0.551))])
trees_new<- gene_trees_pic[c(tree_index)] ## 2082 trees 


###################################### Randomization tests on calibrated gene trees of Dunn et al.#########################################
############ Randomization of trait (Tau) 


## Initializing variables
pval_all<-NULL
pval_2082<-NULL

## Running the loop for 100 times
## We chose 100 times to save time but one can use higher value as the trend does not change
for(a in 1:run)
{
  #print(a)
  ## Initializing variable for each tree
  Tau_randomized_tree<-NULL 
  Tau_randomized_filtered_tree<-NULL
  
  ## We aim to shuffle or permute the trait data for 2 different sets of lists of trees
  Tau_randomized_tree<-mclapply(Dunn_trees_with_duplication,shuffling_tau,mc.cores = cores) ## Using list of 5479 trees with atleast one speciation and duplication event
  Tau_randomized_filtered_tree<-mclapply(trees_new,shuffling_tau,mc.cores = cores)  ## Using list of 2082 trees with strong phylogenetic signal
  
  ## Calculating PIC based on permuted Tau data for both sets of trees
  Tree_pic_randomized<-mclapply(Tau_randomized_tree,contrast_random,mc.cores = cores)
  Filtered_tree_pic_randomized_tau<-mclapply(Tau_randomized_filtered_tree,contrast_random,mc.cores = cores)
  
  ##Summarizing, and combining the data slots for all trees in a single dataframe
  ## Also removing rows with tip data for which contrasts are not available
  summary_tree_pic_randomized<-bind_rows(lapply(Tree_pic_randomized,summary_function))
  summary_filtered_tree_pic_randomized_tau<-bind_rows(lapply(Filtered_tree_pic_randomized_tau,summary_function))
  Random_tau_contrast<-summary_tree_pic_randomized[which(!is.na(summary_tree_pic_randomized$pic_abs_random)),]
  Random_tau_contrast_filtered_tree<-summary_filtered_tree_pic_randomized_tau[which(!is.na(summary_filtered_tree_pic_randomized_tau$pic_abs_random)),]
  
  ## This step is to remove rows with contrast values higher than 0.5 (if any) as used by Dunn et al.(2018)
  Random_tau_contrast<-Random_tau_contrast[which(Random_tau_contrast$pic_abs_random < 0.5),]
  Random_tau_contrast_filtered_tree<-Random_tau_contrast_filtered_tree[which(Random_tau_contrast_filtered_tree$pic_abs_random < 0.5),]
  
  ## Removing NA events, and collecting nodes contrast of all Tau randomized trees
  nodes_contrast_random<-Random_tau_contrast[!is.na(Random_tau_contrast$D),]
  nodes_contrast_random_tau_filtered_tree<-Random_tau_contrast_filtered_tree[!is.na(Random_tau_contrast_filtered_tree$D),]
  
  ##Separating PICs of speciation and duplication events to compare
  speciation_pic.r <- nodes_contrast_random$pic_abs_random[which(nodes_contrast_random$Event=="Speciation")]
  duplication_pic.r <- nodes_contrast_random$pic_abs_random[which(nodes_contrast_random$Event=="Duplication")]
  speciation_pic.rfilter <- nodes_contrast_random_tau_filtered_tree$pic_abs_random[which(nodes_contrast_random_tau_filtered_tree$Event=="Speciation")]
  duplication_pic.rfilter <- nodes_contrast_random_tau_filtered_tree$pic_abs_random[which(nodes_contrast_random_tau_filtered_tree$Event=="Duplication")]
  
  ## Performing two sided Wilcoxon rank test 
  wilcox_2_tailed.r <- wilcox.test(duplication_pic.r,speciation_pic.r)$p.value
  wilcox_2_tailed.rfilter <- wilcox.test(duplication_pic.rfilter,speciation_pic.rfilter)$p.value
  
  ## Collecting P values of 100 run for analyses
  pval_all<-append(pval_all,wilcox_2_tailed.r)
  pval_2082<-append(pval_2082,wilcox_2_tailed.rfilter)
}

################## Randomization of events ("Speciation"/"Duplication"/"NA")  
## Initializing variables
pval_all_events<-NULL
pval_2082_events<-NULL

## Running the loop for 100 times
for(b in 1:run)
{
  #print(b)
  ## Initializing variable for each tree
  Event_randomized_tree<-NULL 
  Event_randomized_filtered_tree<-NULL 
  
  ## We shuffle or permute the internal node events for two different sets of trees
  Event_randomized_tree<-mclapply(Dunn_trees_with_duplication, shuffling_event, mc.cores = cores) ## Using list of 5479 trees with atleast one speciation and duplication event
  Event_randomized_filtered_tree<-mclapply(trees_new,shuffling_event,mc.cores = cores) ## Using list of 2082 trees with strong phylogenetic signal
  
  ##Summarizing and combining the data slots of all trees in a single dataframe, and removing rows with tip data for which contrasts are not available
  summary_event_randomized_pic_tree<-bind_rows(lapply(Event_randomized_tree,summary_function))
  summary_filtered_tree_pic_randomized_event<-bind_rows(lapply(Event_randomized_filtered_tree,summary_function))
  Random_event_contrast<-summary_event_randomized_pic_tree[which(!is.na(summary_event_randomized_pic_tree$pic)),]
  Random_event_contrast_filtered_tree<-summary_filtered_tree_pic_randomized_event[which(!is.na(summary_filtered_tree_pic_randomized_event$pic)),]
  
  ## This step is to check and remove rows with contrast values higher than 0.5 (if any) as used by Dunn et al. 
  Random_event_contrast<-Random_event_contrast[which(Random_event_contrast$pic < 0.5),]
  Random_event_contrast_filtered_tree<-Random_event_contrast_filtered_tree[which(Random_event_contrast_filtered_tree$pic < 0.5),]
  
  ## Removing NA Events and collecting nodes contrasts of newly assigned speciation and duplication events after randomization
  nodes_contrast_random_events<-Random_event_contrast[!is.na(Random_event_contrast$event_new),]
  nodes_contrast_random_event_filtered_tree<-Random_event_contrast_filtered_tree[!is.na(Random_event_contrast_filtered_tree$event_new),]
  
  
  ##Separating PICs of newly assigned speciation and duplication events for both set of trees to compare
  speciation_pic_revent <-abs(nodes_contrast_random_events$pic[which(nodes_contrast_random_events$event_new=="Speciation")])
  duplication_pic_revent <- abs(nodes_contrast_random_events$pic[which(nodes_contrast_random_events$event_new=="Duplication")])
  speciation_pic_revent_2082 <- abs(nodes_contrast_random_event_filtered_tree$pic[which(nodes_contrast_random_event_filtered_tree$event_new=="Speciation")])
  duplication_pic_revent_2082 <- abs(nodes_contrast_random_event_filtered_tree$pic[which(nodes_contrast_random_event_filtered_tree$event_new=="Duplication")])
  
  ## Performing two sided Wilcoxon test, and collecting P values of 100 runs for plotting
  wilcox_2_tailed_revent <- wilcox.test(duplication_pic_revent,speciation_pic_revent)$p.value
  pval_all_events<-append(pval_all_events,wilcox_2_tailed_revent)
  wilcox_2_tailed_revent_2082 <- wilcox.test(duplication_pic_revent_2082,speciation_pic_revent_2082)$p.value
  pval_2082_events<-append(pval_2082_events,wilcox_2_tailed_revent_2082)
}

###################################### Performing diagnostic tests to check for adequate node contrasts standardization #################################
## For this diagnostic tests we used 4288 trees,for which calibration biases are fixed for old duplication nodes
## This diagnostic test function takes a lot of time
##It returns trees for which nodes contrasts are phylogenetically independent
count<-0
standardized_pic_tree<-lapply(trees_of_interest, diagnostic_plot_test) ## Using list of 4288 trees
standardized_pic_tree<-standardized_pic_tree[!is.na(standardized_pic_tree)] ## 2088 tree data passing the diagnostic tests

######################## Adequately contrasts standardized branch length transformed trees #################################

## Adding nodeheight to the sets of trees 
trees_of_interest<-mclapply(trees_of_interest,tree_height,mc.cores=cores) ## list of 4288 trees passed diagnostis tests
standardized_pic_tree<-mclapply(standardized_pic_tree,tree_height,mc.cores=cores) ## list of 2088 trees passed diagnostis tests
trees_new<-mclapply(trees_new,tree_height,mc.cores=cores) ## list of 2082 trees with strong phylogenetic signals
gene_trees_pic<-mclapply(gene_trees_pic,tree_height,mc.cores=cores) ## list of 8520 calibrated trees

##Performing diagnostic tests after branch length transformation
ct<-0
standardized_pic_tree_br_transformed<-lapply(standardized_pic_tree,branch_transform)
standardized_pic_tree_br_transformed<-standardized_pic_tree_br_transformed[ ! sapply(standardized_pic_tree_br_transformed, is.null) ]
standardized_pic_tree_br_transformed<-standardized_pic_tree_br_transformed[!is.na(standardized_pic_tree_br_transformed) ] ## list of 2088 trees passed diagnostic test after branch length transformation

ct<-0
trees_new_br_transformed<-lapply(trees_new,branch_transform)
trees_new_br_transformed<-trees_new_br_transformed[ ! sapply(trees_new_br_transformed, is.null) ]
trees_new_br_transformed<-trees_new_br_transformed[!is.na(trees_new_br_transformed) ] ## list of 2080 trees  passed diagnostic test after branch length transformation

ct<-0
gene_trees_pic_br_transformed<-lapply(gene_trees_pic,branch_transform)
gene_trees_pic_br_transformed<-gene_trees_pic_br_transformed[ ! sapply(gene_trees_pic_br_transformed, is.null) ]
gene_trees_pic_br_transformed<-gene_trees_pic_br_transformed[!is.na(gene_trees_pic_br_transformed) ] ## list of 8417 trees passed diagnostic test after branch length transformation

ct<-0
pic_tree_br_transformed<-lapply(trees_of_interest,branch_transform)
pic_tree_br_transformed<-pic_tree_br_transformed[ ! sapply(pic_tree_br_transformed, is.null) ]
pic_tree_br_transformed<-pic_tree_br_transformed[!is.na(pic_tree_br_transformed) ] ## list of 4190 trees passed diagnostic test after branch length transformation

#save.image("Analyses_TMRR.Rdata") 




