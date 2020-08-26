
##### This script is written for fixing calibration bias for old duplicates, and to perform model-fitting
### To generate dataset with no or less calibration bias for old duplicates, we used the original scripts of Dunn et al. (2018), and modified few steps
### The original scripts of Dunn et al. (2018) are available at https://github.com/caseywdunn/comparative_expression_2017.
### We used our own scripts, and functions to reanalyze the data, and to perform model fitting in a maximum likelihood, and in a Bayesian frameworks


# The following should be installed from github with the specified 
# devtools command. You will need to install devtools first.
library(treeio) # devtools::install_github("GuangchuangYu/treeio")
library(ggtree) # devtools::install_github("GuangchuangYu/ggtree")
library(hutan) # devtools::install_github("caseywdunn/hutan")

# The rest can be installed from CRAN
library(ape)
library(digest)
library(geiger)
library(ggrepel)
library(gridExtra)
library(magrittr)
library(foreach)
library(doParallel)
library(phytools)
library(stringr)
library(tidyverse)
library(caper)
library(mvMORPH)
library(OUwie)
library(bayou)


## Setting working directory
setwd("/Users/admin/Desktop/nwork/simulation/")

## Functions required
source("functions_Dunn.R")	## Original function script provided by Dunn et al.(2018)
source("functions_TM_new.R") ## Functions written for this study
set.seed(23456)

## Setting system computational parameters
cores <- detectCores() - 1
if (cores < 1) {
  cores = 1
}

## Register parallel workers for %dopar%
registerDoParallel(cores)

## To test the ortholog conjecture simulation model, the fold change rates following duplication with respect to speciation
dup_adjust <- 2


## Compara gene trees file was downloaded from ftp://ftp.ensembl.org/pub/release-75/emf/ensembl-compara/homologies/
gene_tree_name <- "Compara.75.protein.nhx.emf"

## Tau values of species were taken from Rscript.R of Kryuchkova-Mostacci and Robinson-Rechavi (2016)
expression_file_names <- c(
  "ChickenBrawandTScomparisonTable_9_6chPC.txt", 
  "ChimpBrawandTScomparisonTable_9_6cmPC.txt", 
  "GorillaBrawandTScomparisonTable_9_6gPC.txt", 
  "HumBrawandTScomparisonTable_9_6hPC.txt", 
  "MacacaBrawandTScomparisonTable_9_6mcPC.txt", 
  "MusBrawandTScomparisonTable_9_6mPC.txt", 
  "OpossumBrawandTScomparisonTable_9_6oPC.txt", 
  "PlatypusBrawandTScomparisonTable_9_6pPC.txt"
)

## Directory for reanalyses
tmrr_directory <- "tmrr/"
expression_file_names %<>% paste(tmrr_directory, ., sep="")
gene_tree_name %<>% paste( tmrr_directory, ., sep="")

## The minimum number of genes with expression data for each gene tree, as used by Dunn et al. (2018)
min_genes_with_expression <- 4

## These clade dates are taken from the study of Kryuchkova-Mostacci and Robinson-Rechavi (2016), which are originally obtained from http://timetree.org/ 
calibration_times <- data.frame(
  age =
    c(20, 92, 167, 9, 42, 74, 296, 535, 104, 937, 29, 722, 441, 65, 
       15, 1215, 414, 371, 162, 25, 74, 77, 86, 7), 
  clade = 
    c("Hominoidea", "Euarchontoglires", "Mammalia", "Homininae", 
      "Simiiformes", "Primates", "Amniota", "Vertebrata", "Eutheria", 
      "Bilateria", "Catarrhini", "Chordata", "Euteleostomi", 
      "Haplorrhini", "Hominidae", "Opisthokonta", "Sarcopterygii", 
      "Tetrapoda", "Theria", "Murinae", "Sciurognathi", "Rodentia", 
      "Glires", "Hominini"
    ),
  stringsAsFactors=FALSE
)

## The 7 speciation time points used by Dunn et al.(2018)
focal_calibrations_clades <- 
  c("Hominini", "Homininae", "Catarrhini", "Euarchontoglires", "Theria", 
     "Mammalia", "Amniota")

## Reading the gene tree file as a vector of strings, one per line
lines <- readLines(gene_tree_name)
tree_lines <- lines[grepl( "^\\(", lines, perl = TRUE )]

## Parse each string to a treeio::treedata object, adding a label column for node names
gene_trees <- foreach(tree_line=tree_lines) %dopar%
  parse_gene_trees(tree_line)

## Total number of gene trees
n_gene_trees <- length(gene_trees)

## We added index to all the gene trees 
Index<-0
for(Index in 1:length(gene_trees))
{
  gene_trees[[Index]]@data$index.tree<-Index
}

## Parsing the tree annotations, these are lines that start with SEQ
tip_annotations <- 
  lines[grepl("^SEQ", lines, perl = TRUE)] %>%	
  str_replace(" \\(\\d+ of \\d+\\)", "") %>%
  str_replace("(([^ ]+ ){8})([^ ]+) ([^ ]+)", "\\1\\3_\\4") %>% 
  str_c(collapse="\n") %>%
  read_delim(delim=" ", col_names = FALSE) %>% 	
  rename(Ensembl.Gene.ID = X8) %>%
  rename(species = X2)

## Selecting required columns
tip_annotations<- tip_annotations[c(2,8,9)]

## Free up memory
rm(lines)
rm(tree_lines)

## Reading the expression data for each species and combine into a single tibble,
expression <- 
  lapply(
    expression_file_names, 
    read.table, 
    stringsAsFactors=FALSE, 
    header=TRUE
  ) %>%
  bind_rows() %>%
  left_join(tip_annotations, by="Ensembl.Gene.ID")

## Adding the expression data to the tips of the trees
## Annotate each tree by joining the corresponding expression data to the 
## @data object. 

add_expression_to_tree <- function(tree){
  tree@data %<>% 
    left_join(expression, by = c("G" = "Ensembl.Gene.ID")) 
  return(tree)
}
gene_trees_annotated <- foreach(tree=gene_trees) %dopar% 
  add_expression_to_tree(tree)

## Cleaning memory
#rm(gene_trees)


############ New ways to time calibrate trees fixed by TMRR ###########

## Fixing clade names following Dunn et al. (2018) before calibration
gene_trees_annotated <- lapply(gene_trees_annotated,fix_hominini)
calibration_times$label<- calibration_times$clade

## Modifyting duplication node labels
#gene_trees_annotated1 = lapply(gene_trees_annotated, modify_label)
gene_trees_annotated1 <- lapply(gene_trees_annotated, tree_nodedepth)

## Time calibrating the trees before pruning unlike Dunn et al. (2018)
count<-0
Time_calibrated_trees <- lapply(gene_trees_annotated1, 
                                tree_calibrate,
                                timeframe=calibration_times,
                                model="correlated")

## Successfully calibrated trees 
gene_trees_calibrated_new<-Time_calibrated_trees[!is.na(Time_calibrated_trees)] 
gene_trees_calibrated_new<-gene_trees_calibrated_new[!sapply(gene_trees_calibrated_new, is.null)]

## Pruning trees, i.e., removing tips from the trees that do not have expression data
gene_trees_prune_new <- foreach( tree=gene_trees_calibrated_new ) %dopar% 
  drop_empty_tips( 
    tree, 
    min_genes_with_expression=min_genes_with_expression 
  )
gene_trees_prune_new<-gene_trees_prune_new[!is.na(gene_trees_prune_new)] 
gene_trees_prune_new<-gene_trees_prune_new[!sapply(gene_trees_prune_new, is.null)]##passed 7336 trees 

## Remodifying the previously modified duplication label 
gene_trees_prune_new <- lapply(gene_trees_prune_new, remodify_label)

## Storing node age using the function provided by Dunn et al. (2018)
gene_trees_prune_calibrated <- foreach(tree=gene_trees_prune_new) %dopar% 
  store_node_age(tree)

gene_trees_calibrated_now<-gene_trees_prune_calibrated ##passed 7336 trees

## Remove trees with no speciation events, i.e., pure duplication trees
gene_trees_calibrated_now <- 
  gene_trees_calibrated_now[ 
    unlist( lapply( gene_trees_calibrated_now, get_n_speciation ) ) > 0 
    ] ## Passed 7297 trees, and removed 39 (7336-7297) pure duplication trees

## Calculating phylogenetic independent contrasts using the function provided by Dunn et al.(2018)
gene_trees_pic_new <- add_pics_to_trees(gene_trees_calibrated_now) ##passed 7297 trees

## Adding node heights to each gene tree
gene_trees_pic_new<-lapply(gene_trees_pic_new,tree_height)
summary_all_trees<-bind_rows(mclapply(gene_trees_pic_new,summary_function,mc.cores = cores))
summary_all_trees<-summary_all_trees[which(!is.na(summary_all_trees$pic)),]
length(summary_all_trees$pic[which(summary_all_trees$Event=="Duplication")]) ##n=15274
length(summary_all_trees$pic[which(summary_all_trees$Event=="Speciation")]) ##n=52081
## Total 52081 speciation events, and 15274 duplication events

## Cleaning memory
rm(Time_calibrated_trees)
rm(gene_trees_prune_calibrated) 
rm(gene_trees_annotated1)
rm(expression)
rm(tip_annotations)
#save.image("Dunn_tree_index_mapping_new.rda") 

## Keeping gene trees with at least one duplication, and one speciation events
## This is necessary to understand the evolutionary trajectories of different events for the same gene
trees_of_interest<-mclapply(gene_trees_pic_new,tree_with_duplication, mc.cores = cores)
trees_of_interest<-trees_of_interest[!is.na(trees_of_interest)]   
trees_of_interest<-lapply(trees_of_interest,tree_height) #passed 4288 trees

## Median number of tips for the retained gene trees
tree_info<-bind_rows(lapply(trees_of_interest, tree_data_statistics,"all"))
median(tree_info$tip_num) ## n= 15

## Summmary statistics of events 
trees_of_interest_PIC <- add_pics_to_trees( trees_of_interest ) ## 4288 trees
summary_4288_trees<-bind_rows(mclapply(trees_of_interest_PIC,summary_function,mc.cores = cores))
summary_4288_trees<-summary_4288_trees[which(!is.na(summary_4288_trees$pic)),]
length(summary_4288_trees$pic[which(summary_4288_trees$Event=="Duplication")]) ##n=15274 duplication
length(summary_4288_trees$pic[which(summary_4288_trees$Event=="Speciation")])  ##n=38882 speciation
length(summary_4288_trees$pic[which(is.na(summary_4288_trees$Event))]) ##n=15201 NA event

## Total 38882 speciation events, and 15274 duplication events (i.e. 13199 speciation events are removed by removal of pure speciation trees)

# Loaded the previously stored data generated by using "manuscript_kernel.R" of Dunn et al.(2018) to compare
load("/Users/admin/Desktop/nwork/simulation/manuscript_dunn.RData") 

################### Model fitting using mvMORPH, and OUwie after painting speciation, young and old duplication branches separately for each tree #####################

trees_of_interest_modified<-lapply(trees_of_interest,modify_tiplabel)  # Changing tip Ensembl protein ID to Ensembl Gene ID

## Painting trees before model fitting
count<-0
trees_of_interest_painted<-lapply(trees_of_interest_modified,paint_tree_mod_new,"Empirical") 
trees_of_interest_painted1<- trees_of_interest_painted[!sapply(trees_of_interest_painted, is.null)]
trees_of_interest_painted1<-trees_of_interest_painted1[!is.na(trees_of_interest_painted1)]


#################### Model fitting analyses #################### 
i<-0
Tree_model_fit_summary<-bind_rows(lapply(trees_of_interest_painted1,model_fitting_genetree_for_Tau,trait="Tau"))

## Adding model parameters to the "phylo" slots of the trees, collected from the best fit model 
i<-0
Tree_added_parameters<-lapply(trees_of_interest_painted1,add_model_parameters_info,Tree_model_fit_summary)
Tree_added_parameters<-mclapply(Tree_added_parameters,tree_added_events_info,mc.cores=cores)

### Counting total speciation and duplication events of 4288 trees
## Classifying duplication events into young (age <= 296 My, maximum speciation age), and old (age > 296 My)  for further analyses
Tree_all<-Tree_added_parameters[c(1:4288)]
table_event<-NULL
for(i in 1:length(Tree_all))
{
  table_event<-bind_rows(table_event,tibble(index=i,
                                            n_speciation=Tree_all[[i]]@phylo$n_speciation,
                                            n_duplication_young=Tree_all[[i]]@phylo$n_duplication_young,
                                            n_duplication_old=Tree_all[[i]]@phylo$n_duplication_old))
}

sum(table_event$n_speciation) #n=38882
sum(table_event$n_duplication_young) #n=8556
sum(table_event$n_duplication_old) ##n=6718
## Total 38882 speciation events, and 15274 (8556 young + 6718 old) duplication events


## Out of 4288 trees, 1630 trees have both young and old duplication events
table_event_shared<-table_event[which(table_event$n_duplication_young > 0 & table_event$n_duplication_old > 0),]
sum(table_event_shared$n_speciation) # n=19698
sum(table_event_shared$n_duplication_young) # n=4653
sum(table_event_shared$n_duplication_old) ## n=4387
## Total 19698 speciation events, and 9040 (4653 young + 4387 old) duplication events for trees having both type of duplicates

## Summary statistics of model fitting by maximum likelihood approach
# BM1:308,
# BMM:704,
# OU1:2874,
# OUM:370,
# 32 trees failed to fit any model 

## save.image("Model_fitting_TMRR.Rdata") 

################ Adaptive optimum (θ) estimation for OUM trees ###########

index_OUM<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model %in% c("OUM"))] 
Tree_OUM<-Tree_added_parameters[c(index_OUM)]##370 trees

table_OUM<-NULL
for(i in 1:length(Tree_OUM))
{
  table_OUM<-bind_rows(table_OUM,tibble(index=i,spe_evol_rate=Tree_OUM[[i]]@phylo$fitted_parameters$sig2_speciation,
                                        spe_alpha=Tree_OUM[[i]]@phylo$fitted_parameters$alpha_speciation,
                                        theta_spe=Tree_OUM[[i]]@phylo$fitted_parameters$theta_speciation,
                                        dupy_evol_rate=Tree_OUM[[i]]@phylo$fitted_parameters$sig2_duplication_young,
                                        dupy_alpha=Tree_OUM[[i]]@phylo$fitted_parameters$alpha_duplication_young,
                                        theta_dupy=Tree_OUM[[i]]@phylo$fitted_parameters$theta_duplication_young,
                                        dupo_evol_rate=Tree_OUM[[i]]@phylo$fitted_parameters$sig2_duplication_old,
                                        dupo_alpha=Tree_OUM[[i]]@phylo$fitted_parameters$alpha_duplication_old,
                                        theta_dupo=Tree_OUM[[i]]@phylo$fitted_parameters$theta_duplication_old,
                                        n_speciation=Tree_OUM[[i]]@phylo$n_speciation,
                                        n_duplication_young=Tree_OUM[[i]]@phylo$n_duplication_young,
                                        n_duplication_old=Tree_OUM[[i]]@phylo$n_duplication_old))
  
}

## Analyses for young duplication 
table_OUM_with_young_dup<-table_OUM[c(2:7,11:12)]
table_OUM_with_young_dup<-table_OUM_with_young_dup[complete.cases(table_OUM_with_young_dup),]
median(table_OUM_with_young_dup$theta_spe) ## 0.4075981 (n=2690)
median(table_OUM_with_young_dup$theta_dupy) ##0.7387739 (n=842)
wilcox.test(table_OUM_with_young_dup$theta_spe,table_OUM_with_young_dup$theta_dupy,paired = T) #p-value = 8.591e-10

## Analyses for old duplication 
table_OUM_with_old_dup<-table_OUM[c(2:4,8:10,11,13)]
table_OUM_with_old_dup<-table_OUM_with_old_dup[complete.cases(table_OUM_with_old_dup),]
median(table_OUM_with_old_dup$theta_spe) ## 0.4221139 (n=4152)
median(table_OUM_with_old_dup$theta_dupo) ##0.9175527 (n=847)
wilcox.test(table_OUM_with_old_dup$theta_spe,table_OUM_with_old_dup$theta_dupo,paired = T) # p-value = 0.0002459


################ Multi rate Brownian trees (BMM) for evolutionary rate (σ2) analyses ###########
index_BMM<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model%in% c("BMM"))]
Tree_BMM<-Tree_added_parameters[c(index_BMM)] ## 704 trees
#bmPlot(Tree_BMM[[2]]@phylo$painted_tree,type="threshold",thresholds=c(0,1))
table_BMM<-NULL
for(i in 1:length(Tree_BMM))
{
  table_BMM<-bind_rows(table_BMM,tibble(index=i,spe_evol_rate=Tree_BMM[[i]]@phylo$fitted_parameters$sig2_speciation,
                                        dupy_evol_rate=Tree_BMM[[i]]@phylo$fitted_parameters$sig2_duplication_young,
                                        dupo_evol_rate=Tree_BMM[[i]]@phylo$fitted_parameters$sig2_duplication_old,
                                        n_speciation=Tree_BMM[[i]]@phylo$n_speciation,
                                        n_duplication_young=Tree_BMM[[i]]@phylo$n_duplication_young,
                                        n_duplication_old=Tree_BMM[[i]]@phylo$n_duplication_old))
  
}

##Analyses for young duplication
table_BMM_with_young_dup<-table_BMM[c(2:3,5:6)]
table_BMM_with_young_dup<-table_BMM_with_young_dup[complete.cases(table_BMM_with_young_dup),]
median(table_BMM_with_young_dup$spe_evol_rate) ## 0.00009075896 (n=4642)
median(table_BMM_with_young_dup$dupy_evol_rate) ##0.0001435775 (n=1742)
wilcox.test(table_BMM_with_young_dup$spe_evol_rate,table_BMM_with_young_dup$dupy_evol_rate,paired = T) #p-value = 5.044e-12

##Analyses for old duplication
table_BMM_with_old_dup<-table_BMM[c(2,4,5,7)]
table_BMM_with_old_dup<-table_BMM_with_old_dup[complete.cases(table_BMM_with_old_dup),]
median(table_BMM_with_old_dup$spe_evol_rate) ##0.0001734904 (n=5356)
median(table_BMM_with_old_dup$dupo_evol_rate) ##2.061154e-09 (n=1295)
wilcox.test(table_BMM_with_old_dup$spe_evol_rate,table_BMM_with_old_dup$dupo_evol_rate,paired = T) #p-value < 2.2e-16


############################ Randomization tests ####################################
set.seed(23456)

################# Analysis on event randomized tree ################

## First we shuffle the internal node events 
Randomized_Event_trees<- mclapply(trees_of_interest_painted,shuffling_event,mc.cores = cores)

## Repainting the event randomized tree according to the newly assigned internal events
count<-0
Event_randomized_trees_painted<-lapply(Randomized_Event_trees,paint_tree_mod_new,"Random")

################### Event randomized BMM trees###################

## Collecting BM trees with randomized events 
df_BMM<-Tree_model_fit_summary[which(Tree_model_fit_summary$best_fit_model %in% c("BMM")),]
index_BMM<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model%in% c("BMM"))]
Tree_randomized_event_BMM<-Event_randomized_trees_painted[c(index_BMM)]## 704 trees

## Adding model parameters by fitting the BMM model
i<-0
Randomized_BMM_tree_added_parameter<-lapply(Tree_randomized_event_BMM,add_model_parameters_info2,df_BMM)
Randomized_BMM_tree_added_parameter2<-Randomized_BMM_tree_added_parameter[!is.na(Randomized_BMM_tree_added_parameter)]
Randomized_BMM_tree_added_parameter3<-mclapply(Randomized_BMM_tree_added_parameter2,tree_added_events_info,mc.cores=cores)

## Analyses
table_randomized_event_BMM<-NULL
for(x in 1:length(Randomized_BMM_tree_added_parameter3))
{
  print(x)
  table_randomized_event_BMM<-bind_rows(table_randomized_event_BMM,tibble(index=x,spe_evol_rate=Randomized_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_speciation,
                                                                          dupy_evol_rate=Randomized_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_duplication_young,
                                                                          dupo_evol_rate=Randomized_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_duplication_old,
                                                                          n_speciation=Randomized_BMM_tree_added_parameter3[[x]]@phylo$n_speciation,
                                                                          n_duplication_young=Randomized_BMM_tree_added_parameter3[[x]]@phylo$n_duplication_young,
                                                                          n_duplication_old=Randomized_BMM_tree_added_parameter3[[x]]@phylo$n_duplication_old))
  
}

## Analyses on young duplicates
table_randomized_event_BMM_young_dup<-table_randomized_event_BMM[c(2:3,5:6)]
table_randomized_event_BMM_young_dup<-table_randomized_event_BMM_young_dup[complete.cases(table_randomized_event_BMM_young_dup),]
median(table_randomized_event_BMM_young_dup$spe_evol_rate) ##  0.0001663995 (n=3215)
median(table_randomized_event_BMM_young_dup$dupy_evol_rate) ## 0.00008534398 (n=1438)
wilcox.test(table_randomized_event_BMM_young_dup$spe_evol_rate,table_randomized_event_BMM_young_dup$dupy_evol_rate,paired = T) #p-value = 0.02783


## Analyses on old duplicates
table_randomized_event_BMM_old_dup<-table_randomized_event_BMM[c(2,4,5,7)]
table_randomized_event_BMM_old_dup<-table_randomized_event_BMM_old_dup[complete.cases(table_randomized_event_BMM_old_dup),]
median(table_randomized_event_BMM_old_dup$spe_evol_rate) ##   0.0001762469 (n=2788)
median(table_randomized_event_BMM_old_dup$dupo_evol_rate) ## 2.061154e-09 (n=800)
wilcox.test(table_randomized_event_BMM_old_dup$spe_evol_rate,table_randomized_event_BMM_old_dup$dupo_evol_rate,paired = T) #p-value < 2.2e-16


################### Event randomized OUM trees###################

## Collecting OUM trees with randomized events using the index of the OUM trees
df_OUM<-Tree_model_fit_summary[which(Tree_model_fit_summary$best_fit_model %in% c("OUM")),]
index_OU<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model%in% c("OUM"))]
Tree_randomized_event_OUM<-Event_randomized_trees_painted[c(index_OUM)]## 370 trees

## Adding model parameters by fitting the OUM model
i<-0
Randomized_OUM_tree_added_parameter<-lapply(Tree_randomized_event_OUM,add_model_parameters_info2,df_OUM)
Randomized_OUM_tree_added_parameter2<-Randomized_OUM_tree_added_parameter[!is.na(Randomized_OUM_tree_added_parameter)]
Randomized_OUM_tree_added_parameter3<-mclapply(Randomized_OUM_tree_added_parameter2,tree_added_events_info,mc.cores=cores)

## Analyses
table_OUM_random<-NULL
for(i in 1:length(Randomized_OUM_tree_added_parameter3))
{
  table_OUM_random<-bind_rows(table_OUM_random,tibble(index=i,spe_alpha=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_speciation,
                                                      spe_evol_rate=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_speciation,
                                                      theta_spe=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_speciation,
                                                      dupy_alpha=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_duplication_young,
                                                      dupy_evol_rate=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_duplication_young,
                                                      theta_dupy=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_duplication_young,
                                                      dupo_alpha=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_duplication_old,
                                                      dupo_evol_rate=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_duplication_old,
                                                      theta_dupo=Randomized_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_duplication_old,
                                                      n_speciation=Randomized_OUM_tree_added_parameter3[[i]]@phylo$n_speciation,
                                                      n_duplication_young=Randomized_OUM_tree_added_parameter3[[i]]@phylo$n_duplication_young,
                                                      n_duplication_old=Randomized_OUM_tree_added_parameter3[[i]]@phylo$n_duplication_old))
  
}

## With young duplication
table_OUM_random_with_young_dup<-table_OUM_random[c(2:7,11:12)]
table_OUM_random_with_young_dup<-table_OUM_random_with_young_dup[complete.cases(table_OUM_random_with_young_dup),]
median(table_OUM_random_with_young_dup$theta_spe) ##   0.5047096 (n=1872)
median(table_OUM_random_with_young_dup$theta_dupy) ## 0.5339038 (n=698)
wilcox.test(table_OUM_random_with_young_dup$theta_spe,table_OUM_random_with_young_dup$theta_dupy,paired = T) #p-value = 0.7483


## With old duplication
table_OUM_random_with_old_dup<-table_OUM_random[c(2:4,8:10,11,13)]
table_OUM_random_with_old_dup<-table_OUM_random_with_old_dup[complete.cases(table_OUM_random_with_old_dup),]
median(table_OUM_random_with_old_dup$theta_spe) ##   0.5093067 (n=2081)
median(table_OUM_random_with_old_dup$theta_dupo) ## 0.660308 (n=482)
wilcox.test(table_OUM_random_with_old_dup$theta_spe,table_OUM_random_with_old_dup$theta_dupo,paired = T) #p-value = 0.7339

#save.image("Model_fitting_TMRR.Rdata")

################# Analysis on Tau randomized tree ################
set.seed( 23456 )

## First we shuffle the tip data 
Randomized_Tau_trees<- mclapply(trees_of_interest_painted,shuffling_tauR,mc.cores = cores)

################### Tau randomized BMM trees###################

## Collecting BM trees with randomized Tau using the index of the BM trees
df_BMM<-Tree_model_fit_summary[which(Tree_model_fit_summary$best_fit_model %in% c("BMM")),]
index_BMM<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model%in% c("BMM"))]
Tree_randomized_Tau_BMM<-Randomized_Tau_trees[c(index_BMM)]## 704 trees

## Adding model parameters by fitting the BMM trees
i<-0
Randomized_Tau_BMM_tree_added_parameter<-lapply(Tree_randomized_Tau_BMM,add_model_parameters_info2,df_BMM)
Randomized_Tau_BMM_tree_added_parameter2<-Randomized_Tau_BMM_tree_added_parameter[!is.na(Randomized_Tau_BMM_tree_added_parameter)]
Randomized_Tau_BMM_tree_added_parameter3<-mclapply(Randomized_Tau_BMM_tree_added_parameter2,tree_added_events_info,mc.cores=cores)

## Analyses
table_randomized_Tau_BMM<-NULL
for(x in 1:length(Randomized_Tau_BMM_tree_added_parameter3))
{
  print(x)
  table_randomized_Tau_BMM<-bind_rows(table_randomized_Tau_BMM,tibble(index=x,spe_evol_rate=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_speciation,
                                                                      dupy_evol_rate=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_duplication_young,
                                                                      dupo_evol_rate=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$fitted_parameters$sig2_duplication_old,
                                                                      n_speciation=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$n_speciation,
                                                                      n_duplication_young=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$n_duplication_young,
                                                                      n_duplication_old=Randomized_Tau_BMM_tree_added_parameter3[[x]]@phylo$n_duplication_old))
  
}

## Analyses on young duplicates
table_randomized_Tau_BMM_young_dup<-table_randomized_Tau_BMM[c(2:3,5:6)]
table_randomized_Tau_BMM_young_dup<-table_randomized_Tau_BMM_young_dup[complete.cases(table_randomized_Tau_BMM_young_dup),]
median(table_randomized_Tau_BMM_young_dup$spe_evol_rate) ##   0.0006957318 (n=4618)
median(table_randomized_Tau_BMM_young_dup$dupy_evol_rate) ## 0.0002236885 (n=1723)
wilcox.test(table_randomized_Tau_BMM_young_dup$spe_evol_rate,table_randomized_Tau_BMM_young_dup$dupy_evol_rate,paired = T,) #p-value = 1.388e-13


## Analyses on old duplicates

table_randomized_Tau_BMM_old_dup<-table_randomized_Tau_BMM[c(2,4,5,7)]
table_randomized_Tau_BMM_old_dup<-table_randomized_Tau_BMM_old_dup[complete.cases(table_randomized_Tau_BMM_old_dup),]
median(table_randomized_Tau_BMM_old_dup$spe_evol_rate) ##0.000908805 (n=5337)
median(table_randomized_Tau_BMM_old_dup$dupo_evol_rate) ##2.529964e-10 (n=1291)
wilcox.test(table_randomized_Tau_BMM_old_dup$spe_evol_rate,table_randomized_Tau_BMM_old_dup$dupo_evol_rate,paired = T) #p-value < 2.2e-16

############## Tau randomized OUM trees ##############

## Collecting OU trees with randomized events using the index of the OU trees
df_OUM<-Tree_model_fit_summary[which(Tree_model_fit_summary$best_fit_model %in% c("OUM")),]
index_OUM<-Tree_model_fit_summary$tree_index[which(Tree_model_fit_summary$best_fit_model%in% c("OUM"))]
Tree_randomized_Tau_OUM<-Randomized_Tau_trees[c(index_OUM)]## 370 trees

## Adding model parameters by fitting the OUM trees
i<-0
Randomized_Tau_OUM_tree_added_parameter<-lapply(Tree_randomized_Tau_OUM,add_model_parameters_info2,df_OUM)
Randomized_Tau_OUM_tree_added_parameter2<-Randomized_Tau_OUM_tree_added_parameter[!is.na(Randomized_Tau_OUM_tree_added_parameter)]
Randomized_Tau_OUM_tree_added_parameter3<-mclapply(Randomized_Tau_OUM_tree_added_parameter2,tree_added_events_info,mc.cores=cores)

table_OUM_Tau_random<-NULL
for(i in 1:length(Randomized_Tau_OUM_tree_added_parameter3))
{
  table_OUM_Tau_random<-bind_rows(table_OUM_Tau_random,tibble(index=i,spe_alpha=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_speciation,
                                                              spe_evol_rate=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_speciation,
                                                              theta_spe=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_speciation,
                                                              dupy_alpha=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_duplication_young,
                                                              dupy_evol_rate=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_duplication_young,
                                                              theta_dupy=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_duplication_young,
                                                              dupo_alpha=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$alpha_duplication_old,
                                                              dupo_evol_rate=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$sig2_duplication_old,
                                                              theta_dupo=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$fitted_parameters$theta_duplication_old,
                                                              n_speciation=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$n_speciation,
                                                              n_duplication_young=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$n_duplication_young,
                                                              n_duplication_old=Randomized_Tau_OUM_tree_added_parameter3[[i]]@phylo$n_duplication_old))
  
}

## With young duplication
table_OUM_Tau_random_with_young_dup<-table_OUM_Tau_random[c(2:7,11:12)]
table_OUM_Tau_random_with_young_dup<-table_OUM_Tau_random_with_young_dup[complete.cases(table_OUM_Tau_random_with_young_dup),]
median(table_OUM_Tau_random_with_young_dup$theta_spe) ##   0.5292699 (n=2690)
median(table_OUM_Tau_random_with_young_dup$theta_dupy) ##  0.5524303 (n=842)
wilcox.test(table_OUM_Tau_random_with_young_dup$theta_spe,table_OUM_Tau_random_with_young_dup$theta_dupy,paired = T) #p-value = 0.9754


## With old duplication
table_OUM_Tau_random_with_old_dup<-table_OUM_Tau_random[c(2:4,8:10,11,13)]
table_OUM_Tau_random_with_old_dup<-table_OUM_Tau_random_with_old_dup[complete.cases(table_OUM_Tau_random_with_old_dup),]
median(table_OUM_Tau_random_with_old_dup$theta_spe) ## 0.5391824 (n=4152)
median(table_OUM_Tau_random_with_old_dup$theta_dupo) ## 1.46883e-11 (n=847)
wilcox.test(table_OUM_Tau_random_with_old_dup$theta_spe,table_OUM_Tau_random_with_old_dup$theta_dupo,paired = T) #p-value = 1.355e-07

#save.image("Model_fitting_TMRR.Rdata")

####### Cleaning memory 
rm(lines)
rm(tree_lines)
rm(i)
rm(x)
rm(Index)
rm(index_BMM)
rm(index_OU)
rm(index_OUM)
rm(df_BMM)
rm(df_OUM)
rm(Randomized_BMM_tree_added_parameter2)
rm(Randomized_BMM_tree_added_parameter3)
rm(Randomized_OUM_tree_added_parameter2)
rm(Randomized_OUM_tree_added_parameter3)
rm(Randomized_Tau_BMM_tree_added_parameter2)
rm(Randomized_Tau_BMM_tree_added_parameter3)
rm(Randomized_Tau_OUM_tree_added_parameter2)
rm(Randomized_Tau_OUM_tree_added_parameter3)

#save.image("Model_fitting_TMRR.Rdata")


################################## Baysian model fitting to identify multi peak OU or OUM trees #################
## We used bayou R package for this purpose
## We consider a strict posterior probability cutoff of >= 0.7 to identify regime, i.e. θ shift

## Considering all OU Trees
count<-0
result_bayou_OU<-NULL
result_bayou_OU<-bind_rows(lapply(Tree_OU,bayou_analysis))

## Considering OUM Trees
count<-0
result_bayou_OUM<-NULL
result_bayou_OUM<-bind_rows(lapply(Tree_OUM,bayou_analysis))

##################### Optima shift analyses on 370 OUM trees ###################

## Analysis considering all branches of the OU trees
result_bayou_OUM_new<-result_bayou_OUM[result_bayou_OUM$pp >= 0.7,]

length(unique(result_bayou_OUM_new$index)) ## 206 out of 370 (55.7%)
Baysian_index<-unique(result_bayou_OUM_new$index)
Trees_OUM_new<-Tree_OUM[c(Baysian_index)]
rm(Baysian_index)

## To obtain total numbers of events
table_OUM_0.7pp<-NULL
for(i in 1:length(Trees_OUM_new))
{
  table_OUM_0.7pp<-bind_rows(table_OUM_0.7pp,tibble(index=i,
                                                      n_speciation=Trees_OUM_new[[i]]@phylo$n_speciation,
                                                      n_duplication_young=Trees_OUM_new[[i]]@phylo$n_duplication_young,
                                                      n_duplication_old=Trees_OUM_new[[i]]@phylo$n_duplication_old))
  
}
rm(i)
sum(table_OUM_0.7pp$n_speciation) ## 2779 speciation events
sum(table_OUM_0.7pp$n_duplication_young) ## 486 young duplication events
sum(table_OUM_0.7pp$n_duplication_old) ## 548 old duplication events

## Regime shift analyses considering the whole tree
OUM_all<-proportion_calculation(Tree_OUM,"all",result_bayou_OUM_new)

median(na.omit(OUM_all$shift_prob_dup)) #0.125
median(na.omit(OUM_all$shift_prob_spe)) #0
mean(na.omit(OUM_all$shift_prob_dup)) #0.509
mean(na.omit(OUM_all$shift_prob_spe)) #0.232
pvalue_regime_shift_all_OUM<- wilcox.test(OUM_all$shift_prob_spe,OUM_all$shift_prob_dup,paired = T)$p.value #p-value = 1.45e-20

## Regime shift analysis considering young duplication and speciation branches of the OUM trees
result_bayoum1<-result_bayou_OUM_new[result_bayou_OUM_new$node_age<=296,]
OUM_young<-proportion_calculation(Tree_OUM,"young",result_bayoum1)

median(na.omit(OUM_young$shift_prob_dup)) #0.083
median(na.omit(OUM_young$shift_prob_spe)) #0.028
mean(na.omit(OUM_young$shift_prob_dup)) #0.138
mean(na.omit(OUM_young$shift_prob_spe)) #0.05
pvalue_regime_shift_OUM_young<- wilcox.test(OUM_young$shift_prob_spe,OUM_young$shift_prob_dup,paired = T)$p.value #p-value = 8.26e-4

## Regime shift analysis considering old duplication and speciation branches of the OUM trees
result_bayoum2<-NULL
result_bayou_OUM_spe<-result_bayou_OUM_new[result_bayou_OUM_new$Event=="Speciation",]
result_bayou_OUM_dup_old<-result_bayou_OUM_new[result_bayou_OUM_new$Event=="Duplication" & result_bayou_OUM_new$node_age > 296,]
result_bayoum2<-rbind(result_bayou_OUM_spe,result_bayou_OUM_dup_old)

OUM_old<-proportion_calculation(Tree_OUM,"old",result_bayoum2)

median(na.omit(OUM_old$shift_prob_dup)) #0.167
median(na.omit(OUM_old$shift_prob_spe)) #0
mean(na.omit(OUM_old$shift_prob_dup)) #0.173
mean(na.omit(OUM_old$shift_prob_spe)) #0.028
pvalue_regime_shift_OUM_old<- wilcox.test(OUM_old$shift_prob_spe,OUM_old$shift_prob_dup,paired = T)$p.value #p-value = 7.26e-40

## Calculating the shift rates/My for OUM trees

## Step 1: Estimating length of the branch where optimum shift was detected
for(i in 1:nrow(result_bayou_OUM_new))
{
  
  tre<-NULL
  index<-NULL
  edgelength<-NULL
  index<-result_bayou_OUM_new$index[i]
  tre<-Tree_OUM[[index]]@phylo
  result_bayou_OUM_new$edgelength[i]<-tre$edge.length[which(tre$edge[,1]==result_bayou_OUM_new$node1[i] & tre$edge[,2]==result_bayou_OUM_new$node2[i])]
}

## Step 2: Based on the relative location provided by Bayesian approach, we estimated the time when optimum shifts took place 
result_bayou_OUM_new$shift_age<-result_bayou_OUM_new$node_age-(result_bayou_OUM_new$edgelength*result_bayou_OUM_new$rel.location)

## Shift rates calculation
result_bayou_OUM_new$shifts_per_my<-1/result_bayou_OUM_new$shift_age

## For young duplicates
result_bayou_OUM_young<-result_bayou_OUM_new[result_bayou_OUM_new$node_age <= 296,]
optima_shifts_per_my_young_duplicates_OUM<-result_bayou_OUM_young$shifts_per_my[which(result_bayou_OUM_young$Event=="Duplication")]
optima_shifts_per_my_speciation_OUM<- result_bayou_OUM_young$shifts_per_my[which(result_bayou_OUM_young$Event=="Speciation")]
median(optima_shifts_per_my_young_duplicates_OUM) #0.032
median(optima_shifts_per_my_speciation_OUM) ##0.0121
wilcox.test(optima_shifts_per_my_young_duplicates_OUM,optima_shifts_per_my_speciation_OUM) #p-value = 2.341e-11

## For old duplicates
result_bayou_OUM_old<-NULL
result_bayou_OUM_old<-result_bayou_OUM_new[result_bayou_OUM_new$Event=="Speciation",]
result_bayou_OUM_old<-rbind(result_bayou_OUM_old,result_bayou_OUM_new[result_bayou_OUM_new$Event=="Duplication" & result_bayou_OUM_new$node_age > 296,])
result_bayou_OUM_old<-na.omit(result_bayou_OUM_old)
optima_shifts_per_my_old_duplicates_OUM<-result_bayou_OUM_old$shifts_per_my[which(result_bayou_OUM_old$Event=="Duplication")]
optima_shifts_per_my_speciation_OUM<- result_bayou_OUM_old$shifts_per_my[which(result_bayou_OUM_old$Event=="Speciation")]
median(optima_shifts_per_my_old_duplicates_OUM) #0.0025
median(optima_shifts_per_my_speciation_OUM) ##0.012
wilcox.test(optima_shifts_per_my_old_duplicates_OUM,optima_shifts_per_my_speciation_OUM) #p-value < 2.2e-16


##################### Optima shift analyses on 3244 OU (OU1+OUM) trees #############

## Analysis considering all branches of the OU trees
result_bayou_OU_new<-result_bayou_OU[result_bayou_OU$pp >= 0.7,]

## Unique OU ID with pp > 0.1
length(unique(result_bayou_OU$index)) ## 2676 out of 3244 (82.49%)

## Unique OU ID with pp >= 0.7
length(unique(result_bayou_OU_new$index)) ## 1101 out of 3244 (33.94%)
Multioptima_Baysian_index<-unique(result_bayou_OU_new$index)
Trees_OUM_Baysian<-Tree_OU[c(Multioptima_Baysian_index)]
rm(Multioptima_Baysian_index)

## To obtain total numbers of events
table_Baysian_OU<-NULL
for(i in 1:length(Trees_OUM_Baysian))
{
  table_Baysian_OU<-bind_rows(table_Baysian_OU,tibble(index=i,
                                                        n_speciation=Trees_OUM_Baysian[[i]]@phylo$n_speciation,
                                                        n_duplication_young=Trees_OUM_Baysian[[i]]@phylo$n_duplication_young,
                                                        n_duplication_old=Trees_OUM_Baysian[[i]]@phylo$n_duplication_old))
  
}
rm(i)
sum(table_Baysian_OU$n_speciation) ## 13824 speciation events
sum(table_Baysian_OU$n_duplication_young) ## 3027 young duplication events
sum(table_Baysian_OU$n_duplication_old) ## 2814 old duplication events

## Regime shift analyses considering all events of a tree
OU_all<-proportion_calculation(Tree_OU,"all",result_bayou_OU_new)

median(na.omit(OU_all$shift_prob_dup)) #0.1
median(na.omit(OU_all$shift_prob_spe)) #0.023
mean(na.omit(OU_all$shift_prob_dup)) #0.128
mean(na.omit(OU_all$shift_prob_spe)) #0.049
pvalue_regime_shift_all_OU<- wilcox.test(OU_all$shift_prob_spe,OU_all$shift_prob_dup,paired = T)$p.value #p-value = 7.56e-55

## Analysis considering young duplication and speciation branches of the OU trees
result_bayou1<-result_bayou_OU_new[result_bayou_OU_new$node_age<=296,]

OU_young<-proportion_calculation(Tree_OU,"young",result_bayou1)

median(na.omit(OU_young$shift_prob_dup)) #0.045
median(na.omit(OU_young$shift_prob_spe)) #0.031
mean(na.omit(OU_young$shift_prob_dup)) #0.147
mean(na.omit(OU_young$shift_prob_spe)) #0.066
pvalue_regime_shift_OU_young<- wilcox.test(OU_young$shift_prob_spe,OU_young$shift_prob_dup,paired = T)$p.value #p-value = 3.43e-12

## Analysis considering old duplication and speciation branches of the OU trees
result_bayou2<-NULL
result_bayou_OU_spe<-result_bayou_OU_new[result_bayou_OU_new$Event=="Speciation",]
result_bayou_OU_dup_old<-result_bayou_OU_new[result_bayou_OU_new$Event=="Duplication" & result_bayou_OU_new$node_age > 296,]
result_bayou2<-rbind(result_bayou_OU_spe,result_bayou_OU_dup_old)

OU_old<-proportion_calculation(Tree_OU,"old",result_bayou2)

median(na.omit(OU_old$shift_prob_dup)) #0.1
median(na.omit(OU_old$shift_prob_spe)) #0.026
mean(na.omit(OU_old$shift_prob_dup)) #0.119
mean(na.omit(OU_old$shift_prob_spe)) #0.039
pvalue_regime_shift_OU_old<- wilcox.test(OU_old$shift_prob_spe,OU_old$shift_prob_dup,paired = T)$p.value #p-value = 7.26e-40

## Calculating the regime shift rates/My

## Step 1: Estimating length of the branch where optimum shift was detected
for(i in 1:nrow(result_bayou_OU_new))
{
  
  tre<-NULL
  index<-NULL
  edgelength<-NULL
  index<-result_bayou_OU_new$index[i]
  tre<-Tree_OU[[index]]@phylo
  result_bayou_OU_new$edgelength[i]<-tre$edge.length[which(tre$edge[,1]==result_bayou_OU_new$node1[i] & tre$edge[,2]==result_bayou_OU_new$node2[i])]
}

## Step 2: Based on the relative location provided by Bayesian approach, we estimated the time when optimum shifts took place 
result_bayou_OU_new$shift_age<-result_bayou_OU_new$node_age-(result_bayou_OU_new$edgelength*result_bayou_OU_new$rel.location)

## Rate calculations
result_bayou_OU_new$shifts_per_my<-1/result_bayou_OU_new$shift_age

## For young duplicates
result_bayou_OU_young<-result_bayou_OU_new[result_bayou_OU_new$node_age <= 296,]
optima_shifts_per_my_young_duplicates<-result_bayou_OU_young$shifts_per_my[which(result_bayou_OU_young$Event=="Duplication")]
optima_shifts_per_my_speciation<- result_bayou_OU_young$shifts_per_my[which(result_bayou_OU_young$Event=="Speciation")]
median(optima_shifts_per_my_young_duplicates) #0.0309
median(optima_shifts_per_my_speciation) ##0.0129
wilcox.test(optima_shifts_per_my_young_duplicates,optima_shifts_per_my_speciation) #p-value = 2.341e-11

## For old duplicates
result_bayou_OU_old<-NULL
result_bayou_OU_old<-result_bayou_OU_new[result_bayou_OU_new$Event=="Speciation",]
result_bayou_OU_old<-rbind(result_bayou_OU_old,result_bayou_OU_new[result_bayou_OU_new$Event=="Duplication" & result_bayou_OU_new$node_age > 296,])
result_bayou_OU_old<-na.omit(result_bayou_OU_old)
optima_shifts_per_my_old_duplicates<-result_bayou_OU_old$shifts_per_my[which(result_bayou_OU_old$Event=="Duplication")]
optima_shifts_per_my_speciation<- result_bayou_OU_old$shifts_per_my[which(result_bayou_OU_old$Event=="Speciation")]
median(optima_shifts_per_my_old_duplicates) #0.0023
median(optima_shifts_per_my_speciation) ##0.0129
wilcox.test(optima_shifts_per_my_old_duplicates,optima_shifts_per_my_speciation) #p-value < 2.2e-16

duplication_shifts_per_my<-result_bayou_OU_new$shifts_per_my[which(result_bayou_OU_new$Event=="Duplication")]
speciation_shifts_per_my<- result_bayou_OU_new$shifts_per_my[which(result_bayou_OU_new$Event=="Speciation")]
median(duplication_shifts_per_my) ## 0.0034
median(speciation_shifts_per_my) ## 0.013
wilcox.test(duplication_shifts_per_my,speciation_shifts_per_my)

#save.image("Model_fitting_TMRR.Rdata")

## To identify trees for which OU1 was the preferred model by Maximum likelihood framework
df<-NULL
for(i in 1:nrow(result_bayou_OU_new))
{
  index<-NULL
  index<-result_bayou_OU_new$index[i]
  
  if(Tree_OU[[index]]@phylo$best_fit_model=="OU1")
  {
    df<-append(df,index)
  }
}
