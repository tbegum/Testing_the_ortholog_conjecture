
############################## Functions written for this study #######################################


##Checking trees with atleast one duplication event
tree_with_duplication<- function(tree)
{
  ## Collecting total internal node data, no. of tips,speciation and duplication nodes
  gene_tree <- tree@phylo
  gene_tree_data <- tree@data
  
  ##Identifying duplication node
  duplication_node <- gene_tree_data$node[which(gene_tree_data$Event=="Duplication")]
  
  ## Now checking for trees with all speciation events
  ## If does not match the criteria returns NA
  if((length(duplication_node)) > 0) 
  {
    return(tree)
  }
  else{
    return(NA)
  }
}

## This function is written to measure phylogenetic signal of the trees
Phylogenetic_signal<-function(tree)
{
  gene_tree<-tree@phylo
  gene_data<-tree@data
  trait_data<-gene_data$Tau[which(is.na(gene_data$pic))]
  names(trait_data)<-gene_data$label[which(is.na(gene_data$pic))]
  Blomberg_K<- tryCatch(phylosig(gene_tree,trait_data,method = "K",test = T), error = function(e) {"Error"})
  if(Blomberg_K=="Error")
  {
    K_Blomberg<-1000
    P_Blomberg<-1000
  }
  if(Blomberg_K!="Error")
  {
    K_Blomberg<-Blomberg_K$K
    P_Blomberg<-Blomberg_K$P
  }
  Pagel_lambda<-tryCatch(phylosig(gene_tree,trait_data,method = "lambda",test = T), error = function(e) {"Error"})
  if(Pagel_lambda=="Error")
  {
    lambda_Pagel<-1000
    P_pagel<-1000
  }
  if(Pagel_lambda!="Error")
  {
    lambda_Pagel<-Pagel_lambda$lambda
    P_pagel<-Pagel_lambda$P
  }
  return(tibble
         (gene=digest(tree),
           K_effecient=K_Blomberg,
           Pvalue_Blomberg=P_Blomberg,
           lambda_estimate=lambda_Pagel,
           Pvalue_pagel=P_pagel))
}

## This function is aimed to summarize the data
summary_function<-function(tree)
{
  if(class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  gene_data<-tree@data 
  gene_data$gene<-digest(tree)
  gene_data$pic<-abs(gene_data$pic)
  tree@data<-gene_data[-c(4,5)]
  return(tree@data)
}

## This function is aimed to summarize the randomized data of BM trees of Dunn et al.(7)
summary_function_BM<-function(tree)
{
  if(class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  gene_data<-tree@data 
  gene_data$gene<-digest(tree)
  gene_data$pic<-abs(gene_data$pic_abs_random)
  tree@data<-gene_data[-c(4,5)]
  return(tree@data)
}

## This function is written to extract speciation contrast
## frame is the dataframe and  x should be mentioned the column name (ex: 'pic' or 'pic_abs_random')
speciation_contrast<-function(frame, x)
{
  contrast_spe <- frame[[x]][which(frame$Event=="Speciation")]
  return(contrast_spe)
}

## This function is written to extract speciation contrast
## frame is the dataframe and  x should be mentioned the column name (ex: 'pic' or 'pic_abs_random')
duplication_contrast<-function(frame, x)
{
  contrast_dup <- frame[[x]][which(frame$Event=="Duplication")]
  return(contrast_dup)
}

## Function to compute Wilcoxon one-tailed test
one_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_one_tailed <- wilcox.test(data1,data2,alternative="greater")$p.value
  
  #star <-stars.pval(wilcox_oc_one_tailed)
  #star0<- stars.pval(0)
  if (wilcox_oc_one_tailed == 0){label_p = paste0("P < 2.2e-16")}
  if ((wilcox_oc_one_tailed != 0))
  {
    wilcox_oc_one_tailed <- format(wilcox_oc_one_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_one_tailed)
  }
  return(label_p)
} 

## Function to compute Wilcoxon two-tailed test
two_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_two_tailed <- wilcox.test(data1,data2,alternative="two.sided")$p.value
  
  #star <-stars.pval(wilcox_oc_two_tailed)
  #star0<- stars.pval(0)
  if (wilcox_oc_two_tailed == 0){label_p = paste0("P < 2.2e-16")}
  if ((wilcox_oc_two_tailed != 0))
  {
    wilcox_oc_two_tailed <- format(wilcox_oc_two_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_two_tailed)
  }
  return(label_p)
} 


### This function allows to perform permutation of tau data as the trait
shuffling_tau<-function(tree)
{
  ## Reading gene tree and its data slot
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting Tau data to randomize
  #Tau.tree<-gene_data$Tau[which(is.na(gene_data$pic))]
  Tau.tree<- gene_data$Tau[ is.tip.nhx( tree ) ]
  
  ## Permuting the trait data and storing the result
  Tau.new<-sample(x=Tau.tree,size = length(Tau.tree), replace = FALSE)
  tree@data$Tau_new<- c(Tau.new,rep(NA, times=gene_tree$Nnode))
  return(tree)
}

### This function is written to perform reshuffle the events of the trees with tau/mean gene expression level as trait data
shuffling_event<-function(tree)
{
  ## Reading gene tree and its data slot
  gene_tree<-tree@phylo
  gene_data<-tree@data
  Internal_node<-gene_data$node[which(!is.na(gene_data$pic))]
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting event data to permute
  event_tree<-gene_data$Event[which(!is.na(gene_data$pic))]
  events_new<-sample(x=event_tree, size = length(event_tree), replace = FALSE )
  Internal_new_event<-as.character(events_new)
  #names(events.new)<-Internal_node
  
  tree@data$event_new<- c(rep(NA,time=length(gene_tree$tip.label)),Internal_new_event)
  tree@data$event_new <- factor(tree@data$event_new, levels=c( "Speciation", "Duplication"))
  return(tree)
}

##This function calculates Phylogenetic indepent contrast for trees based on tau data 
contrast_calc<-function(tree)
{
  ## Tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting tip data 
  Tau_tip<- gene_data$Tau[ is.tip.nhx( tree ) ]
  
  ## Initializing variable
  tree@data$PIC <- NULL
  tree@data$Variance <- NULL
  
  ## Calculating the phylogenetic independent contrasts of tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, gene_tree, var.contrasts=TRUE)
  
  ##Returning the contrast value
  tree@data$PIC<- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,1])
  tree@data$Variance <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,2])
  tree@data$pic<- abs(tree@data$PIC) ## absolute PIC values
  return(tree)
}

##This function calculates Phylogenetic indepent contrast for trees with randomized tau data 
contrast_random<-function(tree)
{
  ## Tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  ## Collecting tip data 
  #Tau_tip<-tree@data$Tau_new[which(is.na(tree@data$pic))]
  Tau_tip<- gene_data$Tau_new[ is.tip.nhx( tree ) ]
  
  ## Initializing variable
  tree@data$PIC_random <- NULL
  tree@data$Variance <- NULL
  
  ## Calculating the phylogenetic independent contrasts of tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, gene_tree, var.contrasts=TRUE)
  
  ##Returning the contrast value
  tree@data$PIC_random <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,1])
  tree@data$Variance <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,2])
  tree@data$pic_abs_random <- abs(tree@data$PIC_random) ## absolute PIC values
  return(tree)
}

## To find out calibrated trees with negative branch lengths
## This function also helps to find out calibrated trees with no duplication events
negative_edgelength<-function(tree) 
{
  index_negative<-vector()
  #species_tree<-vector()
  
  for(i in 1:length(tree))
  {
    calibrated_tree <-tree[[i]]@phylo
    edgelength<-calibrated_tree$edge.length
    negative <- which(edgelength<=0)
    if(length(negative) > 0) 
    {
      index_negative <- append(index_negative,i)
    }
  }
  return(index_negative)
}

## This function helps to identify, and remove trees for which contrast is not properly standarized
## We will obtain contrasts standardized tree of tau
## Crunch method of Caper package will provide the contrast for Tau as well as Mean Expression level
## We just use this function to generate diagnostic plots for tau
## In this function, we use the trait value (tau) without log transformation as it is already a value between 0 to 1
diagnostic_plot_test <- function(tree)
{
  gene_tree<-tree@phylo
  gene_tree$node.label<-NULL
  count<<-count+1
  print(count)
  
  if(class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  data<-tree@data
  data_new<-data[which(!is.na(data$Tau)),] 
  data_new1<-data.frame(label=data_new$label,Tau=data_new$Tau, Mean_Exp=data_new$Mean)
  rownames(data_new1)<-data_new1$label
  treedata<-comparative.data(gene_tree,data_new1,label)
  test<-tryCatch(crunch(Mean_Exp~Tau,data=treedata,equal.branch.length=F), error = function(e) {"Error"})
  if((test=="Error") || (is.na(test$mod[1]$coefficients[1])))
  { print ("Error obtained in crunch")
    return(NA)
  }
  if(test!="Error")
  {
    ## We need to check for diagnostic plots for expected variance and node age from this Caper package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    contrast$height<-data$nodeheight[which(!is.na(data$pic))]
    diagnostic<-caic.diagnostics(test)
    
    ## To added diagnstic test for node age as well
    test_depth<-summary(lm(contrast$Tau~contrast$nodeDepth))
    test_height<-summary(lm(contrast$Tau~contrast$height))
   
    
    ## To remove the trees for which contrast value is not properly standized, 
    ## we need to collect the P values of the diagnostic plot
    p.SDT<-caic.diagnostics(test)[2,4,1]
    p.AgeT<-caic.diagnostics(test)[3,4,1]
    p.depthT<-test_depth$coefficients[2,4]
    p.heightT<-test_height$coefficients[2,4]
    
    ##If contrasts are not properly standized, we return NA
    if((p.SDT < 0.05) || (p.AgeT < 0.05) || (p.depthT < 0.05) || (p.heightT < 0.05)){return (NA) }
    
    ## Collecting the initial pic data in a different variable name
    if((p.SDT >= 0.05) && (p.AgeT >= 0.05) && (p.depthT >= 0.05) && (p.heightT >= 0.05))
    {
      tree@data$pic_Tau<-NULL
      tree@data$variance<-NULL
      
      ## Now we return the trees dataframe
      tree@data$pic_Tau <- c(rep(NA, length(gene_tree$tip.label)), (contrast$Tau))
      tree@data$variance<- c(rep(NA, length(gene_tree$tip.label)), contrast$contrVar)
      return(tree)
    }}
}

## This function helps to identify BM trees by excluding oldest duplication nodes from diagnostic tests
## Nodelimit is used to exclude oldest duplication nodes
## If limit == 0, this refers, we do not exclude any oldest duplication nodes
## Otherwise, we use specific age limit
diagnostic_plot_test_excluding_oldest_dup_nodes <- function(tree,nodelimit)
{
  gene_tree<-tree@phylo
  gene_tree$node.label<-NULL
  count<<-count+1
  print(count)
  
  if(class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  data<-tree@data
  data_new<-data[which(!is.na(data$Tau)),] 
  data_new1<-data.frame(label=data_new$label,Tau=data_new$Tau, Mean_Exp=data_new$Mean)
  rownames(data_new1)<-data_new1$label
  treedata<-comparative.data(gene_tree,data_new1,label)
  test<-tryCatch(crunch(Mean_Exp~Tau,data=treedata,equal.branch.length=F), error = function(e) {"Error"})
  if((test=="Error") || (is.na(test$mod[1]$coefficients[1])))
  { print ("Error obtained in crunch")
    return(NA)
  }
  if(test!="Error")
  {
    ## We need to check for diagnostic plots for expected variance and node age from this Caper package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    if(nodelimit==0){contrast$height<-data$nodeheight[which(data$node_age>0)]}
    if(nodelimit >0){contrast$height<-data$nodeheight[which(!is.na(data$pic))]}
    
    
    ## Considering contrasts data within age limit of 370 My
    if(nodelimit==0){contrast_new<-contrast}
    if(nodelimit>0){contrast_new<-contrast[contrast$nodeAge<=nodelimit,]}

    
    ## length of the dataframe should be higher than 4
    if(nrow(contrast_new) >= 4)
    {
      ## Diagnstic tests for expected variance, node age, node depth and node height as recommended in earlier studies
      test_expvar<-summary(lm(contrast_new$Tau~sqrt(contrast_new$contrVar)))
      test_age<-summary(lm(contrast_new$Tau~log(contrast_new$nodeAge)))
      test_depth<-summary(lm(contrast_new$Tau~contrast_new$nodeDepth))
      test_height<-summary(lm(contrast_new$Tau~contrast_new$height))
    
      if((!(is.na(test_expvar$coefficients[2,4]))) || (!(is.na(test_age$coefficients[2,4]))) || (!(is.na(test_depth$coefficients[2,4]))) || (!(is.na(test_height$coefficients[2,4]))))
      {
        ## To remove the trees for which contrast value is not properly standized, 
        ## we need to collect the P values of the diagnostic plot
        p.SDT<-test_expvar$coefficients[2,4]
        p.AgeT<-test_age$coefficients[2,4]
        p.depthT<-test_depth$coefficients[2,4]
        p.heightT<-test_height$coefficients[2,4]
    
        ##If contrasts are not properly standized, we return NA
        if((p.SDT < 0.05) || (p.AgeT < 0.05) || (p.depthT < 0.05) || (p.heightT < 0.05)){return (NA) }
    
        ## Collecting the initial pic data in a different variable name
        if((p.SDT >= 0.05) && (p.AgeT >= 0.05) && (p.depthT >= 0.05) && (p.heightT >= 0.05))
        {
          tree@data$pic_Tau<-NULL
          tree@data$variance<-NULL
          tree@data$pic<-NULL
      
          ## Now we return the trees dataframe
          tree@data$pic_Tau <- c(rep(NA, length(gene_tree$tip.label)), (contrast$Tau))
          tree@data$variance<- c(rep(NA, length(gene_tree$tip.label)), contrast$contrVar)
          tree@data$pic<- abs(tree@data$pic_Tau)
          return(tree)
        }}}}
}

##This function calculates Phylogenetic indepent contrast for trees with simulated tau data 
contrast_sim<-function(tree)
{
  ## Tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting tip data 
  Tau_tip = tree@data$Tau[is.tip.nhx(tree)]
  
  ## Initializing variable
  tree@data$PIC_sim <- NULL
  tree@data$Var_sim <- NULL
  
  ## Calculating the phylogenetic independent contrasts of tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, gene_tree, var.contrasts=TRUE)
  
  ##Returning the contrast value
  tree@data$PIC_sim  <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,1])
  tree@data$Var_sim <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,2])
  tree@data$pic<- abs(tree@data$PIC_sim) ## absolute PIC values
  return(tree)
}

##This function is written to plot histogram on P values of randomized data
## Frame should be dataframe and xName should be 'pval_tau' or 'pval_event'
##groupName should be 'treeset'
histo_plot<-function(frame, xName,groupName)
{ 
  ggplot2.histogram(data=frame, xName=xName,
  groupName=groupName,legendPosition="none",groupColors = "#7CAE00",
  alpha=0.8, addDensity=F,
  addMeanLine=F,meanLineColor = "#C77CFF", meanLineSize=0.5, meanLineType = "solid") +
  theme_classic() +
  xlab(expression(bold(paste(bolditalic(P)," value"))))+ 
  ylab("Frequency") +
  theme(legend.title=element_blank(),legend.position=c(0.7,0.8),legend.text=element_text(size=8)) +
  theme(legend.position="none")+
  theme(plot.title = element_text (face="bold", 
                                   colour="black", lineheight=1.0),
        axis.title.x=element_text( face="bold",
                                   colour="black", hjust=0.5),
        axis.title.y=element_text(face="bold",
                                  colour="black", vjust=0.5, angle=90))
}


##This function is written to collect median data required for jitter plot
## Frame should be dataframe,x and y should be defined (for example pic and Event) 
median_jitter<-function(frame,x,y)
{
  median_data<-NULL
  median_data<- aggregate(frame[[x]]~frame[[y]], frame, median) ## aggregate function format
  median_data$pic.round<-paste0("Median = ",round(median_data$`frame[[x]]`,4),sep="")
  return(median_data)
} 


##This function is written to generate jitter plot
## Frame should be dataframe
## df_median is the dataframe with median data
jitter_plot_tau<-function(frame,df_median,pvalue)
{
  ggplot(frame,aes(x=Event, y=pic_Tau, fill=Event)) + 
    geom_jitter(aes(colour = Event),position=position_jitter(0.02), alpha=0.5) +
    geom_text(data = df_median, aes(label = pic.round),fontface = 2,size = 3, 
              hjust=0.5,vjust =-2.5)+
    geom_crossbar(data=df_median, aes(ymin = pic_Tau, ymax = pic_Tau),
                  size=0.2,col= "black", width = .2)+
    xlab( NULL ) +
    ylab(expression(bold("PIC"))) +
    ylim(0, 0.2) +
    theme_classic()+
    theme(legend.title=element_blank(),legend.position=c(1.9,1.9)) +
    theme(axis.text=element_text(size=10,face="bold")) +
    theme(legend.text=element_text(size=10,face="bold")) +
    annotate("text", x = 1.5, y = 0.17, label= pvalue, fontface = 4)
}

## This function is written to calculate the proportion of duplication, and speciation events
tree_data_collection<-function(tree)
{
  ##Collecting gene tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  ntips <-length(gene_tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1)
  root_event<-gene_data$Event[root]
  root_age<-gene_data$node_age[root]
  Internal_event_number<-as.numeric(length(gene_data$node[which(!is.na(gene_data$pic))]))
  
  ## initializing variables
  pdup<-0
  pspe<-0
  pNA<-0
  
  ##Since, many internal node events are assigned as "NA", sum of proportion of speciation and duplication may not be equal to 1
  dup<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication"))])) 
  pdup<-round(dup/Internal_event_number,2) ## proportion of duplication
  spe<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Speciation"))])) 
  pspe<-round(spe/Internal_event_number,2) ## proportion of speciation
  NA_event<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (is.na(gene_data$Event)))])) 
  pNA<-round(NA_event/Internal_event_number,2) ## proportion of speciation
  spe_var<-round(median(gene_data$var_exp[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Speciation"))]),0) ## variance of speciation
  dup_var<-round(median(gene_data$var_exp[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication"))]),0) ## variance of speciation
  if(pdup==0){dup_var<-NA}
  return(tibble(tip_num=ntips,
                internal_events=Internal_event_number,
                dup_num=dup,
                dup_prop=pdup,
                spe_prop=pspe,
                NA_prop=pNA,
                dup_var=dup_var,
                spe_var=spe_var,
                root_age=root_age,
                root_event=root_event))
  
}

## This function additionally tests the Brownian motion of contrasts distribution of tau
post_diagnostic_test_analysis <- function(tree)
{
  count<<-count+1
  print(count)
  #Collecting data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  #correlation analyis to check for significant relation
  tau_contrast<-abs(gene_data$pic_Tau)[which(is.na(gene_data$Tau))]
  standard_deviation<-sqrt(gene_data$var_exp[which(is.na(gene_data$Tau))])
  age_node<-gene_data$node_age[which(is.na(gene_data$Tau))]
  depth_node<-gene_data$node_depth[which(is.na(gene_data$Tau))]
  
  ## If 4 or more observations are available, we considered those trees 
  if(length(na.omit(tau_contrast))>=4)
  {
    #correlation
    pval_abscontrast_stdev<-cor.test(tau_contrast,standard_deviation)$p.value
    pval_abscontrast_age<-cor.test(tau_contrast,age_node)$p.value
    pval_abscontrast_depth<-cor.test(tau_contrast, depth_node)$p.value
    
    if(pval_abscontrast_stdev >= 0.05 && pval_abscontrast_age >= 0.05 && pval_abscontrast_depth >= 0.05)
    {
      
      # retuning values
      return(tibble(gene=count,
                    pval_abscontrast_stdev=pval_abscontrast_stdev,
                    pval_abscontrast_age=pval_abscontrast_age,
                    pval_abscontrast_depth=pval_abscontrast_depth))
    }  
  }
}
  
## This function adds node height to our dataframe of tree@data for further use in building time calibration matrix
tree_height <- function(tree)
{
  gene_tree <- tree@phylo
  gene_data <- tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Identifying nodes
  gene_tree_nodes<-gene_data$node
  height_node<-NULL
  
  ##Computing node height and returned it to the tree into "@data" slot
  height_node<-ape::node.height(gene_tree)
  tree@data$nodeheight<-height_node
  return(tree)
}

## This function is written for null and OC simulation of empirical trees
simulation<-function(tree,dup_adjust)
{

  ## Considering real tree data
  gene_tree<-tree@phylo
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  gene_data<-tree@data
  
  ## Identifiying duplication nodes and edges to paint them to simulate them at specified trait evolutionary rate
  dup_nodes <- gene_data$node[which(gene_data$Event=="Duplication")]
  dup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% dup_nodes), 2])
  gene_tree$edge.length[dup_edges] = gene_tree$edge.length[dup_edges]*dup_adjust
  tree_painted <- paintBranches (gene_tree, edge=dup_edges, "D", anc.state="S")
  
  ## Values
  anc<-gene_tree$model_bm$opt$z0
  sig2<-gene_tree$model_bm$opt$sigsq
  
  ## Declaring vector for trait simulation
  #sigma_vector <- c(sig2, sig2*dup_adjust) ##Trait evolution rate differs when we consider OC
  #names(sigma_vector) <- c("S", "D")
  
  ## Trait simulation
  
  
  trait_simulated<-fastBM(tree_painted, a=anc, sig2,bounds=c(0,1))
  #trait_simulated<-sim.rates_bound(tree_painted,sigma_vector,anc=anc)
  
  ## Replacing the original tau data with simulated one 
  ## and return data to the tree
  ## Initializing variable
  tree@data$pic_null_sim <- NULL
  tree@data$pic_null_sim_abs <- NULL
  
  ## Calculating the phylogenetic independent contrasts of tree
  ## Returning the results to "data" frame of the tree 
  pic_null <- ape::pic(trait_simulated, gene_tree, var.contrasts=TRUE)
  
  ##Returning the contrast value
  tree@data$pic_null_sim<- c(rep(NA, length(gene_tree$tip.label)), pic_null[,1])
  tree@data$pic_null_sim_abs<- abs(tree@data$pic_null_sim) ## absolute PIC values
  return(tree)
}

## This function is written for null brownian model simulation of empirical trees
simulation.new<-function(tree,dup_adjust)
{
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  ## Considering real tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  tau <- gene_data$Tau[which(is.tip.nhx(tree))]
  names(tau) <- gene_tree$tip.label
  
  ## Identifiying duplication nodes and edges to paint them to simulate them at specified trait evolutionary rate
  dup_nodes <- gene_data$node[which(gene_data$Event=="Duplication")]
  dup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% dup_nodes), 2])
  tree_painted <- paintBranches (gene_tree, edge=dup_edges, "D", anc.state="S")
  
  ##
  fit.bm<-fitContinuous(tree_painted, tau, model = "BM")
  anc<-fit.bm$opt$z0
  sig2<-fit.bm$opt$sigsq
  
  ## Declaring vector for trait simulation
  sigma_vector <- c(sig2, sig2*dup_adjust) ##Trait evolution rate differs when we consider OC
  names(sigma_vector) <- c("S", "D")
  
  ## Trait simulation
  trait_simulated<-sim.rates_bound(tree_painted, sig2=sigma_vector, anc=anc)
  
  ## Replacing the original tau data with simulated one 
  ## and return data to the tree
  tree@data$Tau<-c(trait_simulated,rep(NA, time=gene_tree$Nnode))
  return(tree)
}

## The following function was taken from easyGgplot2 package and modified to add median trendline in the histogram plot

ggplot2.histogram.mod<-function (data, xName = NULL, groupName = NULL, position = c("identity", 
                                                             "stack", "dodge"), addMedianLine = FALSE, medianLineColor = NULL, 
          medianLineType = "dashed", medianLineSize = 1, addDensityCurve = FALSE, 
          densityFill = "#FF6666", densityAlpha = 0.2, densityLineType = "solid", 
          densityLineColor = "#2F2F2F", scale = c("frequency", "density"), 
          groupColors = NULL, brewerPalette = NULL, fill = "black", 
          color = "black", linetype = "solid", ...) 
{
  .hist_params <- function(...){
    x <- list(...)
    res <- list()
    res$binwidth <- x$binwidth
    res$bins <- x$bins
    return(res)
  }
  pms <- .hist_params(...)
  alpha <- ifelse(is.null(groupName), 1, 0.5)
  if (is.null(xName) & !is.numeric(data)) 
    stop("xName is missing or NULL. In this case data should be a numeric vector")
  else if (is.numeric(data)) {
    data = cbind(x = data, grp = rep(1, length(data)))
    xName = "x"
  }
  data = data.frame(data)
  if (is.null(groupName)) 
    p <- ggplot(data = data, aes_string(x = xName))
  else {
    data[, groupName] = factor(data[, groupName])
    p <- ggplot(data = data, aes_string(x = xName, fill = groupName, 
                                        colour = groupName))
  }
  if (addDensityCurve) {
    if (!is.null(groupName)) {
      p <- ggplot(data = data, aes_string(x = xName, fill = groupName, 
                                          colour = groupName))
      p <- p + geom_histogram(aes_string(y = "..density.."), 
                              position = position[1], binwidth = pms$binwidth, 
                              bins = pms$bins, alpha = alpha)
      p <- p + geom_density(alpha = densityAlpha, linetype = densityLineType)
    }
    else {
      p <- p + geom_histogram(aes_string(y = "..density.."), 
                              position = position[1], binwidth = pms$binwidth, 
                              bins = pms$bins, color = color, fill = fill, 
                              linetype = linetype)
      p <- p + geom_density(fill = densityFill, alpha = densityAlpha, 
                            linetype = densityLineType, colour = densityLineColor)
    }
  }
  else {
    if (scale[1] == "density") 
      p <- p + geom_histogram(aes_string(y = "..density.."), 
                              position = position[1], binwidth = pms$binwidth, 
                              bins = pms$bins)
    else p <- p + geom_histogram(aes_string(y = "..count.."), 
                                 position = position[1], binwidth = pms$binwidth, 
                                 bins = pms$bins, alpha = alpha)
  }
  if (addMedianLine) {
    if (is.null(groupName)) {
      if (is.null(meanLineColor)) 
        medianLineColor = "red"
      m = median(data[, xName], na.rm = T)
      p <- p + geom_vline(aes_string(xintercept = m), color = medianLineColor, 
                          linetype = medianLineType, size = medianLineSize)
    }
    else {
      df <- data.frame(grp = factor(data[, groupName]), 
                       x = data[, xName])
      df.m <- stats::aggregate(df[, "x"], by = list(grp = df[, 
                                                             "grp"]), median)
      names(df.m) <- c(groupName, "x.median")
      if (is.null(medianLineColor)) 
        p <- p + geom_vline(data = df.m, aes_string(xintercept = "x.median", 
                                                    colour = groupName), linetype = medianLineType, 
                            size = medianLineSize)
      else p <- p + geom_vline(data = df.m, aes_string(xintercept = "x.median"), 
                               linetype = medianLineType, color = medianLineColor, 
                               size = medianLineSize)
    }
  }
  if (!is.null(groupColors)) {
    p <- p + scale_fill_manual(values = groupColors)
    p <- p + scale_colour_manual(values = groupColors)
  }
  else if (!is.null(brewerPalette)) {
    p <- p + scale_fill_brewer(palette = brewerPalette)
    p <- p + scale_colour_brewer(palette = brewerPalette, 
                                 guide = "none")
  }
  p <- ggplot2.customize(p, ...)
  p
}

