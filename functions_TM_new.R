
############################## Functions written for this study #######################################


##Selecting gene trees with atleast one duplication event
tree_with_duplication<- function(tree)
{
  ## Obtaining total internal node data, number of tips,speciation and duplication nodes for each tree
  gene_tree <- tree@phylo
  gene_tree_data <- tree@data
  
  ##Identifying duplication nodes
  duplication_node <- gene_tree_data$node[which(gene_tree_data$Event=="Duplication")]
  
  ## If does not match the criteria returns NA
  if((length(duplication_node)) > 0) 
  {
    return(tree)
  }
  else{
    return(NA)
  }
}

## This function is written to estimate phylogenetic signals (with Blomberg's K and with Pagel's lambda methos) for each gene tree
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

## This function is written to modify the labels of duplication nodes
## Otherwise duplication nodes with the same label name as of speciation nodes will be calibrated using the speciation time points or will interfere during  the time calibration process  
modify_label <- function(tree)
{
  gene_tree <- tree@phylo
  gene_tree_data <- tree@data
  
  ## Collecting internal clade labels of duplication nodes to modify them
  Internal_clade_label <- gene_tree_data$label[which(!is.tip.nhx(tree))]
  dup_node<-gene_tree_data$node[which(gene_tree_data$Event=="Duplication")]
  dup_node_data_label<-gene_tree_data$label[dup_node]
  tree@data$label[dup_node] <- paste(dup_node_data_label,"_d",sep = "")
  return(tree)
}

## This function is to time calibrate gene trees on the basis of the speciation time points
## Modified from Dunn's script
## Maintaining the original topology of the trees
tree_calibrate <- function(tree, timeframe, model=model)
{ 
  
  gene_tree <- tree@phylo
  count<<-count+1
  print(count)
  
    ## Create calibration matrix for speciation nodes
     calibration_matrix <- 
      tree@data[ !is.tip.nhx( tree ), ] %>%
      filter( D == "N" ) %>%   
      left_join( timeframe, c( "label" = "clade" ) ) %>%
      mutate( age.min = age ) %>%
      mutate( age.max = age ) %>% 
      mutate( soft.bounds = NA )
    
    calibration_matrix <- calibration_matrix[c("node", "age.min","age.max","soft.bounds")]
    
    ## Time calibrating trees
    calibrate_trees <- try(ape::chronos(gene_tree, lambda = 0, calibration = calibration_matrix, model = model))
    
    ##Trees those are not time calibrated can not be used further
    ##To avoid error due to non calibrated tree we did the following
    if( "phylo" %in% class(calibrate_trees))
    {
      class(calibrate_trees) <- "phylo"
      tree@phylo <- calibrate_trees
      return(tree)
    }
    
    else{
      return(NA)
    }
}

## This function adds node depth to our dataframe of tree@data for further use in building time calibration matrix
tree_nodedepth <- function(tree)
{
  gene_tree <- tree@phylo
  
  ##Computing node depth and returning it to the tree into "@data" slot
  nodedepth <- ape::node.depth(gene_tree)
  tree@data$node_depth<-nodedepth
  return(tree)
}

## This function adds heights of each node to the '@data' slot
tree_height <- function(tree)
{
  gene_tree <- tree@phylo
  gene_data <- tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Identifying nodes of a gene tree
  gene_tree_nodes<-gene_data$node
  height_node<-NULL
  
  ## The purpose of using node height in this study is to use it for one diagnostic test analysis, although using this node height estimate of the ape R package is not very essential
  ## It is recommended to use the nodeheight () of the phytools R package
  ## Even for the diagnostic tests, one can use the nodeheight() of the phytools R package, i.e. the root to the node distance in each tree
  ## We did not consider nodeheight() of the phytools R package in this case because we already included the node age, i.e., the tip to the node distance, in our diagnostic tests. This means that it will not add anything extra in our diagnostic tests.
  ## Hence, we added node.height() of the ape R package. If node age is already included, node height can be excluded from diagnostic tests if some one wants. 
  ## The aim for us is to exclude all possible trees with phylogenetic dependence of any internal node parameters, hence we also kept node height but from the ape R package, not from the phytools R package. In any case (with the node height option of the phytools/ape R package or excluding the  node height), the trends of results will not change
    
  ##Computing node height and returned it to the tree into "@data" slot
  height_node<-ape::node.height(gene_tree)
  tree@data$nodeheight<-height_node
  return(tree)
}

## This function is written to change the modified labels of duplication nodes to the original one 
remodify_label <- function(tree)
{
  gene_tree <- tree@phylo
  gene_tree_data <- tree@data
  
  
  ## Collecting internal clade labels of duplication nodes to modify them
  Internal_clade_label <- gene_tree_data$label[which(!is.tip.nhx(tree))]
  Internal_node<-gene_tree_data$node[which(!is.tip.nhx(tree))]
  New_internal_clade_label<-sapply(Internal_clade_label,
                                   function(x) 
                                   {x <- unlist(strsplit(toString(x), split='_d', fixed=TRUE))[1]})
  tree@data$label[Internal_node] <-as.character(New_internal_clade_label)
  #tree@data$S[dup_node] <- paste(dup_node_data_S,".d",sep = "")
  return(tree)
}

##This function is written to permute Tau at the tips of the trees 
shuffling_tauR<-function(tree)
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
  tree@data$Tau<- c(Tau.new,rep(NA, times=gene_tree$Nnode))
  return(tree)
}

## This function is writting to modify Ensembl protein id to gene id
modify_tiplabel<-function(tree)
{
  ## Reading data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  ## Modifying tip labels
  tree@phylo$tip.label<-tree@data$G[is.tip.nhx(tree)]
  
  ## Modifying tree data labels
  tree_labels<-tree@data$label
  internal_labels<-tree@data$label[!is.tip.nhx(tree)]
  tree@data$label<-c(tree@phylo$tip.label,internal_labels)
  
  return(tree)
}
## This function is to summarize the data of all trees
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

## This function allows to perform permutation of the trait data (tau here)
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
  Tau.tree<- gene_data$Tau[ is.tip.nhx(tree)]
  
  ## Permuting the trait data, and returning the result
  Tau.new<-sample(x=Tau.tree,size = length(Tau.tree), replace = FALSE)
  tree@data$Tau_new<- c(Tau.new,rep(NA, times=gene_tree$Nnode))
  return(tree)
}

## This function is written to perform permutation of the internl node events of gene trees with any trait data
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

##This function calculates Phylogenetic Independent Contrasts (PICs) for the gene trees for any trait data (tau here) 
contrast_calc<-function(tree)
{
  ## Collecting data for each tree
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  if (class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  
  ## Collecting trait data at the tips of each tree  
  Tau_tip<- gene_data$Tau[ is.tip.nhx( tree ) ]
  
  ## Initializing variable
  tree@data$PIC <- NULL
  tree@data$Variance <- NULL
  
  ## Calculating the PICs of the corresponding gene tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, gene_tree, var.contrasts=TRUE)
  tree@data$PIC<- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,1])
  tree@data$Variance <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,2])
  tree@data$pic<- abs(tree@data$PIC) ## absolute PIC values
  return(tree)
}

##This function calculates PICs for each tree with randomized trait data (Tau here)
contrast_random<-function(tree)
{
  ## Tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  
  ## Collecting trait data at the tips of each tree 
  Tau_tip<- gene_data$Tau_new[is.tip.nhx(tree)]
  
  ## Initializing variable
  tree@data$PIC_random <- NULL
  tree@data$Variance <- NULL
  
  ## Calculating the phylogenetic independent contrasts for each tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, gene_tree, var.contrasts=TRUE)
  tree@data$PIC_random <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,1])
  tree@data$Variance <- c(rep(NA, length(gene_tree$tip.label)), pic_tree[,2])
  tree@data$pic_abs_random <- abs(tree@data$PIC_random) ## absolute PIC values
  return(tree)
}

## Identification of calibrated time trees with negative branch lengths
negative_edgelength<-function(tree) 
{
  ## Initialization
  index_negative<-vector()
  
  ## Returning the indices of the lists of gene trees with negative edge lengths
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

## This function helps to identify, and to exclude gene trees for which contrast is not properly standarized
## crunch () of the Caper R package is used for this purpose
## We used two traits here (Tau and Mean expression level) to perform the phylogenetic regression for each tree
## The output for each tree also provides the contrasts for both the traits
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
    ## We need to check for diagnostic plots for the expected variance and node age with the absolute PICs from the Caper R package
    contrast<-caic.table(test)
    contrast$node<-as.numeric(rownames(contrast))
    contrast$height<-data$nodeheight[which(!is.na(data$pic))]
    diagnostic<-caic.diagnostics(test)
    
    ## Adding the node depth and the node height to the diagnstic tests
    test_depth<-summary(lm(contrast$Tau~contrast$nodeDepth))
    test_height<-summary(lm(contrast$Tau~contrast$height))
   
    
    ## We aim to remove the trees for which PICs show any phylogenetic dependence, assessed by the significant correlation(s) between the absolute PICs and either of the four internal parameters considered in this study 
    ## Collecting the P values of the diagnostic tests
    p.SDT<-caic.diagnostics(test)[2,4,1]
    p.AgeT<-caic.diagnostics(test)[3,4,1]
    p.depthT<-test_depth$coefficients[2,4]
    p.heightT<-test_height$coefficients[2,4]
    
    ##If PICs show a significant correlation with any of the four internal phylogenetic parameters, we return 'NA'
    if((p.SDT < 0.05) || (p.AgeT < 0.05) || (p.depthT < 0.05) || (p.heightT < 0.05)){return (NA) }
    
    ## Collecting the PICs in a different variable name
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

## This function is written to transform the branch lengths of the gene trees
branch_transform<-function(tree)
{
  ct<<-ct+1
  print(ct)
  ##Collecting gene tree and their edge lengths
  gene_tree<-tree@phylo
  gene_data<-tree@data
  gene_tree$edge.length[which(gene_tree$edge.length <= 1)] <- 1 
  edge_length<-gene_tree$edge.length
  ntips <-length(gene_tree$tip.label)
  scale<-0 ##Initialization for each tree
  
  if(class(tree) == "treedata") 
  {
    tree@data <- tidytree::as_tibble(tree@data)
  }
  ## Initializing variables
  tree@data$pic_transformed<-NULL
  tree@data$var_transformed<-NULL
  
  ## Transforming the branch length
  tree@phylo$edge.length.new<-log(edge_length,10)+edge_length**scale
  
  ##Computing the PICs for the branch length transformed trees
  transformed_tree<-contrast_br_transformed(tree)
  
  if(ntips >= 3)
  {
    ## Performing diagnostic test to check for the non-correlation between standard deviations and absolute PICs for each tree
    ## If there is a significant correlation, we need to continue the branch length transformation for that tree
    ## Calculating the correlation between absolute PICs and standard deviations
    new_tree<-transformed_tree@phylo 
    new_data<-transformed_tree@data
    abs_pic<-new_data$pic_transformed[which(!is.na(new_data$pic))]
    var<-new_data$var_transformed[which(!is.na(new_data$pic))]
    varr<-sqrt(var)
    group<-new_data$Event[which(!is.na(new_data$pic))]
    
    ## When all the node contrasts are zero, we return the tree 
    length_pic<-length(abs_pic)
    length_pic_zero <-length(abs_pic[abs_pic==0])
    if(length_pic == length_pic_zero)
    {
      return(transformed_tree)
    }
    else
    {
      #p.SD<-cor.test(abs_pic,varr)$p.value
      var.sd<-summary(lm(abs_pic~varr))
      p.SD<-var.sd$coefficients[2,4]
      
      
      ## If no significant correlation is found, the corresponding tree is returned
      if(p.SD > 0.05) 
      {
        print(paste0("passed:",scale))
        return(transformed_tree)
      }
      
      ## Else, we perform recursive transformtion
      if(p.SD <= 0.05) 
      {
        if(scale < 2)
        {
          scale<-scale + 0.1
          print(paste0("scale:",scale))
          recursive_transformation(tree,scale)
        }
        else
        {
          return(NA)
        }
      }
    }
  }
  else
  {
    return(NA)
  }
}

#This function calculates PICs for the branch length transformed trees
contrast_br_transformed<-function(tree1)
{
  
  ## Tree data
  gtree<-tree1@phylo
  gdata<-tree1@data
  tree1@phylo$edge.length<-tree1@phylo$edge.length.new
  
  if(class(tree1) == "treedata") 
  {
    tree1@data <- tidytree::as_tibble(tree1@data)
  }
  
  ## Collecting trait data at the tips 
  Tau_tip<-tree1@data$Tau[which(is.na(tree1@data$pic))]

  ## Calculating the PICs for each tree
  ## Returning the results to "data" frame of the tree 
  pic_tree <- ape::pic(Tau_tip, tree1@phylo, var.contrasts=TRUE)
  
  absolute_pic <- abs(pic_tree[,1]) ## absolute PIC values
  
  ##Returning PICs for the branch transformed trees
  tree1@data$pic_transformed <-c(rep(NA, length(gtree$tip.label)), absolute_pic)
  tree1@data$var_transformed <- c(rep(NA, length(gtree$tip.label)), pic_tree[,2])
  
  return(tree1)
}

##This function is to perform recursive branch length transformations
recursive_transformation <- function(phy,scalelimit)
{
  ##Collecting gene trees and their edge lengths
  data_tree<-phy@phylo
  data_edge_length<-data_tree$edge.length
  data_edge_length[which(data_edge_length <= 1)] <- 1
  
  ## Transforming the branch lengths
  phy@phylo$edge.length.new<-log(data_edge_length,10)+data_edge_length**scalelimit

  ## Collecting the PICs in a different variable name
  phy@data$pic_transformed<-NULL
  phy@data$var_transformed<-NULL
  
  ##Computing PICs for the branch length transformed tree
  transformed_tree_new<-contrast_br_transformed(phy)
  
  ## Performing diagnostic test to check for a non-significant correlation between standard deviations and absolute PICs for each tree
  ## Calculating the correlation between absolute PICs and standard deviations
  new_tree<-transformed_tree_new@phylo 
  new_data<-transformed_tree_new@data
  pic<-new_data$pic_transformed[which(!is.na(new_data$pic))]
  var_sqrt<-sqrt(new_data$var_transformed[which(!is.na(new_data$pic))])
  group.new<-new_data$Event[which(!is.na(new_data$pic))]
  var.sd.new<-summary(lm(pic~var_sqrt))
  p.sd<-var.sd.new$coefficients[2,4]
  
  ## If no significant correlation is found, the tree is returned
  if(p.sd > 0.05) 
  {
    print(paste0("passed recur:",scalelimit))
    return(transformed_tree_new)
  }
  
  ## If the P value is still significant, we perform recursive transformtion
  if(p.sd <= 0.05) 
  {
    if(scalelimit <= 2)
    {
      print(paste0("scale:",scalelimit))
      scalelimit<-scalelimit + 0.1
      recursive_transformation(phy,scalelimit)
    }
  }
}

## This function is written to extract phylogenetic indepedent contrasts of the speciation events for the empirical or for the randomized trees 
## frame refers to the dataframe and x should is the column name (ex: 'pic' or 'pic_abs_random' for the empirical or for therandomized data))
speciation_contrast<-function(frame, x)
{
  contrast_spe <- frame[[x]][which(frame$Event=="Speciation")]
  return(contrast_spe)
}

##This function is written to extract phylogenetic indepedent contrasts of the duplication events for the empirical or for the randomized trees
## frame refers to the dataframe and x should is the column name (ex: 'pic' or 'pic_abs_random' for the empirical or for the randomized data))
duplication_contrast<-function(frame, x)
{
  contrast_dup <- frame[[x]][which(frame$Event=="Duplication")]
  return(contrast_dup)
}

##This function is written to extract phylogenetic indepedent contrasts of the young duplication events for the empirical or for the randomized trees
## frame refers to the dataframe and x should is the column name (ex: 'pic' or 'pic_abs_random' for the empirical or for the randomized data))
young_duplication_contrast<-function(frame, x)
{
  contrast_dup <- frame[[x]][which(frame$Event=="Duplication" & frame$node_age <= 296)]
  return(contrast_dup)
}

##This function is written to extract phylogenetic indepedent contrasts of old duplication events for the empirical or for the randomized trees
## frame refers to the dataframe and x should is the column name (ex: 'pic' or 'pic_abs_random' for the empirical or for the randomized data))
old_duplication_contrast<-function(frame, x)
{
  contrast_dup <- frame[[x]][which(frame$Event=="Duplication" & frame$node_age > 296)]
  return(contrast_dup)
}

## Function to compute the Wilcoxon one-tailed test
one_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_one_tailed <- wilcox.test(data1,data2,alternative="greater")$p.value
  
  #star <-stars.pval(wilcox_oc_one_tailed)
  #star0<- stars.pval(0)

  if((wilcox_oc_one_tailed != 0) & (wilcox_oc_one_tailed > 2.2e-16))
  {
    wilcox_oc_one_tailed <- format(wilcox_oc_one_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_one_tailed)
  }
  
  else{label_p = paste0("P < 2.2e-16")}

  return(label_p)
} 

## Function to compute the Wilcoxon two-tailed test
two_tailed_wilcox<-function (data1, data2)
{
  wilcox_oc_two_tailed <- wilcox.test(data1,data2,alternative="two.sided")$p.value
  
  #star <-stars.pval(wilcox_oc_two_tailed)
  #star0<- stars.pval(0)
  if((wilcox_oc_two_tailed != 0) & (wilcox_oc_two_tailed > 2.2e-16))
  {
    wilcox_oc_two_tailed <- format(wilcox_oc_two_tailed, digits= 3, scientific = TRUE)
    label_p =paste0("P = ",wilcox_oc_two_tailed)
  }
  else{label_p = paste0("P < 2.2e-16")}
  return(label_p)
} 

## Function to perform two-sided Wilcoxon test on the branch transformed trees
Wilcoxon_2_sided_transformed<-function(dataframe)
{
  speciation_contrast_tau <- abs(speciation_contrast(dataframe,"pic"))
  duplication_contrast_tau <-  abs(duplication_contrast(dataframe,"pic"))

  ## Performing test 
  wilcox_test_tau_br_transformed <- two_tailed_wilcox(duplication_contrast_tau,speciation_contrast_tau)
  return(wilcox_test_tau_br_transformed)
}

##This function is written to generate boxplots
boxplot_new<-function(dataframe,pval,med,type)
{
  dodge <- position_dodge(width = 0.51)
  
  if(type=="type1") ## "type1" for the data without branch length transformation
  {
    plot<-ggplot(dataframe,aes(x=group, y=pic, fill=Event ) ) + 
      guides( colour = guide_legend( override.aes = list( shape = 16 ) ) ) +
      geom_boxplot( width=0.5,outlier.colour=NA, position = dodge, notch = T) +
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ",tau)))) +
      coord_cartesian(ylim=c(0, 0.05)) +
      theme_classic()+
      theme(legend.title=element_blank(),legend.position=c(0.9,0.9)) +
      theme(axis.text=element_text(size=10,face="bold")) +
      theme(legend.text=element_text(size=10,face="bold"))+
      theme(plot.title = element_text(face = "bold"))+
      annotate("text", x = 1.00, y = 0.04, label= pval, fontface = 4)
     # annotate("text", x=0.87, y=0.01, label=med$pic.round[1], fontface=2)+
    #  annotate("text", x=1.13, y=0.01, label=med$pic.round[2], fontface=2)
    return(plot)
  }
  
  if(type=="type2") ## "type2" for the  branch length transformed data
  {
    plot<-ggplot(dataframe,aes(x=group, y=pic, fill=Event ) ) + 
      guides( colour = guide_legend( override.aes = list( shape = 16 ) ) ) +
      geom_boxplot( width=0.5,outlier.colour=NA, position = dodge, notch = T) +
      xlab( NULL ) +
      ylab(expression(bold(paste("PICs of ",tau)))) +
      coord_cartesian(ylim=c(0, 0.22)) +
      theme_classic()+
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank())+
      theme(legend.title=element_blank(),legend.position=c(0.85,0.7)) +
      #theme(legend.text=element_text(size=10,face="bold"))+
      theme(plot.title = element_text(face = "bold"))+ 
      annotate("text", x = 1.00, y = 0.2, label= pval, fontface = 4)
    return(plot)
  }
}

##This function is written to plot histogram on the P values of the randomized data
## Frame should be the dataframe and xName should be 'pval_tau' or 'pval_event'
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

## The following function is taken from easyGgplot2 R package (https://github.com/kassambara/easyGgplot2)
## and modified to add median trendline in the histogram plot instead of adding meanline
ggplot2_histogram_mod<-function (data, xName = NULL, groupName = NULL, position = c("identity", 
                                                                                    "stack", "dodge"), addMedianLine = FALSE, medianLineColor = NULL, 
                                 medianLineType = "dashed", medianLineSize = 1, addDensityCurve = FALSE, 
                                 densityFill = "#FF6666", densityAlpha = 0.2, densityLineType = "solid", 
                                 densityLineColor = "#2F2F2F", scale = c("frequency", "density"), 
                                 groupColors = NULL, brewerPalette = NULL, fill = "black", 
                                 color = "black", linetype = "solid", ...) 
{
  ## The following function .hish_params() taken from easyGgplot2 package
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
      m = median(data[, xName], na.rm = T) ## Getting median value
      p <- p + geom_vline(aes_string(xintercept = m), color = medianLineColor, 
                          linetype = medianLineType, size = medianLineSize)
    }
    else {
      df <- data.frame(grp = factor(data[, groupName]), 
                       x = data[, xName])
      df.m <- stats::aggregate(df[, "x"], by = list(grp = df[, 
                                                             "grp"]), median) ## Adding median value to different groups
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

## This function generates a dataframe with median data
median_df<-function(Frame)
{
  median_data<-Frame %>% 
    group_by(Event) %>% 
    summarise(Median=median(pic),
              count=n())
  median_data$pic.round<-paste0(round(median_data$Median,4))
  return(median_data)
}

## This function is written to generate dataframe for plotting
plot_frame<-function(dataframe,estimate,type)
{
  plot_df <- data.frame(pic=NA,group=NA, Event=NA) ## declaring data frame
  
  if(type=="type1" && estimate=="pic_Tau")
  {
    plot_df <- rbind(plot_df, data.frame(pic=abs(dataframe$pic_Tau), group="Age <= 296 My", Event=dataframe$Event))
  }
  if(type=="type2" && estimate=="pic")
  {
    plot_df <- rbind(plot_df, data.frame(pic=abs(dataframe$pic), group="No age limit", Event=dataframe$Event))
  }
  
  plot_df <- plot_df[-1,]
  plot_df$Event <- factor(plot_df$Event, levels=c("Speciation", "Duplication"))
  return(plot_df)
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

## This function is written to calculate the proportions of duplication, and speciation events
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

## This function is written to calculate the proportions of duplication (according to age), speciation and NA events after the model fitting
tree_data_statistics<-function(tree, duplication_type)
{
  ##Reading gene tree data slot
  gene_tree<-tree@phylo
  gene_data<-tree@data
  ntips <-length(gene_tree$tip.label) ## Number of tips
  root<-as.numeric(ntips+1) 
  root_event<-gene_data$Event[root]
  root_age<-gene_data$node_age[root]
  Internal_event_number<-as.numeric(length(gene_data$node[which(!is.na(gene_data$pic))])) ##Length of internal nodes excluding nodes assigned to tips
  
  ## Initializing variables
  pdup<-0
  pspe<-0
  pNA<-0
  
  ##Since many internal node events are assigned as "NA", sum of proportion of speciation and duplication may not be equal to 1
  if(duplication_type=="young")
  { 
    dup<-NULL
    dup<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication") &(gene_data$node_age <= 296) )])) 
  }
  if(duplication_type=="old")
  { 
    dup<-NULL
    dup<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication") &(gene_data$node_age > 296))])) 
  }
  if(duplication_type=="all")
  { 
    dup<-NULL
    dup<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication"))])) 
  }
  pdup<-round(dup/Internal_event_number,2) ## Proportion of duplication event
  spe<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Speciation"))])) 
  pspe<-round(spe/Internal_event_number,2) ## Proportion of speciation events
  NA_event<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (is.na(gene_data$Event)))])) 
  pNA<-round(NA_event/Internal_event_number,2) ## Proportion of NA event
  
  
  ## Returning the tibble
  return(tibble(tip_num=ntips,
                internal_events=Internal_event_number,
                dup_num=dup,
                dup_prop=pdup,
                spe_num=spe,
                spe_prop=pspe,
                NA_prop=pNA,
                root_age=root_age,
                root_event=root_event))
  
}


## This function is written to paint trees using plotSimmap R function for three different states (speciation, duplication and NA)
paint_tree_mod<-function(tree)
{
  
  ## Considering empirical gene tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  ntips <- length(gene_tree$tip.label)
  Internal_node_data <- nrow(gene_data)-ntips
  Duplication_node_length <- length(gene_data$Event[which(gene_data$Event == "Duplication")])
  NA_node_length <- length(gene_data$Event[which(is.na(gene_data$Event))])
  Speciation_node_length <- length(gene_data$Event[which(gene_data$Event == "Speciation")])

  ## Now checking for trees with at least a duplication and a speciation events
  ## If does not match the criteria returns NA
  if((Internal_node_data - Duplication_node_length == 0) | (Internal_node_data - Speciation_node_length == 0))
  {
    return(NA)
  }
  else
  {
    if(NA_node_length > 0)
    {
      ##Identifying duplication nodes and edges to paint them
      dup_nodes <- gene_data$node[which(gene_data$Event=="Duplication")]
      NA_nodes <- gene_data$node[which(is.na(gene_data$Event))]
      NA_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% NA_nodes), 2])
      tree_painted <- paintBranches (gene_tree, edge=NA_edges, "NA", anc.state="S")
    }
    if(NA_node_length == 0)
    {

      ## Identifiying duplication nodes and edges to paint them
      dup_nodes <- gene_data$node[which(gene_data$Event=="Duplication")]
      dup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% dup_nodes), 2])
      tree_painted<-paintBranches (gene_tree, edge=dup_edges, "D", anc.state="S")
    }
    if("phylo" %in% class(tree_painted))
    {
      return(list(tree_painted, gene_tree, gene_data))
    }
  }
}

## This function is written to paint trees using the plotSimmap() for four different states (speciation, old duplication, young duplication and NA)
paint_tree_mod_new<-function(tree,event)
{
  count <<-count+1
  print(count)
  ## Considering real tree data
  gene_tree<-tree@phylo
  gene_data<-tree@data
  ntips <- length(gene_tree$tip.label) 
  Internal_node_data <- nrow(gene_data)-ntips
  Event <- NULL
  
  ## Checking for the empirical or the randomized data events
  if(event == "Empirical")
  {
    Event<-gene_data$Event
    gene_data$Event<-Event
  }
  if(event == "Random")
  {
    Event<-gene_data$event_new
    gene_data$Event<-Event
  }
  Duplication_node_length <- length(gene_data$Event[which(gene_data$Event == "Duplication")])
  Young_Duplication_node_length <- length(gene_data$Event[which(gene_data$Event == "Duplication" & gene_data$node_age <=296)])
  Old_Duplication_node_length <- length(gene_data$Event[which(gene_data$Event == "Duplication" & gene_data$node_age > 296)])
  NA_node_length <- length(gene_data$Event[which((is.na(gene_data$Event)) & (gene_data$node_age > 0))])
  Speciation_node_length <- length(gene_data$Event[which(gene_data$Event == "Speciation")])
  
  ## Now checking for trees with at least a duplication and a speciation events
  ## If it does not match the criteria, returns 'NA'
  if((Internal_node_data - Duplication_node_length == 0) | (Internal_node_data - Speciation_node_length == 0))
  {
    return(NA)
  }
  else
  {
    if(NA_node_length > 0) 
    {
      NA_nodes <- gene_data$node[which(is.na(gene_data$Event))]
      NA_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% NA_nodes), 2])
      tree_painted <- paintBranches (gene_tree, edge=NA_edges, "S_omit", anc.state="S")
      
      ## When both young and old duplication nodes are present
      if(Young_Duplication_node_length > 0 & Old_Duplication_node_length > 0)
      {
        ## Identifiying both types of duplication nodes to paint the edges 
        young_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age <=296)]
        old_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age > 296)]
        youngdup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% young_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=youngdup_edges, "DY", anc.state="S")
        olddup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% old_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=olddup_edges, "DO", anc.state="S")
      }
      
      ## If only young duplication nodes are present
      if(Young_Duplication_node_length > 0 & Old_Duplication_node_length == 0)
      {
        # Identifiying young duplication nodes to paint the edges 
        young_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age <=296)]
        youngdup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% young_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=youngdup_edges, "DY", anc.state="S")
      }
      
      ## If only old duplication nodes are present 
      if(Young_Duplication_node_length == 0 & Old_Duplication_node_length > 0)
      {
        # Identifiying old duplication nodes to paint the edges 
        old_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age > 296)]
        olddup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% old_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=olddup_edges, "DO", anc.state="S")
      }
    }
    if(NA_node_length == 0)
    {
      tree_painted<-gene_tree
      ## When both young and old duplication nodes are present
      if(Young_Duplication_node_length > 0 & Old_Duplication_node_length > 0)
      {
        # Identifiying both types of duplication nodes to paint the edges 
        young_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age <=296)]
        old_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age > 296)]
        youngdup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% young_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=youngdup_edges, "DY", anc.state="S")
        olddup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% old_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=olddup_edges, "DO", anc.state="S")
      }
      
      ## If only young duplication nodes are present
      if(Young_Duplication_node_length > 0 & Old_Duplication_node_length == 0)
      {
        # Identifiying young duplication nodes to paint the edges 
        young_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age <=296)]
        youngdup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% young_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=youngdup_edges, "DY", anc.state="S")
      }
      
      ## If only old duplication nodes are present 
      if(Young_Duplication_node_length == 0 & Old_Duplication_node_length > 0)
      {
        # Identifiying old duplication nodes to paint the edges 
        old_dup_nodes <- gene_data$node[which(gene_data$Event == "Duplication" & gene_data$node_age > 296)]
        olddup_edges <- unique(gene_tree$edge[which(gene_tree$edge[,1] %in% old_dup_nodes), 2])
        tree_painted<-paintBranches (tree_painted, edge=olddup_edges, "DO", anc.state="S")
      }
    }
    if("phylo" %in% class(tree_painted))
    {
      tree@phylo$painted_tree<- tree_painted
      #return(list(tree_painted, gene_tree, gene_data))
      return(tree)
    }
    else{return(NA)}
  }
}

## Functions required for plotting the painted trees
plot_phylogeny<-function(tree,fsize)
{
  ## Setting colors
  colSimmap <- c("#F8766D","#00BFC4","#C77CFF")
  names(colSimmap) <- c("S", "D","NA")
  p1<-plotSimmap(tree,col=colSimmap,fsize = fsize,ftype="b",mar=c(4,1,1,1))+
    axisPhylo(cex=1,font=2)
  return(p1)
}

## This function is written for model fitting on all the gene trees 
## We use the diagnostic tests included in the fitted models to avoid poor model choice
## We consider model with highest AICc weight as the best fit model
model_fitting_genetree_for_Tau<-function(input_tree,trait)
{
  i<<-i+1
  print(i)
  
  ## Extracting data
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau<-NULL
  
  ## Trait from empirical data
  if(trait=="Tau")
  {
    Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
    rownames(Tau)<-data$G[1:length(tree$tip.label)]
  }
  
  ## Trait data after permutations 
  if(trait=="rTau")
  {
    Tau <- as.matrix(data$Tau_new[1:length(tree$tip.label)])
    rownames(Tau)<-data$G[1:length(tree$tip.label)]
  }
  
  ## Generating dataframe for using the OUwie R package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ## Initializing all variables
  fit_OU1<-NULL
  fit_OUM<-NULL
  fit_BM1<-NULL
  fit_BMM<-NULL
  
  model_AICc<-NULL
  names(model_AICc)<-NULL
  
  model_support_package<-NULL
  
  ## Model fitting along with checking for the reliability of the model parameters 
  ## If the estimated parameter is available in both the "mvMORPH" and "OUwie" R packages, we consider the model fitting result from the mvMORPH R package. Otherwise, from the OUwie R package.
  ## We include a model only if the estimated model parameters are reliable
  
  ## Single optima OU model fitting
  fit_OU1<-tryCatch(mvOU(paint_tree, Tau, model="OU1",scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
  if((fit_OU1!="Error") && (fit_OU1$convergence == 0) && (fit_OU1$hess.values == 0))
  {
    fit_OU1$AICc<-fit_OU1$AICc
    names(fit_OU1$AICc)<-"OU1"
    model_AICc<-append(model_AICc,fit_OU1$AICc)
    model_support_package<-append(model_support_package,"mvMORPH")
  }
  if((fit_OU1=="Error") || (fit_OU1$convergence > 0) || (fit_OU1$hess.values > 0))
  { 
    fit_OU1<-NULL
    fit_OU1<-tryCatch(OUwie(paint_tree, data_OUWie,model="OU1",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F,diagn = T), error = function(e) {"Error"})
    if((fit_OU1!="Error") && (all(fit_OU1$eigval > 0)))
    {
      fit_OU1$AICc<-fit_OU1$AICc
      names(fit_OU1$AICc)<-"OU1"
      model_AICc<-append(model_AICc,fit_OU1$AICc)
      model_support_package<-append(model_support_package,"OUwie")
      # H0<-append(H0,fit_OU1$AICc)
      
      ## Adding for estimating model averaging
      #ou.results[["OU1"]]<-fit_OU1
      
    }
    if((fit_OU1 == "Error") || (any(fit_OU1$eigval < 0)))
    {
      fit_OU1<-NULL
    }
  }
  
  ## Single rate BM model fitting
  fit_BM1<-tryCatch(mvBM(paint_tree, Tau, model="BM1",scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
  if((fit_BM1!="Error") && (fit_BM1$convergence == 0) && (fit_BM1$hess.values == 0))
  {
    fit_BM1$AICc<-fit_BM1$AICc
    names(fit_BM1$AICc)<-"BM1"
    model_AICc<-append(model_AICc,fit_BM1$AICc)
    model_support_package<-append(model_support_package,"mvMORPH")
  }
  if((fit_BM1=="Error") || (fit_BM1$convergence > 0) || (fit_BM1$hess.values > 0))
  { 
    fit_BM1<-NULL
    fit_BM1<-tryCatch(OUwie(paint_tree, data_OUWie,model="BM1",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F,diagn = T), error = function(e) {"Error"})
    if((fit_BM1!="Error") && (all(fit_BM1$eigval > 0)))
    {
      fit_BM1$AICc<-fit_BM1$AICc
      names(fit_BM1$AICc)<-"BM1"
      model_AICc<-append(model_AICc,fit_BM1$AICc)
      model_support_package<-append(model_support_package,"OUwie")
    }
    if((fit_BM1 == "Error") || (any(fit_BM1$eigval < 0)))
    {
      fit_BM1<-NULL
    }
  }
  
  ## Multi rate BM model fitting
  fit_BMM<-tryCatch(mvBM(paint_tree, Tau, model="BMM",scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
  if((fit_BMM!="Error") && (fit_BMM$convergence == 0) && (fit_BMM$hess.values == 0))
  {
    fit_BMM$AICc<-fit_BMM$AICc
    names(fit_BMM$AICc)<-"BMM"
    model_AICc<-append(model_AICc,fit_BMM$AICc)
    model_support_package<-append(model_support_package,"mvMORPH")
  }
  if((fit_BMM=="Error") || (fit_BMM$convergence > 0) || (fit_BMM$hess.values > 0))
  { 
    fit_BMM<-NULL
    fit_BMM<-tryCatch(OUwie(paint_tree, data_OUWie,model="BMS",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F,diagn = T), error = function(e) {"Error"})
    if((fit_BMM!="Error") && (all(fit_BMM$eigval > 0)))
    {
      fit_BMM$AICc<-fit_BMM$AICc
      names(fit_BMM$AICc)<-"BMM"
      model_AICc<-append(model_AICc,fit_BMM$AICc)
      model_support_package<-append(model_support_package,"OUwie")
    }
    if((fit_BMM == "Error") || (any(fit_BMM$eigval < 0)))
    {
      fit_BMM<-NULL
    }
  }
  
  ## Multi optima OU model fitting
  fit_OUM<-tryCatch(mvOU(paint_tree, Tau, model="OUM",scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
  if((fit_OUM!="Error") && (fit_OUM$convergence == 0) && (fit_OUM$hess.values == 0))
  {
    fit_OUM$AICc<-fit_OUM$AICc
    names(fit_OUM$AICc)<-"OUM"
    model_AICc<-append(model_AICc,fit_OUM$AICc)
    model_support_package<-append(model_support_package,"mvMORPH")
  }
  if((fit_OUM=="Error") || (fit_OUM$convergence > 0) || (fit_OUM$hess.values > 0))
  { 
    fit_OUM<-NULL
    fit_OUM<-tryCatch(OUwie(paint_tree, data_OUWie,model="OUM",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F,diagn = T), error = function(e) {"Error"})
    if((fit_OUM!="Error") && (all(fit_OUM$eigval > 0)))
    {
      fit_OUM$AICc<-fit_OUM$AICc
      names(fit_OUM$AICc)<-"OUM"
      model_AICc<-append(model_AICc,fit_OUM$AICc)
      model_support_package<-append(model_support_package,"OUwie")
    }
    if((fit_OUM == "Error") || (any(fit_OUM$eigval < 0)))
    {
      fit_OUM<-NULL
    }
  }
  
  ## Now we need to check if "model_AICc" is not empty
  if(length(model_AICc) > 0){
    
    ## Creating a dataframe of AICc weight
    model_AICc_df<-data.frame(Model=names(model_AICc),AICc=model_AICc,model_supprt_package=model_support_package)
    model_AICc_df$delta_model_AICc<-model_AICc_df$AICc-min(model_AICc_df$AICc)
    model_AICc_df$exp_AICc <- exp(-0.5 * model_AICc_df$delta_model_AICc)
    model_AICc_df$weights_AICc <- round(model_AICc_df$exp_AICc/sum(model_AICc_df$exp_AICc),3)
    
    ## Ordering dataframe based on AICc weights
    model_AICc_df<-model_AICc_df[order(-model_AICc_df$weights_AICc),] 
    
    ## Identifying best fit model 
    best_fit_model<-model_AICc_df$Model[1]
    
    return(tibble(tree_index=i,best_fit_model=best_fit_model,best_fit_model_AICc=model_AICc_df$AICc[1],best_fit_model_AICc_weight=model_AICc_df$weights_AICc[1],model_support_package=model_AICc_df$model_supprt_package[1]))
  }
  
  if(length(model_AICc) == 0){
    return(tibble(tree_index=i,best_fit_model="None",best_fit_model_AICc=NA,best_fit_model_AICc_weight=NA,model_support_package=NA))
  }
}

## Adding best fit model parameters into the '@phylo' slots of the trees
add_model_parameters_info <- function(input_tree,dataframe)
{
  i<<-i+1
  print(i)
  
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Prapering dataframe for using OUwie package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ## Checking for successful model fitting 
  if(dataframe$best_fit_model[i] != "None")
  {
    input_tree@phylo$best_fit_model<-dataframe$best_fit_model[i]
    input_tree@phylo$best_fit_model_AICc<-dataframe$best_fit_model_AICc[i]
    input_tree@phylo$best_fit_model_AICc_weight<-dataframe$best_fit_model_AICc_weight[i]
    source<-dataframe$model_support_package[i]
    
    ## Extracting parameters 
    if(source == "OUwie" && dataframe$best_fit_model[i]=="BM1")
    {
      fit_model<-NULL
      fit_model<-OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_BM1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="BM1")
    {
      fit_model<-NULL
      fit_model<-mvBM(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_BM1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="BMM")
    {
      fit_model<-NULL
      fit_model<-OUwie(paint_tree, data_OUWie,model="BMS",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_BMM_parameters(input_tree,"BMM")
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="BMM")
    {
      fit_model<-NULL
      fit_model<-mvBM(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_BMM_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="OU1")
    {
      fit_model<-NULL
      fit_model<-OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_OU1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="OU1")
    {
      fit_model<-NULL
      fit_model<-mvOU(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_OU1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="OUM")
    {
      fit_model<-NULL
      fit_model<-OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_OUM_parameters(input_tree,"OUM")
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="OUM")
    {
      fit_model<-NULL
      fit_model<-mvOU(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE)
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_OUM_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)
    }
    
  }
  
  if(dataframe$best_fit_model[i] == "None")
  {
    input_tree@phylo$best_fit_model<-dataframe$best_fit_model[i]
    input_tree@phylo$best_fit_model_AICc<-NA
    input_tree@phylo$best_fit_model_AICc_weight<-NA
    source<-NA
    input_tree@phylo$fitted_model<-NA
    input_tree@phylo$fitted_parameters<-NA
    return(input_tree)
  }
}


## This function is written to extract parameters from the fitted OU1 model of the OUwie R package
OUwie_OU1_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Prapering dataframe for using the OUwie R package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-OUwie(paint_tree, data_OUWie,model="OU1",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=colnames(fit_model$solution),alpha=fit_model$solution[1,],
             sig2=fit_model$solution[2,],theta=fit_model$theta[,1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=as.numeric(df$alpha[1]),
                  alpha_duplication_old= as.numeric(df$alpha[1]),
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=as.numeric(df$theta[1])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=as.numeric(df$alpha[1]),
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=NA,
                  alpha_duplication_old= as.numeric(df$alpha[1]),
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[1])))
  }
}

## This function is written to extract parameter values from the fitted OU1 model of the mvMORPH R package
mvMORPH_OU1_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-mvOU(paint_tree, Tau, model="OU1",scale.height = F,echo=FALSE)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(alpha=fit_model$alpha,sig2=fit_model$sigma,theta=fit_model$theta[1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=as.numeric(df$alpha[1]),
                  alpha_duplication_old= as.numeric(df$alpha[1]),
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=as.numeric(df$theta[1])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=as.numeric(df$alpha[1]),
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[1]),
                  alpha_duplication_young=NA,
                  alpha_duplication_old= as.numeric(df$alpha[1]),
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[1])))
  }
}

## This function is written to extract parameter values from the fitted OUM model of the OUwie R package
OUwie_OUM_parameters<-function(input_tree,model)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Prapering dataframe for using the OUwie R package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  if(model %in% "OUM")
  {
    fit_model<-NULL
    fit_model<-OUwie(paint_tree, data_OUWie,model="OUM",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
  }
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=colnames(fit_model$solution),alpha=fit_model$solution[1,],
             sig2=fit_model$solution[2,],theta=fit_model$theta[,1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=as.numeric(df$alpha[df$event %in% "DY"]),
                  alpha_duplication_old= as.numeric(df$alpha[df$event %in% "DO"]),
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=as.numeric(df$alpha[df$event %in% "DY"]),
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=NA,
                  alpha_duplication_old= as.numeric(df$alpha[df$event %in% "DO"]),
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }
}


## This function is written to extract data from the fitted OUM model of the mvMORPH R package
mvMORPH_OUM_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-mvOU(paint_tree, Tau, model="OUM",scale.height = F,echo=FALSE)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=rownames(fit_model$theta),alpha=fit_model$alpha,sig2=fit_model$sigma,theta=fit_model$theta)
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=as.numeric(df$alpha[df$event %in% "DY"]),
                  alpha_duplication_old= as.numeric(df$alpha[df$event %in% "DO"]),
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=as.numeric(df$alpha[df$event %in% "DY"]),
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=as.numeric(df$alpha[df$event %in% "S"]),
                  alpha_duplication_young=NA,
                  alpha_duplication_old= as.numeric(df$alpha[df$event %in% "DO"]),
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }
}


## This function is written to extract parameter values from the fitted BM1 model of the OUwie R package
OUwie_BM1_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Prapering dataframe for using the OUwie R package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-OUwie(paint_tree, data_OUWie,model="BM1",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=colnames(fit_model$solution),alpha=fit_model$solution[1,],
             sig2=fit_model$solution[2,],theta=fit_model$theta[1,1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=as.numeric(df$theta[1])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[1])))
  }
}

## This function is written to extract parameters from the fitted BM1 model of the mvMORPH R package
mvMORPH_BM1_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-mvBM(paint_tree, Tau, model="BM1",scale.height = F,echo=FALSE)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(sig2=fit_model$sigma,theta=fit_model$theta[1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=as.numeric(df$theta[1])))
  }  
  
  ## If only the new duplication events are present, but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=as.numeric(df$sig2[1]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=as.numeric(df$theta[1]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present, but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[1]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[1]),
                  theta_speciation=as.numeric(df$theta[1]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[1])))
  }
}


## This function is written to extract data from the fitted BMM model of the OUwie R package
OUwie_BMM_parameters<-function(input_tree,model)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Preparing dataframe for using the OUwie R package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  if(model == "BMM")
  {
    fit_model<-NULL
    fit_model<-OUwie(paint_tree, data_OUWie,model="BMS",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F)
  }
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=colnames(fit_model$solution),alpha=fit_model$solution[1,],sig2=fit_model$solution[2,],theta=fit_model$theta[1,1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }  
  
  ## If only the new duplication events are present but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }
}

## This function is written to extract parameters from the fitted BMM model of the mvMORPH R package
mvMORPH_BMM_parameters<-function(input_tree)
{
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ##Checking for young duplication (<= old speciation age) and old duplication events (> old speciation age)
  Young_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age <=296)])
  Old_Duplication_node_length <- length(data$Event[which(data$Event == "Duplication" & data$node_age >296)])
  
  fit_model<-NULL
  fit_model<-mvBM(paint_tree, Tau, model="BMM",scale.height = F,echo=FALSE)
  
  ## Preparing dataframe
  df<-NULL
  df<-tibble(event=names(fit_model$sigma[,,]),sig2=fit_model$sigma[,,],theta=fit_model$theta[1,1])
  
  ## If both the old and new duplication events are present
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }  
  
  ## If only the new duplication events are present but not old
  if(Young_Duplication_node_length > 0 && Old_Duplication_node_length == 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=as.numeric(df$sig2[df$event %in% "DY"]),
                  sig2_duplication_old=NA,
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=as.numeric(df$theta[df$event %in% "DY"]),
                  theta_duplication_old=NA))
  }
  
  ## If only the old duplication events are present but not new
  if(Young_Duplication_node_length == 0 && Old_Duplication_node_length > 0)
  {
    return(tibble(alpha_speciation=NA,
                  alpha_duplication_young=NA,
                  alpha_duplication_old= NA,
                  sig2_speciation=as.numeric(df$sig2[df$event %in% "S"]),
                  sig2_duplication_young=NA,
                  sig2_duplication_old=as.numeric(df$sig2[df$event %in% "DO"]),
                  theta_speciation=as.numeric(df$theta[df$event %in% "S"]),
                  theta_duplication_young=NA,
                  theta_duplication_old=as.numeric(df$theta[df$event %in% "DO"])))
  }
}

## This function is written to add the duplication (young and old), speciation, and NA events information for each tree
tree_added_events_info<-function(tree)
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
  ndup<-0
  ndup_young<-NA
  ndup_old<-NA
  nspe<-0
  nNA<-0
  nNA_young<-NA
  nNA_old<-NA
  
  ##Since, many internal node events are assigned as "NA", sum of the proportions of speciation and duplication events may not be equal to 1
  ndup<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication"))])) 
  ndup_young<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication") & (gene_data$node_age <= 296))]))
  ndup_old<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Duplication") & (gene_data$node_age > 296))])) 
  nspe<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (gene_data$Event %in% "Speciation"))])) 
  nNA<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (is.na(gene_data$Event)))]))
  nNA_young<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (is.na(gene_data$Event)) & (gene_data$node_age <= 296))]))
  nNA_old<-as.numeric(length(gene_data$node[which((!is.na(gene_data$pic)) & (is.na(gene_data$Event)) & (gene_data$node_age > 296))]))
  tree@phylo$n_speciation<- nspe
  tree@phylo$n_duplication_young<- ndup_young
  tree@phylo$n_duplication_old<- ndup_old
  tree@phylo$n_NA_young<- nNA_young
  tree@phylo$n_NA_old<- nNA_old
  return(tree)
  
}

## Adding best fit model parameters into the '@phylo' slots of the randomized trees
add_model_parameters_info2 <- function(input_tree,dataframe)
{
  i<<-i+1
  print(i)
  
  ## Extracting tree info
  paint_tree<-input_tree@phylo$painted_tree
  tree<-input_tree@phylo
  data<-input_tree@data
  Tau <- as.matrix(data$Tau[1:length(tree$tip.label)])
  rownames(Tau)<-data$G[1:length(tree$tip.label)]
  
  ## Preparing dataframe for using OUwie package
  trait <- Tau
  tip<-rownames(trait)
  x<-getStates(paint_tree,"tips")
  data_OUWie <- data.frame(tip, times=x,trait)
  root.age<-max(data$node_age)
  
  ## Checking for the successful model fitting 
  if(dataframe$best_fit_model[i] != "None")
  {
    input_tree@phylo$best_fit_model<-dataframe$best_fit_model[i]
    input_tree@phylo$best_fit_model_AICc<-dataframe$best_fit_model_AICc[i]
    input_tree@phylo$best_fit_model_AICc_weight<-dataframe$best_fit_model_AICc_weight[i]
    source<-dataframe$model_support_package[i]
    
    ## Extracting fitted parameters 
    if(source == "OUwie" && dataframe$best_fit_model[i]=="BM1")
    {
      fit_model<-NULL
      fit_model<-tryCatch(OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_BM1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="BM1")
    {
      fit_model<-NULL
      fit_model<-tryCatch(mvBM(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_BM1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="BMM")
    {
      fit_model<-NULL
      fit_model<-tryCatch(OUwie(paint_tree, data_OUWie,model="BMS",simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_BMM_parameters(input_tree,"BMM")
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="BMM")
    {
      fit_model<-NULL
      fit_model<-tryCatch(mvBM(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_BMM_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="OU1")
    {
      fit_model<-NULL
      fit_model<-tryCatch(OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_OU1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="OU1")
    {
      fit_model<-NULL
      fit_model<-tryCatch(mvOU(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_OU1_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "OUwie" && dataframe$best_fit_model[i]=="OUM")
    {
      fit_model<-NULL
      fit_model<-tryCatch(OUwie(paint_tree, data_OUWie,model=dataframe$best_fit_model[i],simmap.tree=T,root.age = root.age,scaleHeight =F,quiet = T,warn = F), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-OUwie_OUM_parameters(input_tree,"OUM")
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
    if(source == "mvMORPH" && dataframe$best_fit_model[i]=="OUM")
    {
      fit_model<-NULL
      fit_model<-tryCatch(mvOU(paint_tree, Tau, model=dataframe$best_fit_model[i],scale.height = F, diagnostic=T, echo=FALSE), error = function(e) {"Error"})
      if(fit_model!="Error"){
      input_tree@phylo$fitted_model<-fit_model
      fitted_parameters<-NULL
      fitted_parameters<-mvMORPH_OUM_parameters(input_tree)
      input_tree@phylo$fitted_parameters<-fitted_parameters
      return(input_tree)}
      else{return(NA)}
    }
    
  }
  
  if(dataframe$best_fit_model[i] == "None")
  {
    input_tree@phylo$best_fit_model<-dataframe$best_fit_model[i]
    input_tree@phylo$best_fit_model_AICc<-NA
    input_tree@phylo$best_fit_model_AICc_weight<-NA
    source<-NA
    input_tree@phylo$fitted_model<-NA
    input_tree@phylo$fitted_parameters<-NA
    return(input_tree)
  }
}

## This function is written to use bayou R package to identify regime shifts in a Baysian framework
## We here use a posterior probability cutoff >0.1 to detect and to export all the regime shifts
## Among them, we consider a statistical support for regime shift in the posterior probability of >= 0.7 in our manuscript
bayou_analysis<-function(tree)
{
  count<<-count+1
  print(count)
  tre<-tree
  OUTree<-tre@phylo
  OUTree <- reorder(OUTree, "postorder")

  #To remove trees with negative branch length
  negative<-0
  negative<- OUTree$edge.length[OUTree$edge.length < 0]
  if(length(negative)==0){

  OUdata<-tre@data$Tau[is.tip.nhx(tre)]
  names(OUdata)<-tre@data$label[is.tip.nhx(tre)]


  ## Setting prior
  priorOU <-NULL
  priorOU <- make.prior(OUTree , dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dsb="dsb", dk="cdpois", dtheta="dnorm"), param=list(dalpha=list(scale=1), dsig2=list(scale=1), dk=list(lambda=15, kmax=200), dsb=list(bmax=1,prob=1), dtheta=list(mean=mean(OUdata), sd=2)))


  ## Running MCMC
  mcmcOU <-NULL
  mcmcOU <- tryCatch(bayou.makeMCMC(OUTree, OUdata, SE=0, model="OU", prior=priorOU, plot.freq=5000, ticker.freq=1000),error = function(e) {"Error"}) # Set up the MCMC

  if(mcmcOU!="Error"){
  mcmcOU$run(100000) # Running the MCMC

  #Loading result
  chainOU <- mcmcOU$load()

  # Setting a "burnin" parameter that tells the package **coda** to discard the first 30% of the chain.
  chainOU <- set.burnin(chainOU, 0.3)

  ## Required dataframe
  BayOU.output<-NULL
  BayOU.output<-summary(chainOU)["branch.posteriors"]
  BayOU.output<-BayOU.output$branch.posteriors
  BayOU.output$branch<-as.numeric(rownames(BayOU.output))
  BayOU.output<-BayOU.output[BayOU.output$pp>0.1,] ## Considering data with posterior probability >= 0.1
  if(nrow(BayOU.output)>0){
        BayOU.output$node1<-OUTree$edge[BayOU.output$branch,1]
        BayOU.output$node2<-OUTree$edge[BayOU.output$branch,2]
        BayOU.output$Event<-tre@data$Event[BayOU.output$node1]
        BayOU.output$node_age<-tre@data$node_age[BayOU.output$node1]
        BayOU.output$index<-as.numeric(count)
        return(BayOU.output)}}}
}

## This function is written to calculate proportions of adaptive regime shifts for the speciation and the duplication events in the trees
## datatype can be "all", "young","old"
proportion_calculation<-function(Treeset,datatype,dataframe)
{
  ##Numbers of speciation and duplication events for trees with regime shifts
  regime_shift_info<-NULL
  regime_shift_info<-bind_rows(lapply(Treeset, tree_data_statistics,datatype))
  regime_shift_info$index<-as.numeric(rownames(regime_shift_info))
  #regime_shift_info$spe_num<-regime_shift_info$internal_events-regime_shift_info$dup_num
  regime_shift_info<-regime_shift_info[c(10,5,3,1,2)]

  ## Counting frquency of regime shift  for "speciation" and "duplication" events
  regime_shift_duplication<-
    dataframe  %>%
    filter(Event=="Duplication") %>%
    .$index %>% 
    table() %>% 
    as.data.frame()
  colnames(regime_shift_duplication)<-c("index","regime_shift_duplication")

  regime_shift_speciation<-
    dataframe%>%
    filter( Event=="Speciation" ) %>%
    .$index %>% 
    table() %>% 
    as.data.frame()
  colnames(regime_shift_speciation)<-c("index","regime_shift_speciation") 

  ## Merging datasets to get actual statististics of regime shift
  merged_regime_shift<-full_join(regime_shift_speciation,regime_shift_duplication)
  info_regime_shift<-merge(regime_shift_info,merged_regime_shift, by = c("index"))
  info_regime_shift[is.na(info_regime_shift)] <- 0
  info_regime_shift_final<-info_regime_shift[which(info_regime_shift$dup_num>0 & info_regime_shift$spe_num>0),]

  ## Proportion of regime-shifts
  info_regime_shift_final$shift_prob_dup<-(info_regime_shift_final$regime_shift_duplication/(info_regime_shift_final$dup_num*2))
  info_regime_shift_final$shift_prob_spe<-(info_regime_shift_final$regime_shift_speciation/(info_regime_shift_final$spe_num*2))

return(info_regime_shift_final) 
}

