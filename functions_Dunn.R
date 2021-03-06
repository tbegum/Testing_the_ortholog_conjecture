
#' Create a string from a p.value
#' 
#' @param x htest object, like that returned from cor.test()
#' @param digits the number of significant digits
#' @return A character object
#' @export
format_p = function ( x, digits=3 ){
		fp <- format.pval( x$p.value, digits = max(1L, digits - 3L) )
	out = if(substr(fp, 1L, 1L) == "<") fp else paste("=",fp)
		return( out )
}

#' Parses nhx text to a treeio::treedata object
#' 
#' @param tree_text Character string repressenting nhx tree
#' @return A treeio::treedata object
#' @export
parse_gene_trees = function( tree_text ){
	tree_tc = textConnection( tree_text )
	tree = treeio::read.nhx( tree_tc )
	close( tree_tc )
	
	# Parse clade labels from phylo object into the data frame
	tree@data$label = c( tree@phylo$tip.label, tree@phylo$node.label )
	
	# Compara trees sometimes have speciation nodes whose descendants have the 
	# same clade name. This isn't biologically possible, so change such nodes 
	# from speciation events to NA so they don't interfere with tree calibration 
	# and are not used for calculating speciation contrasts
	
	n_nodes = nrow( tree@data )
	n_tips = length( tree@phylo$tip.label )
	internal_nodes = ( n_tips + 1 ):n_nodes
	is_speciation = tree@data$D == "N"
	is_speciation[ is.na( is_speciation ) ] = FALSE
	internal_speciation_nodes = tree@data$node[ ( tree@data$node > n_tips ) & is_speciation ]
	
	# Create a vector of clade names of speciation events, with NA for all other nodes
	speciation_names = rep( NA, n_nodes )
	speciation_names[ internal_speciation_nodes ] = tree@data$label[ internal_speciation_nodes ]
	
	# Loop over the internal speciation nodes
	for( i in internal_speciation_nodes ){
		# get descendent internal nodes
		descendants = hutan::descendants( tree@phylo, i )
		descendants = descendants[ descendants %in% internal_nodes ]
		
		descendant_names = speciation_names[ descendants ]
		
		if( speciation_names[i] %in% descendant_names ){
			tree@data$D[i] = NA
		}
	}
	
	
	# Create a human readable Event column
	tree@data$Event = NA
	tree@data$Event[ tree@data$D == "N" ] = "Speciation"
	tree@data$Event[ tree@data$D == "Y" ] = "Duplication"
	tree@data$Event = factor( tree@data$Event, levels=c( "Speciation", "Duplication" ) )
	
	return( tree )
}


#' Adds some random noise from a normal distribution to calibration
#' times
#' 
#' @param calibration_times A dataframe with two columns: age.min is the age 
#' of a clade, and clade is the name of the clade
#' @param sd_fraction this is multiplied by each value to set sd of that value 
#' @return A dataframe as above, but with noise
#' @export
add_noise_calibration_times = function( calibration_times, sd_fraction ) {
	age = calibration_times$age
	new_age = NA
	
	# Keep trying until dates are coherent with tree structure
	# This is a ladder tree, so that means checking that they 
	# are still in ascending order
	while( NA %in% new_age ) {
		x = sapply( 
			age,
			function( i ){
				rnorm( 1, mean = i, sd = ( sd_fraction * i ) )
			}
		)
		
		if( ( ! is.unsorted( x ) ) & all( x > 0 ) ){
			new_age = x
		}
		
	}
	
	calibration_times$age = new_age
	return( calibration_times )
}	


#' Get a boolean vector corresponding to all tips and internal nodes of 
#' a tree, with value TRUE for tips and FALSE for internal nodes
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object
#' @return A boolean vector
#' @export
is.tip.nhx = function( nhx ) {
	is.tip = rep( FALSE, nrow( nhx@data ) )
	is.tip[ 1:length( nhx@phylo$tip.label ) ] = TRUE
	is.tip
}


#' Fix ensembl clade names.
#' Homininae is pan + homo + gorilla
#' Hominini is pan + homo
#' ENSEMBL Compara trees lislabel Hominini as Homininae
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object
#' @return A phylogenetic tree as a treeio::treedata object
#' @export
fix_hominini = function( nhx ) {
	# Get the numbers for all nodes labeled as homininae
	homininae_nodes = nhx@data$node[ nhx@data$label == "Homininae" ]
	
	# Check which have no gorilla as a descendent
	no_gorilla = sapply(
		homininae_nodes,
		function( x ){
			tips = hutan::tip_descendants( nhx@phylo, x )
			tip_species = nhx@data$species[ tips ]
			return( ! "gorilla_gorilla" %in% tip_species )
		}
	)
	
	if( any( no_gorilla ) ){
		hominini_nodes = homininae_nodes[ no_gorilla ]
		nhx@data$label[ hominini_nodes ] = "Hominini"
	}
	
	return( nhx )
}


#' Record the distance from a tip for each node in a tree in 
#' the data slot
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata
#' @return A treeio::treedata object, with an additional node_age column
#' in the data slot, or the original object if it was not of class treedata
#' @export
store_node_age = function( nhx ) {
	
	if ( class( nhx ) != "treedata" ) {
		return( nhx )
	}
	
	node_age = hutan::distance_from_tip( nhx@phylo )
	
	# make sure the dataframe is ordered by consecutive nodes
	stopifnot( all( nhx@data$node == 1:length( nhx@data$node ) ) )
	
	nhx@data$node_age = node_age 
	
	return( nhx )
}

#' Get the number of speciation events in an nhx tree
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object, with a 
#' column @data$D that has value N for speciation events
#' @return An integer indicating the number of speciation events
#' @export
get_n_speciation = function( nhx ){
	sum( nhx@data$D[ ! is.tip.nhx( nhx ) ] == "N", na.rm = TRUE )
}

#' Drop tips without expression data from a tree. 
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object, with a 
#' column @data$Tau containing expression data
#' @param min_genes_with_expression The minimum number of tips with 
#' expression data for a tree to be retained
#' @return A treeio::treedata object, or NA if less than 
#' min_genes_with_expression
#' @export
drop_empty_tips = function( nhx, min_genes_with_expression ) {
	
	# Identify the tips without expression data
	to_drop = which( is.na( nhx@data[ 1:length( nhx@phylo$tip.label ), ]$Tau ) )
	
	# Return the pruned tree if it meets sampling criteria, otherwise return NA
	remaining = length( nhx@phylo$tip.label ) - length( to_drop )
	if( remaining >= min_genes_with_expression ){
		pruned = treeio::drop.tip( nhx, to_drop )
		
		node_depth = ape::node.depth( pruned@phylo )
		
		# Add node depth
		pruned@data$node_depth = node_depth
		
		return( pruned )
	}
	else{
		return( NA )
	}
}


#' Adjust branch lengths of a phylogenetic tree to make it ultrametric and 
#' to time calibrate the speciation nodes to fixed values
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object, with speciation 
#' nodes annotated with a value of "N" in @data$D
#' @param calibration_times A dataframe with two columns: age is the age 
#' of a clade, and clade is the name of the clade
#' @param ... Any additional arguments to pass to ape::chronos()
#' @return A treeio::treedata object if successfully calibrated
#' @export
calibrate_tree = function ( nhx, calibration_times, ... ) {
	
	# Create calibration matrix for speciation nodes
	calibration = 
		nhx@data[ !is.tip.nhx( nhx ), ] %>%
		filter( D == "N" ) %>%
		left_join( calibration_times, c( "label" = "clade" ) ) %>%
		mutate( age.min = age ) %>%
		mutate( age.max = age ) %>% 
		select( node, age.min, age.max ) %>%
		mutate( soft.bounds = NA )
	
	tree = try( 
		ape::chronos( nhx@phylo, calibration=calibration, ... ) 
	)
	
	if( "phylo" %in% class( tree ) ){
		class( tree ) = "phylo"
		nhx@phylo = tree
		return( nhx )
	}
	else{
		return( NA )
	}
}

#' Calculate phylogenetic independent contrasts (PIC) with the tree and 
#' character data in a treeio::treedata object
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object, with a 
#' column @data$Tau containing expression data
#' @param model_method The model of trait evolution. Can be one of c("BM", "OU")
#' @return A treeio::treedata object, with new @data columns pic and var_exp
#' @export
pic.nhx = function( nhx, model_method="BM" ) {
	
	tau = nhx@data$Tau[ is.tip.nhx( nhx ) ]
	
	if( any( is.na(tau) ) ){
		stop( "A trait value is NA." )
	}
	
	# Calculate the contrasts
	if ( model_method=="BM" ){
		p = ape::pic( 
			tau, 
			nhx@phylo, 
			var.contrasts=TRUE 
		)
	} else if ( model_method=="OU" ){
		p = hutan::picx( 
			tau, 
			nhx@phylo, 
			var.contrasts=TRUE , 
			model_method="OU", 
			model_parameters=nhx@phylo$model_ou
		)
	} else {
		stop("ERROR: Invalid model_method")
	}
	
	# Remove previous values if they exist
  nhx@data$pic = NULL
  nhx@data$var_exp = NULL

	# Add the results back to the @data slot, padding the rows that correspond to 
	# tips with NA
	#contrast<-NULL
	#contrast<-as.data.frame(p)
	#contrast$label<-row.names(contrast)
	#colnames(contrast)<-c("pic","var_exp","label")
	#nhx@data <-as.data.frame(left_join(nhx@data,contrast,by=c("label")))
  nhx@data$pic <- c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,1 ] )
	nhx@data$var_exp <- c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,2 ] )
	
	return( nhx )
}

# 
# pic.nhx.pre = function( nhx ) {
#   
#   tau = nhx@data$Tau[ is.tip.nhx( nhx ) ]
#   
#   if( any( is.na(tau) ) ){
#     stop( "A trait value is NA." )
#   }
#   
#   # Calculate the contrasts
#   p = pic( tau, nhx@phylo, var.contrasts=TRUE )
#   
#   # Add the results back to the @data slot, padding the rows that correspond to 
#   # tips with NA
#   nhx@data$pic = c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,1 ] )
#   nhx@data$var_exp = c( rep( NA, length( nhx@phylo$tip.label ) ), p[ ,2 ] )
#   
#   return( nhx )
# }

#' Add phylogenetic independent contrasts (PIC) with the tree and 
#' character data in a treeio::treedata object
#' 
#' @param gene_trees_calibrated A list of time calibrated phylogenetic trees 
#' as treeio::treedata objects 
#' @param model_method The model of trait evolution. Can be one of c("BM", "OU")
#' @return A list of time phylogenetic trees with pics 
#' as treeio::treedata objects
#' @export
add_pics_to_trees = function( gene_trees_calibrated, model_method="BM" ) {
	
	# Calculate independent contrasts for tau on each tree, storing the results 
	# back into the @data slot of the tree objects
	gene_trees_pic = foreach( tree=gene_trees_calibrated ) %dopar% 
		pic.nhx( tree, model_method=model_method )
	
	return( gene_trees_pic )
}


#' Add phylogenetic independent contrasts (PIC) with the tree and 
#' character data in a treeio::treedata object
#' 
#' @param gene_trees_pic A list of phylogenetic trees with pics
#' as treeio::treedata objects 
#' @return A tibble with combined data, including phylogenetic independent 
#' contrasts, for internal nodes of all trees
#' @export
summarize_contrasts = function( gene_trees_pic ) {
	# Combine the @data slots across all trees. These have labels, node annotations,
	# expression data, contrast results, etc...
	# Also take the absolute value of the contrast, since we only consider magnitude.
	nodes_all = 
		lapply( 
			gene_trees_pic, 
			function( nhx ){
				tags = nhx@data
				tags$gene = digest( nhx ) # Creates a hash that is unique to each gene tree
				tags %<>% select( -( B ) ) # B has inconsistent types, remove it
				return( tags )
			}
		) %>%
		bind_rows() %>%
		mutate( pic = abs( pic ) )
	
	# Create a subset that corresponds to nodes with contrasts. This excludes tips.
	nodes_contrast = 
		nodes_all %>%
		filter( ! is.na( pic ) ) %>%
		filter( ! is.na( D ) )
	
	return( nodes_contrast )
}

#' Summarize edge data and properties of a single tree
#' 
#' @param gene_trees A phylogenetic tree
#' as treeio::treedata object
#' @return A tibble with edge data
summarize_nhx_edges = function( nhx ){
	tags = nhx@data
	phy = nhx@phylo
	edge_length = phy$edge.length
	
	tau_tips = nhx@data$Tau[1:length(phy$tip.label)]
	tau_internal = ace( tau_tips, phy, method="pic" )$ace
	tau = c( tau_tips, tau_internal )
	tau_on_edges = sapply( phy$edge, function(x){ tau[x] })
	dim( tau_on_edges ) = dim( phy$edge )
	tau_parent = tau_on_edges[ , 1 ]
	tau_child = tau_on_edges[ , 2 ]
	
	event_type = as.character( tags$Event )
	event_type[1:length(phy$tip.label)] = "Tip"
	events_on_edges = sapply( phy$edge, function(x){ event_type[x] })
	dim( events_on_edges ) = dim( phy$edge )
	
	scaled_change = (tau_child - tau_parent) / edge_length
	
	tibble(
		gene = digest( nhx ),
		edge_length = edge_length,
		tau_parent = tau_parent,
		tau_child = tau_child,
		event_parent = events_on_edges[ , 1 ],
		event_child = events_on_edges[ , 2 ],
		scaled_change = scaled_change
	)

}

#' Summarize edge data and properties of a single tree
#' 
#' @param gene_trees A list of phylogenetic trees 
#' as treeio::treedata objects 
#' @return A tibble with combined edge data
#' @export
summarize_edges = function( gene_trees ){

	edge_summaries = foreach( tree=gene_trees ) %dopar% 
		summarize_nhx_edges( tree )

	return( edge_summaries %>% bind_rows() )

}


#' Prepare summary statistics for each tree
#' 
#' @param gene_trees A list of phylogenetic trees with pics
#' as treeio::treedata objects 
#' @return A tibble with combined data
#' @export
summarize_trees = function( gene_trees_pic ) {
	
	tree_summary = 
		lapply( 
			gene_trees_pic, 
			function( nhx ){
				tags = nhx@data
				phy = nhx@phylo
				x = nhx@data$Tau[1:length(phy$tip.label)]
				names( x ) = phy$tip.label
				
				tibble(
					gene = digest( nhx ),
					n_tips = length( x ),
					tau_mean = mean( x ),
					tau_var = var( x ),
					K = phytools::phylosig( phy, x, method="K" ),
					z0_bm = nhx@phylo$model_bm$opt$z0,
					z0_ou = nhx@phylo$model_ou$opt$z0,
					sigsq_bm = nhx@phylo$model_bm$opt$sigsq,
					sigsq_ou = nhx@phylo$model_ou$opt$sigsq,
					alpha_ou = nhx@phylo$model_ou$opt$alpha,
					aic_bm = nhx@phylo$model_bm$opt$aic,
					aic_ou = nhx@phylo$model_ou$opt$aic,
					lnL_bm = nhx@phylo$model_bm$opt$lnL,
					lnL_ou = nhx@phylo$model_ou$opt$lnL
				)
			}
		) %>%
		bind_rows()
	
	return( tree_summary )
}


#' Calibrate a list of gene trees given a set of speciation node dates
#' 
#' @param gene_trees_pruned A list of phylogenetic trees as treeio::treedata objects 
#' with associated character data
#' @param calibration_times A dataframe with two columns: age.min is the age 
#' of a clade, and clade is the name of the clade
#' @param ... Any additional arguments to pass to ape::chronos()
#' @return A list of calibrated phylogenetic trees as treeio::treedata objects. May be 
#' shorter than gene_trees_pruned if some calibrations fail.
#' @export
calibrate_trees = function( gene_trees_pruned, calibration_times, ... ) {
	# Make the trees ultrametric and calibrate the speciation nodes to 
	# specified times	
	gene_trees_calibrated = foreach( tree=gene_trees_pruned ) %dopar% 
		calibrate_tree( tree, calibration_times=calibration_times, ... )

	# Remove trees that could not be successfully calibrated
	gene_trees_calibrated = gene_trees_calibrated[ ! is.na( gene_trees_calibrated ) ]
	
	# Parse the calibrated node ages from the internal @phylo object and store them with the 
	# corresponding rows in the @data object
	gene_trees_calibrated = foreach( tree=gene_trees_calibrated ) %dopar% 
		store_node_age( tree )
	
	return( gene_trees_calibrated )
}


#' Wilcoxon test of the ortholog conjecture, specifically that contrasts
#' associated with duplication nodes are greater than those for 
#' speciation nodes. A significant result indicates rejection of the null
#' hypothesis that they are not greater.
#' 
#' @param nodes_contrast A tibble of independent contrasts
#' @return Wilcoxon test p value
#' @export
wilcox_oc = function( nodes_contrast ) {
	
	p = wilcox.test( 
		nodes_contrast %>% filter( D=="Y")  %>% .$pic, 
		nodes_contrast %>% filter( D=="N" ) %>% .$pic, 
		alternative="greater" 
	)$p.value
	
	return( p )
}


#' Estimate tau evolution parameters from a tree::treedata tree
#' and add them back into the tree object.
#' 
#' @param nhx A phylogenetic tree and associated tau values 
#' as a treeio::treedata object
#' @param ... Additional arguments to pass to geiger::fitContinuous()
#' @return A phylogenetic trees as a treeio::treedata object, 
#' with estimated brownian motion model parameters stored 
#' in @phylo$model
#' @export
add_model_parameters = function( nhx, ... ) {
	phy = nhx@phylo
	tau_original = nhx@data$Tau [ 1:length( phy$tip.label ) ]
	names( tau_original ) = phy$tip.label

	# Set bounds on fitContinuous per manual, as unbond searches
	# can get stuck. These bounds were selected according to 
	# observed estimates on runs when no fitContinuous calls
	# got stuck. Set ncores=1 so that can wrap in parallel 
	# code without problems.

	brownian_model = fitContinuous( 
		phy, 
		tau_original, 
		model="BM", 
		bounds=list(sigsq=c(0.0, 1.0)),
		ncores=1,
		... 
	)
	nhx@phylo$model_bm$opt = brownian_model$opt

	ou_model = fitContinuous( 
		phy, 
		tau_original, 
		model="OU", 
		bounds=list(sigsq=c(0.0, 1.0), alpha=c(0.0, 3.0)),
		ncores=1,
		... 
	)
	nhx@phylo$model_ou$opt = ou_model$opt	
	
	return( nhx )
}


#' Simulate tau on a tree::treedata. Existing observed values are
#' fit to a brownian model to estimate parameters. These parameters
#' are then used to replace the original values with similated 
#' values.
#' 
#' @param nhx A phylogenetic tree and associated tau values 
#' as a treeio::treedata object
#' @param dup_adjust A multiplier for adjusting the branch lengths 
#' following duplication to effectively change the rate following duplication
#' @param a ancestral state for tau. If not specified, use estimated value 
#' @return A phylogenetic trees and simulated tau values as a 
#' treeio::treedata object
#' @export
sim_tau = function( nhx, dup_adjust=1, a=NA ) {
	
	phy = nhx@phylo
	
	## Modified this part to run the code
	if(class(nhx) == "treedata") 
	{
	  nhx@data <- tidytree::as_tibble(nhx@data)
	}
	
	# Get parameter estimates
	brownian_model = nhx@phylo$model_bm
	
	# Adjust length of branches that descend from duplication events. 
	# This is how we implement heterogeneous rates after duplication 
	# rather than speciation. Extending the branch lengths is equivalent 
	# to using a higher rate along them.
	dup_nodes = nhx@data$node[ nhx@data$D == "Y" ]
	dup_edges = phy$edge[,1] %in% dup_nodes
	phy$edge.length[ dup_edges ] = phy$edge.length[ dup_edges ] * dup_adjust
	
	
	if( is.na( a ) ){
		a = brownian_model$opt$z0
	}
	
	# Simulate trait given the tree and parameter estimates
	x = fastBM( 
		phy, 
		a = a, 
		sig2 = brownian_model$opt$sigsq, 
		bounds = c(-Inf,Inf) 
	) %>% abs()
	
	# Adjust the results so they fall between 0 and 1
	x[ x > 1 ] = 1
	x[ is.nan( x ) ] = 1
	
	
	names( x ) = NULL
	nhx@data$Tau = c( x, rep( NA, nrow( nhx@data ) - length (x ) ) )
	return( nhx )
}


#' Simulate a set of pics from a gamma distribution
#' 
#' @param scale Gamma scale parameter
#' @param n Number of observations per event type	
#' @return A dataframe
#' @export
sim_expectation = function( scale = 1, n = 10000 ){
	expectation = 
		bind_rows(
			data.frame( 
				pic = rgamma( n, 1, scale = scale ),
				Event = "Duplication",
				stringsAsFactors = FALSE
			),
			data.frame( 
				pic = rgamma( n, 1, scale = 1 ),
				Event = "Speciation",
				stringsAsFactors = FALSE
			)
		)
	
	event_levels = c( "Speciation", "Duplication" )
	expectation$Event = 
		factor( expectation$Event, levels=event_levels )
	
	return( expectation )
}

#' Summarize pairwise comparisons between tips of a tree
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object 
#' 	
#' @return A dataframe
#' @export
get_pairwise_summary = function ( nhx ) {
	
	# Create a tibble with one row for each pairwise combination of tip nodes
	ntips = length( nhx@phylo$tip.label )
	D = 1:ntips %>% combn( 2 ) %>% t() %>% as_tibble()
	names( D ) = c( "tip_a", "tip_b" )
	
	# Add a column with the gene tree hash 
	D %<>% mutate( gene_tree = digest( nhx ) )
	
	# Add the tip names
	D %<>% mutate ( name_a = nhx@data$G[ D$tip_a ] )
	D %<>% mutate ( name_b = nhx@data$G[ D$tip_b ] )
	
	# Add the species
	D %<>% mutate ( species_a = nhx@data$species[ D$tip_a ] )
	D %<>% mutate ( species_b = nhx@data$species[ D$tip_b ] )
	
	# Add the values at the tips
	D %<>% mutate ( tau_a = nhx@data$Tau[ D$tip_a ] )
	D %<>% mutate ( tau_b = nhx@data$Tau[ D$tip_b ] )
	D %<>% mutate ( difference = abs( tau_a - tau_b ) )
	
	# For each row, get the most recent common ancestor of the pairwise tip combination
	D$mrca = 
		apply(
			D, 
			1, 
			function( y ) ape::getMRCA( 
				nhx@phylo, 
				c( as.integer( y['tip_a'] ), as.integer( y['tip_b'] ) )
			)
		) %>%
		as.integer()
	
	# Merge data of this mrca to the node
	D %<>% left_join(
		nhx@data %>%
			select( D, Event, node, node_depth, node_age, label, pic ),
		by = c( "mrca" = "node")
		
	)
	
	# Calculate the distance between the tips, which for an ultrametric tree 
	# is twice the age of their most recent common ancestor
	D %<>% mutate( distance = 2 * node_age )
	
	return( D )
	
}

#' Create a ggplot summary of independent contrast values
#' 
#' @param contrasts The contrasts to be displayed 
#' @param mode The type of plot to produce 
#' 	
#' @return A ggplot object
#' @export
make_contrast_plot = function( nodes_contrast, mode="diff" ){
	
	wilcox_test_result = wilcox_oc( nodes_contrast )
	label_p = str_c( "Wilcoxon P=", signif( wilcox_test_result, 2 ) )
	
	if( mode == "freqpoly" ){
		ggcontrasts = nodes_contrast %>%
			ggplot( aes( x=pic, y=..density.., col=Event ) ) + 
				geom_freqpoly( binwidth=.001, position="identity" ) +
				xlab( "Tau Phylogenetic Contrast" ) +
				ylab( "Density" ) +
				xlim( 0, 0.05 ) + 
				ylim( 0, 75 ) +
				theme_classic() + 
				annotate( "text", x = Inf, y = 20, hjust="right", vjust="top", label=label_p ) +
				theme( legend.position="none" )
	} else if( mode == "boxplot" ) {
		ggcontrasts = nodes_contrast %>%
			filter( pic < 0.05 ) %>%
			ggplot( aes(factor(Event), pic ) ) +
			geom_boxplot( notch = TRUE ) +
			coord_flip()
	} else if( mode == "diff" ){
		binwidth = 0.001
		
		D = 
			nodes_contrast %>% 
			group_by(Event) %>% 
			# calculate densities for each group over same range; store in list column
			summarise(d = list(density(pic, from = 0, to = 0.05))) %>% 
			# make a new data.frame from two density objects
			do(data.frame(x = .$d[[1]]$x,    # grab one set of x values (which are the same)
										y = .$d[[1]]$y - .$d[[2]]$y))    # and subtract the y values
		
		D$Event[ D$y >= 0 ] = "Speciation"
		D$Event[ D$y <= 0 ] = "Duplication"
		D$Event = factor( D$Event, levels=c( "Speciation", "Duplication" ) )
			
		ggcontrasts = D %>%
			ggplot(aes(x, y)) +    # now plot
				geom_line() +
				xlab( "Tau Phylogenetic Contrast" ) +
				ylab( "Density Difference" ) +
				xlim( 0, 0.05 ) + 
				ylim( -30, 30 ) +
				geom_area( aes( fill=Event ) ) +
				geom_hline( yintercept = 0 ) +
				annotate( "text", x = Inf, y = -20, hjust="right", vjust="top", label=label_p ) +
				theme_classic() +
				theme(legend.position="none")
		
	}
	
	ggcontrasts
	
}

#' Extends each terminal branch by specified length
#' 
#' @param nhx A phylogenetic tree as a treeio::treedata object
#' @param x Amount to extend each branch by
#' @return A phylogenetic tree as a treeio::treedata object
#' @export
extend_nhx = function( nhx, x ){
  nhx@phylo = hutan::extend_terminal_branches( nhx@phylo, x )
  return( nhx )
}


