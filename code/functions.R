#' Computes pdf of duplication event for time s given tree origin time t,
#' birth rate lambda, death rate mu
#'
#' @param s time of duplication event
#' @param t time of tree origin
#' @param lambda birth rate in duplication events/unit of time
#' @param mu death rate in loss events/unit of time
dup_pdf <- function(s,t,lambda,mu) {
  #if(s > t) { return(0) }
  #else {
  num <- (lambda-mu)^2*exp(-(lambda-mu)*s)*(lambda - mu*exp(-(lambda-mu)*t))
  denom <- (lambda-mu*exp(-(lambda-mu)*s))^2*(1-exp(-(lambda-mu)*t))
  f <- num/denom
  return(f)
  #}
}

#' Cdf of duplication event for time s given tree origin time t,
#' birth rate lambda, death rate mu - uniform (for sampling)
#'
#' @param s time of duplication event
#' @param t time of tree origin
#' @param lambda birth rate in duplication events/unit of time
#' @param u sample from uniform distribution
#' @param mu death rate in loss events/unit of time
dup_cdf_sampler <- function(s,t,lambda,mu,u) {
  if(s > t) { return(1-u) }
  else {
  num <- (1-exp(-(lambda-mu)*s))*(lambda - mu*exp(-(lambda-mu)*t))
  denom <- (lambda-mu*exp(-(lambda-mu)*s))*(1-exp(-(lambda-mu)*t))
  f <- num/denom
  return(f-u)
  }
}

dup_cdf <- function(s,t,lambda,mu) {
    num <- (1-exp(-(lambda-mu)*s))*(lambda - mu*exp(-(lambda-mu)*t))
    denom <- (lambda-mu*exp(-(lambda-mu)*s))*(1-exp(-(lambda-mu)*t))
    f <- num/denom
    return(f)
}

#' Helper function to read in newick trees; will skip over a non-existent
#' tree rather than return an error
#'
#' @param string Path to tree file
#' @return tree_text - Character string representing nhx tree
#' @export
processTree = function(x) tryCatch({ readLines(x) },
                                   error = function(e) { NULL })

#' Annotates a tree from the agalma 7 taxa Siphonophora test dataset
#' with internal node clade names
#'
#' @param phylo Phylo object containing only tips labeled with the 7
#' taxa from the test dataset
#' @return phylo Same phylo object, now with internal node annotations
annotate_siph = function( phylo ) {
  # get vector of internal node numbers
  n_tips = length(phylo$tip.label)
  internal_nodes = n_tips+1:phylo$Nnode

  # go thru vector & assign clade label based on
  # set of descendant tips
  clade_labels = sapply(internal_nodes, function(x) suppressWarnings(siph_help(phylo, x)))

  # go thru internal nodes and make sure there are no assignments of
  # speciation events with ancestor/child nodes as the same speciation event
  # (i.e. (Cnidaria, Cnidaria):Cnidaria)
  new_clade_labels = sapply(internal_nodes, function(x) remove_ancestor(phylo, x, clade_labels, internal_nodes))

  phylo$node.label = new_clade_labels

  return(phylo)
}

#' Checks if for 2 transcripts either both the trinity clustering and corset clustering
#' are both the same or if they are both different.
#'
#' @param row_ids transcript IDs
#' @param corset_clustering data frame with trinity & corset gene clusterings
#' @return 1 if trinity & corset clustering are both the same/different, 0 else
check_cluster_identity = function(row_ids, clustering) {
  rows <- clustering[row_ids,]
  unique_rows <- length(unique(rows$Trinity.gene)) + length(unique(rows$Corset.gene))
  # if trinity_gene & corset_gene are both same, unique_rows = 2
  # if trinity_gene & corset_gene are both different, unique_rows = 4
  # else = 3
  if(unique_rows == 2 || unique_rows == 4) return(1)
  else return(0)
}

#' Computes the proportion of transcript pairs that share the same gene assignment
#' between trinity and corset. A pair of transcripts shares the same gene assignment
#' if they are assigned to the same gene by both trinity & corset, or if they are
#' assigned to different genes by both trinity & corset.
#'
#' @param corset_clustering data frame with trinity & corset gene clusterings; from corset output
#' @return proportion of transcript pairs sharing same gene assignment
compute_corset_trinity_similarity = function(corset_clustering) {
  all_pairs <- combn(1:nrow(clustering),2,simplify=FALSE)
  gene_assignments <- mclapply(all_pairs,function(x) check_cluster_identity(x, clustering))
  return(sum(gene_assignments)/length(all_pairs))
}

#' Returns cluster size distribution
#'
#' @param clustering vector with gene clusterings
#' @return data frame with columns size (=cluster size) and freq (=frequency of that cluster)
cluster_size_distribution = function(clustering) {
  df<-data.frame(table(table(clustering)))
  colnames(df) <- c("size","freq")
  return(df)
}

#' Goes through corset clustering file and filters out all singleton clusters
#'
#' @param corset_clustering data frame with triniy & corset gene clusterings
#' @return row numbers of transcripts that have different corset and trinity clusterings
filter_singletons = function(corset_clustering) {
  n_occur1 <- data.frame(table(corset_clustering$Trinity.gene))
  n_occur2 <- data.frame(table(corset_clustering$Corset.gene))
  filter1 <- corset_clustering[corset_clustering$Trinity.gene %in% n_occur1$Var1[n_occur1$Freq>1],]
  filter2 <- corset_clustering[corset_clustering$Corset.gene %in% n_occur2$Var1[n_occur2$Freq>1],]
  corset_clustering[rowSums(corset_clustering)>1,] #think this should do the trick but need to test
  # actually needs slightly more thought because needs to be true for both corset & trinity columns
}

#' Go through internal nodes and make sure there is no assignment of
#' speciation events with ancestor/child nodes as the same speciation event
#' (i.e. (Cnidaria, Cnidaria):Cnidaria)
#'
#' @param phylo Phylo object containing only tips labeled with the 7
#' taxa from the test dataset
#' @param x node number corresponding to internal speciation node in phylo
#' @param clade_labels Vector of clade labels
#' @param internal_nodes Vector of internal node
#' @return annotation either NA or the original clade label
remove_ancestor = function(phylo, x, clade_labels, internal_nodes) {
  # get descendent internal nodes
  descendants = hutan::descendants( phylo, x )
  descendants = descendants[ descendants %in% internal_nodes ]
  n_tips = length(phylo$tip.label)

  descendant_names = clade_labels[descendants-n_tips]

  if( clade_labels[x-n_tips] %in% descendant_names ){ return(NA) }
  else { return(clade_labels[x-n_tips]) }
}

#' Helper function for annotate_siph
#'
#' @param phylo Phylo object containing only tips labeled with the 7
#' taxa from the test dataset
#' @param x node number corresponding to internal node in phylo
#' (note: we do not check that it is an internal node)
#' @return clade_label Clade annotation for node x
siph_help = function(phylo, x) {
  # get unique tip descendants
  descendants = sort(unique(sapply(strsplit(names(hutan::tip_descendants(phylo, x)), '_'), '[[', 1)))
  if (all(descendants == c("Abylopsis", "Craseoa"))) { return("Calycophora") }
  if (all(descendants == c("Agalma", "Nanomia"))) { return("Agalmatidae") }
  if (all(descendants == c("Abylopsis", "Agalma", "Craseoa", "Nanomia"))) { return("Chodorophora") }
  if (all(descendants == c("Abylopsis", "Agalma", "Craseoa", "Nanomia", "Physalia"))) { return("Siphonophora") }
  if (all(descendants == c("Abylopsis", "Agalma", "Craseoa", "Hydra", "Nanomia", "Physalia"))) { return("Hydrozoa") }
  if (all(descendants == c("Abylopsis", "Agalma", "Craseoa", "Hydra", "Nanomia", "Nematostella", "Physalia"))) { return("Cnidaria") }
  else { return(NA) }
}

#' Parses nhx text from gene trees and returns branch lengths
#' for duplication events
#'
#' @param tree_text Character string representing nhx tree
#' @return A list of branch lengths
#' @export
parse_gene_trees = function( tree_text ){
  if (is.null(tree_text)) { return(NULL) }
  tree_tc = textConnection( tree_text )
  tree = treeio::read.nhx( tree_tc )
  close( tree_tc )

  tree@phylo = annotate_siph(tree@phylo)
  tree@data$label = c( tree@phylo$tip.label, tree@phylo$node.label )

  n_nodes = nrow( tree@data )
  n_tips = length( tree@phylo$tip.label )
  internal_nodes = ( n_tips + 1 ):n_nodes
  is_speciation = tree@data$Ev == "S"
  is_speciation[ is.na( is_speciation ) ] = FALSE
  internal_speciation_nodes = tree@data$node[ ( tree@data$node > n_tips ) & is_speciation ]

  return( tree )
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

#' Adjust branch lengths of a phylogenetic tree to make it ultrametric and
#' to time calibrate the speciation nodes to fixed values.
#'
#' @param nhx A phylogenetic tree as a treeio::treedata object, with speciation
#' nodes annotated with a value of "D" in @data$Ev
#' @param calibration_times A dataframe with two columns: age is the age
#' of a clade, and clade is the name of the clade
#' @param ... Any additional arguments to pass to ape::chronos()
#' @return A treeio::treedata object if successfully calibrated
#' @export
calibrate_tree = function ( nhx, calibration_times, ... ) {

  # Create calibration matrix for speciation nodes
  calibration =
    nhx@data[ !is.tip.nhx( nhx ), ] %>%
    filter(Ev=="S" & !is.na(label)) %>%
    #filter( Ev == "S" ) %>%
    left_join( calibration_times, c( "label" = "clade" ) ) %>%
    mutate( age.min = age ) %>%
    mutate( age.max = age ) %>%
    select( node, age.min, age.max ) %>%
    mutate( soft.bounds = NA )

  tree = try(
    ape::chronos( nhx@phylo, calibration=calibration)
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

#' Calibrate a list of treeio::treedata objects
#'
#' @param trees A list of treeio::treedata objects with branch lengths
#'  in expected numbers of substitutions
#' @param cores Number of cores to use for mclapply
#' @return A list of calibrated treeio::treedata objects with
#'  branch lengths in units of time
calibrate_trees = function(trees, cores) {
  return(mclapply(trees, function(x) calibrate_tree(x, calibration_times), mc.cores=cores))
}

#' Annotate a list of treeio::treedata objects with subtree heights
#'
#' @param trees A list of treeio::treedata objects
#' @param cores Number of cores to use for mclapply
#' @return A list of calibrated treeio::treedata objects with
#'  annotated heights
heights_trees = function(trees, cores) {
  mclapply(trees, function(x) {
    if(is.na(x)) { return(NA) }
    else { return(heights(x)) }
  }, mc.cores=cores)
}

#' Compute inner node heights from a nhx object where @data contains a column corresponding
#' to duplication events.
#'
#' @param nhx A phylogenetic tree as a treeio::treedata object, with speciation
#' nodes annotated with a value of "N" in @data$D
#' @return A phylogenetic tree with subtree heights annotated
#' @export
heights = function(nhx) {
  phylo = nhx@phylo
  h = hutan::distance_from_tip(phylo)
  nhx@data[,"height"] = h
  return(nhx)
}

#' Calibrates trees and extracts duplication times from phyldog gene trees.
#'
#' @param k indices of phyldog gene trees to be extracted
#' @param ResultFiles path to phyldog ResultFiles directory
#' @param calibration_times Calibration times for the gene trees
#' @return A vector of duplication times
dt_phyldog = function(k, ResultFiles, calibration_times) {
  cores = parallel::detectCores()
  if (cores < 1) { cores = 1 }

  trees <- mclapply(k, function(x) parse_gene_trees(processTree(paste0(ResultFiles, x, ".ReconciledTree"))), mc.cores = cores)
  trees <- trees[which(!unlist(mclapply(trees, is.null)))]

  calibrated = calibrate_trees(trees, cores)
  annotated = heights_trees(calibrated, cores)
  no_na = annotated[which(!is.na(annotated))]

  # filter out trees with height > 1
  no_heights = mclapply(no_na, function(x) {
    h = x@data %>% select(height) %>% max
    if(h > 1) { return(NA) }
    else { return(x) }}, mc.cores=cores)
  no_heights = no_heights[which(!is.na(no_heights))]
  dt = unlist(mclapply(no_heights, function(x) unlist(x@data %>% filter(Ev=="D") %>% select(height))))
  dt = dt[which(!is.na(dt))]

  return(list(trees, annotated, dt))
}
