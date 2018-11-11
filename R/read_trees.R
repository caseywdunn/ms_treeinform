#' Reads all newick trees in a path and parses with treedata structure so that
#' we can annotate with things like subtree length and height.
#'
#' @param path path to the newick trees
#' @param extension extension for newick trees
#' @param cores number of cores for mclapply
#' @return A list of treedata objects
#' @export
read_trees = function(path, extension, cores) {
  files <- list.files(path,pattern=paste0("*.",extension),full.names=TRUE)
  trees <- mclapply(files, function(x) {
    t <- parse_gene_trees(processTree(x))
    names(t@data) <- c("node","label")
    return(t)
    }, mc.cores=cores)
  return(trees)
}

#' Parses nhx (& newick) text from gene trees and returns branch lengths
#' for duplication events as a treedata object. Also annotates with Siphonophore
#' speciation events.
#'
#' @param tree_text Character string representing nhx tree
#' @return A treedata object with branch lengths as part of data
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

#' Helper function to read in newick trees; will skip over a non-existent
#' tree rather than return an error
#'
#' @param string Path to tree file
#' @return tree_text - Character string representing nhx tree
processTree = function(x) tryCatch({ readLines(x) },
                                   error = function(e) { NULL })