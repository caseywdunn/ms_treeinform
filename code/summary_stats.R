#' Return only species labels for each tip in a tree
#'
#' @param phylo phylo object
#' @return vector of species labels
model_ids = function(phylo) {
  sapply(strsplit(phylo$tip.label, "@"), function(x) x[1])
}

#' Count number of transcripts for each species in a tree
#'
#' @param tree Phyldog-annotated nhx tree
#' @return list of model ID and species taxa tuple
count_transcripts = function(tree) {
  phylo_labels = model_ids(tree)
  table(phylo_labels)
}

#' Count number of transcripts for each species in a tree for a list of trees
#'
#' @param trees list of Phyldog-annotated nhx objects
#' @return list of count of transcripts for each species for each tree
multi_count_transcripts = function(trees, cores) {
  mclapply(trees, count_transcripts, mc.cores=cores)
}

#' Count number of transcripts for each species summed across all trees
#'
#' @param phylos list of phylo objects
#' @param cores number of cores for parallel
#' @return list of count of transcripts for each species summed across all trees
count_all_transcripts = function(phylos, cores) {
  all_transcripts = unlist(mclapply(phylos, model_ids, mc.cores=cores))
  table(all_transcripts)
}
