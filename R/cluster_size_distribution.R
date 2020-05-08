#' Returns cluster size distribution
#'
#' @param clustering vector with gene names/IDs
#' @return data frame with columns size (=cluster size) and freq (=frequency of that cluster)
#' 
#' @example 
#' # gene vector is based off of Trinity annotation but can be any format as long as it is gene names/IDs
#' # for example Corset annotation would be of format Cluster-1.0, Cluster-2.0, etc.
#' genes <- c("DN1_c0_g1","DN1_c0_g1","DN1_c0_g2","DN2_c0_g1")
#' cluster_size_distribution(genes)
#' # output is data frame with 2 clusters of size 1 and 1 cluster of size 2
cluster_size_distribution = function(clustering) {
  df<-data.frame(table(table(clustering)))
  colnames(df) <- c("size","freq")
  df$size <- as.numeric(size)
  return(df)
}

#' Returns number of pairs there should be if all transcripts were correctly assigned.
#' mutate(test_csd, numpairs=choose(size, 2)*freq) %>% summarise(sum(numpairs))