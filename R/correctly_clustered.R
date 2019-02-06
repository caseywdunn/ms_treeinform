#' Computes number of pairs in a matrix from `table` that
#' are in the same subset. Intended to identify transcripts
#' that have been mapped to the same gene clusters between
#' either two different methods or one method and the true clustering.
#' 
#' Note: gene clusters consisting of a single transcript will not be identified
#' in this way. Use `identify_singletons` in combination with this.
#' 
#' @param tab Object from `table`.
#' @return Number of pairs that are in the same subset.
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' tab <- table(df$col1, df$col2)
#' ss <- same_subset(tab)
#' @export
same_subset <- function(tab) {
  return(sum(choose(tab,2)))
}

#' Computes number of pairs in a matrix from `table` that
#' are in the different subsets.
#' 
#' @param tab Object from `table`.
#' @return Number of pairs that are in different subsets.
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' tab <- table(df$col1, df$col2)
#' ds <- different_subset(tab)
#' @export
different_subset <- function(tab) {
  return(sum(choose(rowSums(tab),2))-a)
}