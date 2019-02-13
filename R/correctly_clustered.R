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
#' are in the same subset in X but different subsets in Y.
#' 
#' @param tab Object from `table`.
#' @return Number of pairs that are in the same subset in X & different in Y.
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' tab <- table(df$col1, df$col2)
#' ds <- same_X_diff_Y(tab)
#' @export
same_X_diff_Y <- function(tab) {
  return(sum(choose(rowSums(tab),2))-same_subset(tab))
}

#' Computes number of pairs in a matrix from `table` that
#' are in the same subset in Y but different subsets in X.
#' 
#' @param tab Object from `table`.
#' @return Number of pairs that are in the same subset in X & different in Y.
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' tab <- table(df$col1, df$col2)
#' ds <- same_Y_diff_X(tab)
#' @export
same_Y_diff_X <- function(tab) {
  return(sum(choose(colSums(tab),2))-same_subset(tab))
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
  return(choose(sum(tab),2)-sum(choose(colSums(tab),2))-sum(choose(rowSums(tab),2))+sum(choose(tab,2)))
}

#' Given a data frame with two columns consisting of clusterings,
#' gives number of entries that have been correctly clustered, assuming
#' that the 2nd column is the true clustering.
#' 
#' @param df Data frame with 2 columns of clusterings
#' @return Number of entries that have been correctly clustered
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' cc <- correctly_clustered(df)
#' @export
correctly_clustered <- function(df) {
  tab <- table(df[,1],df[,2])
  return(same_subset(tab) + singletons(df))
}

#' Computes number of singleton clusters.
#' 
#' @param df Data frame with two columns containing only clusterings (do not need same label)
#' @return row numbers of transcripts that have singleton clusters and the number of singleton clusterings
#' @examples
#' df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
#' s <- singletons(df)
#' @export
singletons <- function(df) {
  c1 <- df[,1]
  c2 <- df[,2]
  d <- duplicated(c1) | duplicated(c1, fromLast=TRUE) | duplicated(c2) | duplicated(c2, fromLast=TRUE)
  list(df[d,],sum(d))
}

#' Computes Adjusted Rand Index taking into account count of filtered out
#' singletons. Most of source code directly lifted from mclust package.
#'
#' @param f number of filtered out singletons
#' @param c corset & trinity filtered clustering data frame
#' @return adjusted Rand Index
ari <- function(f, c) {
  x<-as.vector(c$Trinity.gene)
  y<-as.vector(c$Corset.gene)
  if(length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x,y)
  if(all(dim(tab)==c(1,1))) return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab)+f, 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d)) /
    ((a + b + a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}