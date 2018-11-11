#' Creates a single data frame from a list of cluster size distribution data frames
#' with Method column indicating method (examples: Trinity, Corset, treeinform) and
#' ID column indicating sample (examples: SRX5, melanogaster).
#' 
#' @param list_cluster_dfs List of cluster size distribution data frames
#' @param methods Vector of methods. Must be same length as list_cluster_dfs
#' @param ids Vector of IDs. Must be same length as list_cluster_dfs
#' @return Data frame with columns: size, freq, Method, ID
#' @importFrom assertthat assert_that
#' @examples
#' df1 <- data.frame(size=c(1,2), freq=c(10,20))
#' df2 <- data.frame(size=c(1,2), freq=c(5,5))
#' methods <- c("method1", "method2")
#' ids <- c("id1", "id1")
#' single_cluster_df(list(df1,df2),methods,ids)
#' @export
single_cluster_df <- function(list_cluster_dfs, methods, ids) {
  l <- length(list_cluster_dfs)
  assert_that(length(methods) == l)
  assert_that(length(ids) == l)
  new_dfs <- lapply(1:l, function(x) {
    list_cluster_dfs[[x]]$Method<-methods[[x]]
    list_cluster_dfs[[x]]$ID<-ids[[x]]
    return(list_cluster_dfs[[x]])
    })
  df<-do.call(rbind,new_dfs)
  return(df)
}

#' Plots cluster sizes and counts. Plot does not come with title; user has to add their own.
#' User is also recommended to adjust breaks with scale_x_discrete as well.
#' 
#' @param df Data frame with columns: size, freq, Method, ID
#' @export
plot_clusters <- function(df) {
  ggplot(data=df, aes(x=size, y=freq, color=Method, group=Method)) +
    geom_point(alpha=0.7) +
    scale_y_log10() +
    xlab("Cluster Size") +
    ylab("Count") +
    theme_classic() +
    geom_line() +
    facet_wrap(~ID, ncol=2)
}