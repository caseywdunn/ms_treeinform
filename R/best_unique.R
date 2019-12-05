#' Best unique value from a column in a (BLAST) dataframe.
#' 
#' Takes the unique result from column col with the lowest evalue then highest percent_identity. If there are still
#' multiple results from column col after the first one is selected.
#' Intended to select most significant BLAST results.
#' 
#' @param df Data frame
#' @param col column to take unique results from
#' @return Data frame filtered to only unique values of A with lowest evalue & highest % identity
#' @examples
#' df <- data.frame(col1=c('a','a', 'a'), evalue=c(1,2,1), percent_identity=c(94,99,94))
#' df2 <- best_unique(df,'col1')
#' @importFrom dplyr group_by filter
#' @export
best_unique <- function(df, col) {
  df_filter <- df %>% group_by_(col) %>% filter(evalue==min(evalue)) %>%
    filter(percent_identity==max(percent_identity))
  return(df_filter[!duplicated(df_filter[,col]),])
}