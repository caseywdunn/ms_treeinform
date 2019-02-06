#' Replaces values in one column of df1 using df2. If value doesn't
#' exist in df2 then will replace with NA. Primarily intended for use
#' in replacing transcript IDs with gene IDs.
#' 
#' @param df1 Data frame 1, with values to be replaced
#' @param df2 Data frame 2, with replacement values. Should be only 2 columns.
#' @param col_id Column ID of df1 to be replaced.
#' @return Data frame 1 but with values in col_id replaced.
#' @examples
#' df1 <- data.frame(col1=c(1,2), col2=c('a','b'))
#' df2 <- data.frame(ID1=c('a','b'), ID2=c('c','c'))
#' rv <- replace_values(df1, df2, 2)
#' @export
replace_values <- function(df1, df2, col_id) {
  df1[,col_id] <- df2[,2][match(df1[,col_id],df2[,1])]
  return(df1)
}