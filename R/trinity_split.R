#' Split up trinity IDs
#' 
#' @param df Data frame
#' @return Df but with trinity column split into trinity and isoform
#' @export
trinity_split <- function(df) {
  trinity <- as.character(df$trinity)
  trinity_str <- strsplit(trinity, "_")
  df$gene <- sapply(trinity_str, function(z) Reduce(function(x,y) paste(x,y,sep="_"), z[2:4]))
  df$isoform <- sapply(trinity_str, function(x) as.numeric(substr(x[[5]], 2, nchar(x))))
  return(select(df,gene,CDS,isoform))
}