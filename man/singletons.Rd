% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/correctly_clustered.R
\name{singletons}
\alias{singletons}
\title{Computes singleton clusters.}
\usage{
singletons(df)
}
\arguments{
\item{df}{Data frame with two columns containing only clusterings (do not need same label)}
}
\value{
vector of TRUE/FALSE for transcripts that have singleton clusters
}
\description{
Computes which rows in the data frame have singleton clusters for both column 1 and column 2.
}
\examples{
df <- data.frame(col1=c(1,1,2,2,2,2,3,4,5),col2=c(1,1,2,3,3,3,4,4,5))
s <- singletons(df)
}
