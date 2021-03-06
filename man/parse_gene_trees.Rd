% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_trees.R
\name{parse_gene_trees}
\alias{parse_gene_trees}
\title{Parses nhx (& newick) text from gene trees and returns branch lengths
for duplication events as a treedata object. Also annotates with Siphonophore
speciation events.}
\usage{
parse_gene_trees(tree_text)
}
\arguments{
\item{tree_text}{Character string representing nhx tree}
}
\value{
A treedata object with branch lengths as part of data
}
\description{
Parses nhx (& newick) text from gene trees and returns branch lengths
for duplication events as a treedata object. Also annotates with Siphonophore
speciation events.
}
