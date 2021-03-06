% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggPhyloTree.R
\name{ggPhyloTree}
\alias{ggPhyloTree}
\alias{ggPhyloTree.Mirror}
\title{Phylogenetic tree colored with metadata.}
\usage{
ggPhyloTree(
  tree.newick.file,
  combinedDataset,
  groupLabelName,
  groupVarNames,
  colorSet,
  tipped = F
)

ggPhyloTree.Mirror(
  tree.newick.file,
  combinedDataset,
  groupLabelNames,
  groupVarNamesList,
  colorSetList,
  tipped = F
)
}
\arguments{
\item{tree.newick.file}{A phylotree file name in a Newick format.}

\item{combinedDataset}{A combined dataframe returned by the dataImportAndCombine function.}

\item{groupLabelName}{The name of the grouping variable.}

\item{groupVarNames}{The names of the categories of the grouping variable.}

\item{colorSet}{The colors to be used.}

\item{tipped}{Logical. Whether the output phylogenetic trees should be tipped with corresponding sequence names.}

\item{groupLabelNames}{A character vector of the names of the grouping variable.}

\item{groupVarNamesList}{A list of the names of the categories of the grouping variable.}

\item{colorSetList}{A list containing the color sets to be used.}
}
\description{
Phylogenetic tree colored with metadata.
}
