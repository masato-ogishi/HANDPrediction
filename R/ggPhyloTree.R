# Phylogenetic tree colored with metadata.
#' @title Phylogenetic tree colored with metadata.
#' @param tree.newick.file A phylotree file name in a Newick format.
#' @param combinedDataset A combined dataframe returned by the dataImportAndCombine function.
#' @param groupLabelName The name of the grouping variable.
#' @param groupVarNames The names of the categories of the grouping variable.
#' @param colorSet The colors to be used.
#' @param groupLabelNames A character vector of the names of the grouping variable.
#' @param groupVarNamesList A list of the names of the categories of the grouping variable.
#' @param colorSetList A list containing the color sets to be used.
#' @param tipped Logical. Whether the output phylogenetic trees should be tipped with corresponding sequence names.
#' @importFrom ape read.tree
#' @importFrom purrr is_empty
#' @importFrom ape drop.tip
#' @importFrom ggtree groupOTU
#' @importFrom ggtree ggtree
#' @importFrom ggtree theme_tree2
#' @importFrom ggtree geom_tiplab
#' @import ggplot2
#' @export
#' @rdname ggPhyloTree
ggPhyloTree <- function(tree.newick.file, combinedDataset, groupLabelName, groupVarNames,
                        colorSet, tipped=F){
  tree <- ape::read.tree(tree.newick.file)
  accession_by_group <- lapply(split(combinedDataset, combinedDataset[[groupLabelName]]), function(d){d[["SequenceID"]]})
  accessions_drop <- setdiff(tree$"tip.label", as.character(unlist(accession_by_group)))
  if(!purrr::is_empty(accessions_drop)) tree <- ape::drop.tip(tree, accessions_drop)
  tree <- ggtree::groupOTU(tree, accession_by_group, group_name=groupLabelName)
  attr(tree, groupLabelName) <- factor(attr(tree, groupLabelName), levels=groupVarNames)
  ptree <- ggtree::ggtree(tree, aes_string(color=groupLabelName), branch.length="none") +
    scale_color_manual(values=colorSet) +
    labs(color=NULL) + guides(color=guide_legend(nrow=2)) +
    ggtree::theme_tree2() +
    theme(legend.position="top")
  if(tipped) ptree <- ptree + ggtree::geom_tiplab(size=1)
  return(ptree)
}

#' @importFrom gridExtra grid.arrange
#' @importFrom ggplot2 scale_x_reverse
#' @import grid
#' @export
#' @rdname ggPhyloTree
ggPhyloTree.Mirror <- function(tree.newick.file, combinedDataset, groupLabelNames, groupVarNamesList,
                               colorSetList, tipped=F){
  ptree.1 <- ggPhyloTree(
    tree.newick.file=tree.newick.file,
    combinedDataset=combinedDataset,
    groupLabelName=groupLabelNames[[1]],
    groupVarNames=groupVarNamesList[[1]],
    colorSet=colorSetList[[1]], tipped=tipped
  )
  ptree.2 <- ggPhyloTree(
    tree.newick.file=tree.newick.file,
    combinedDataset=combinedDataset,
    groupLabelName=groupLabelNames[[2]],
    groupVarNames=groupVarNamesList[[2]],
    colorSet=colorSetList[[2]], tipped=tipped
  ) + ggplot2::scale_x_reverse()
  gridExtra::grid.arrange(ptree.1, ptree.2, nrow=1)
}
