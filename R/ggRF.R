# Visualize Random Forest model in ggplot
#' @title Phylogenetic tree colored with metadata.
#' @param caret_rf_model The final model of the random forest model from caret
#' @param num_nodes The number of nodes for visualization. Either "smallest", "largest", "median", or an integer.
#' @param brewerPalName The color palette name in the Brewer palettes.
#' @param lastColorGrey Logical. Whether the last category should be colored as grey.
#' @importFrom randomForest getTree
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom igraph V
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph delete_vertices
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggraph ggraph
#' @importFrom ggraph geom_edge_link
#' @importFrom ggraph geom_node_point
#' @importFrom ggraph geom_node_text
#' @importFrom ggraph geom_node_label
#' @importFrom grid unit
#' @import ggplot2
#' @export
#' @rdname ggRF
#' @author Shirin Glander
ggRF <- function(caret_rf_model, num_nodes="median",
                 brewerPalName="Set1", lastColorGrey=F){

  # get tree by index
  tree_num <- num_nodes
  tree_nd_set <- caret_rf_model$forest$ndbigtree
  if(num_nodes=="smallest")
   tree_num <- which(tree_nd_set==min(tree_nd_set))[[1]]
  if(num_nodes=="largest")
    tree_num <- which(tree_nd_set==max(tree_nd_set))[[1]]
  if(num_nodes=="median")
    tree_num <- which(tree_nd_set==median(tree_nd_set))[[1]]
  tree <- randomForest::getTree(caret_rf_model,
                                k = tree_num,
                                labelVar = TRUE) %>%
    tibble::rownames_to_column() %>%
    # make leaf split points to NA, so the 0s won't get plotted
    mutate(`split point` = ifelse(is.na(prediction), `split point`, NA))

  # prepare data frame for graph
  graph_frame <- data.frame(from = rep(tree$rowname, 2),
                            to = c(tree$`left daughter`, tree$`right daughter`))

  # convert to graph and delete the last node that we don't want to plot
  graph <- igraph::graph_from_data_frame(graph_frame) %>% igraph::delete_vertices("0")

  # set node labels
  igraph::V(graph)$node_label <- gsub("_", " ", as.character(tree$`split var`))
  igraph::V(graph)$leaf_label <- as.character(tree$prediction)
  igraph::V(graph)$split <- as.character(round(tree$`split point`, digits = 2))

  # plot
  col_lab <- unique(as.character(tree$prediction))
  col_lab <- col_lab[!is.na(col_lab)]
  cols <- ggsci::pal_npg()(length(col_lab))
  if(lastColorGrey) cols[length(cols)] <- "grey25"
  plot <- ggraph(graph, 'dendrogram') +
    theme_bw() +
    geom_edge_link() +
    geom_node_point() +
    geom_node_text(aes(label = node_label), repel = T) +
    geom_node_label(aes(label = split), vjust = 2.5, fill = "white", repel = F) +
    geom_node_label(aes(label = leaf_label, fill = leaf_label),
                    repel = F, colour = "white", fontface = "bold", show.legend = F) +
    scale_fill_manual(values = cols) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          panel.background = element_blank(),
          plot.background = element_rect(fill = "white"),
          panel.border = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 18))
  print(plot)
}
## This function was originally posted on the blog by Dr Shirin Glander.
## https://shiring.github.io/machine_learning/2017/03/16/rf_plot_ggraph
