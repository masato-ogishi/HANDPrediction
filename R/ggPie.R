# Pie chart in ggplot.
#' @title Pie chart in ggplot.
#' @param labelSet A vector of categorical variables.
#' @param colorSet The color set to be used.
#' @import ggplot2
#' @export
#' @rdname ggPie
#' @author stack overflow
ggPie <- function (labelSet, colorSet) {
  # prepare name
  deparse(substitute(labelSet)) -> name

  # prepare percents for legend
  table(factor(labelSet)) -> tmp.count1
  prop.table(tmp.count1) * 100 -> tmp.percent1
  paste(round(tmp.percent1, digits=1), " %", sep = "") -> tmp.percent2
  as.vector(tmp.count1) -> tmp.count1

  # find breaks for legend
  rev(tmp.count1) -> tmp.count2
  rev(cumsum(tmp.count2) - (tmp.count2 / 2)) -> tmp.breaks1

  # prepare data
  data.frame(vector1 = tmp.count1, names1 = names(tmp.percent1)) -> tmp.df1

  # plot data
  tmp.graph1 <-
    ggplot(tmp.df1, aes(x = 1, y = vector1, fill = names1)) +
    geom_bar(stat = "identity", color = "black") +
    guides(fill = guide_legend(title = NULL, override.aes = list(colour = NA))) +
    coord_polar(theta = "y") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),
      panel.grid = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(vjust = 1, size = 16, family = "Helvetica", color = "black"),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      legend.text = element_text(size=16, family = "Helvetica", color = "black"),
      plot.title = element_text(hjust = 0.5, vjust = 0.5)
    ) +
    scale_y_continuous(breaks = tmp.breaks1, labels = tmp.percent2) +
    scale_fill_manual(values = colorSet)

  return(tmp.graph1)
}
