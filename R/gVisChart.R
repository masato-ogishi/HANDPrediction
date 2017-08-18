# Comparison of reference and prediction in sankey plot using googleVis API.
#' @title Comparison of reference and prediction in sankey plot using googleVis API.
#' @param reference A set of reference labels.
#' @param prediction A set of predicted labels.
#' @param colors_link A set of colors.
#' @importFrom DescTools Sort
#' @importFrom googleVis gvisSankey
#' @export
#' @rdname gvisChart
gvSankey <- function(reference, prediction, colors_link){
  d <- data.frame(From=reference,
                  To=paste0(prediction, " "), # To avoid label overlapping among From and To sets...
                  Weight=1)
  d <- DescTools::Sort(x=d, ord=c("From","To"), factorsAsCharacter = F)
  colors_array <- paste0("[", paste0("'", colors_link, "'", collapse = ','), "]")
  opts <- paste0("{
                 link: { colorMode: 'gradient',
                 colors: ", colors_array ," },
                 node: { colors: ", colors_array ,",
                 label: { fontName: 'Helvetica',
                 fontSize: 15,
                 color: 'black',
                 bold: true } } }" )
  g <- gvisSankey(d, from="From", to="To", weight="",
                  options=list(sankey=opts))
  plot(g)
}
