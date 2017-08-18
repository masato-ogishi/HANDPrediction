# Generate a summary table from the combined data
#' @title Generate a summary table from the combined data
#' @param combinedDataset A combined dataframe returned by the dataImportAndCombine function.
#' @param strataColName A column name representing the variable for stratification.
#' @param varColNamesString A string representing the column names of variables in a summary table.
#' @param factorColNamesString A string representing the column names of categorical (factorized) variables in a summary table.
#' @param nonNormalColNameString A string representing the column names of numerical variables not normally distributed in a summary table.
#' @importFrom dplyr %>%
#' @importFrom stringr str_split
#' @importFrom stringr fixed
#' @importFrom tableone CreateTableOne
#' @export
#' @rdname combinedDatasetSummaryTable
combinedDatasetSummaryTable <- function(
  combinedDataset,
  strataColName="Clinical.Status",
  varColNamesString="Patient.Sex|Georegion|Sample.Tissue.Category|Viral.load|CD4.count|Number.of.patient.seqs",
  factorColNamesString="Patient.Sex|Georegion|Sample.Tissue.Category",
  nonNormalColNameString="Viral.load|CD4.count|Number.of.patient.seqs"
){
  varColNames <- str_split(varColNamesString, fixed("|")) %>% unlist()
  factorColNames <- str_split(factorColNamesString, fixed("|")) %>% unlist()
  nonNormalColNames <- str_split(nonNormalColNameString, fixed("|")) %>% unlist()
  tableOne <- suppressMessages(suppressWarnings(
    CreateTableOne(data=combinedDataset,
                   strata=strataColName,
                   vars=varColNames,
                   factorVars=factorColNames,
                   includeNA=T)
  ))
  tableOne <- as.data.frame(print(tableOne, nonnormal=nonNormalColNames))
  return(tableOne)
}

# Export the summary table in a Word format
#' @title Export the summary table in a Word format
#' @param tableOneDF A "tableone" object returned by the combinedDatasetSummaryTable function.
#' @param outputFileName A file name for the output Word file.
#' @importFrom dplyr %>%
#' @importFrom ReporteRs docx
#' @importFrom ReporteRs addFlexTable
#' @importFrom ReporteRs FlexTable
#' @importFrom ReporteRs cellProperties
#' @importFrom ReporteRs textBold
#' @importFrom ReporteRs writeDoc
#' @export
#' @rdname combinedDatasetSummaryTable
combinedDatasetSummaryTable.Export <- function(tableOneDF, outputFileName="tableOne.docx"){
  docx( ) %>%
    addFlexTable(tableOneDF %>%
                   FlexTable(header.cell.props=cellProperties(background.color="#003366"),
                             header.text.props=textBold(color="white"),
                             add.rownames=T)) %>%
    writeDoc(file=outputFileName)
  ## Adopted from: http://www.r-bloggers.com/table-1-and-the-characteristics-of-study-population/
}
