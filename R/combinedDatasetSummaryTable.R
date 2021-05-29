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

