# Stepwise feature reduction.
#' @title Stepwise feature reduction.
#' @param featureSet A set of features as a starting point. Can be either a character vector or a list.
#' @param ml.data A machine-learning-ready dataframe returned by the machineLearning.DataFormat function.
#' @param outcomeLabelName A column name for neurological outcomes.
#' @param ... Parameters for the machineLearning.Comparison.Bootstrap function.
#' @importFrom dplyr %>%
#' @importFrom dplyr select
#' @importFrom dplyr matches
#' @importFrom dplyr mutate
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom DescTools Sort
#' @export
machineLearning.StepwiseFeatureReduction <- function(featureSet, ml.data, outcomeLabelName, ...){
  # Stepwise feature reduction
  vars <- featureSet
  statDFList <- list()
  statSummaryDFList <- list()
  remove_NA <- function(l){l[!is.na(l)]}
  for(j in 1:(length(featureSet)-1)){
    statDFList[[j]] <- lapply(1:length(vars), function(i){
      ml <- try(machineLearning.Comparison.Bootstrap(
        ml.data=dplyr::select(ml.data, dplyr::matches(paste0(c(outcomeLabelName, vars[-i]), collapse="|"))),
        outcomeLabelName=outcomeLabelName,
        ...
      ), silent=T)
      if(class(ml)!="try-error"){
        ml$"MLStatDF" %>% dplyr::mutate("RemovedFeature"=vars[[i]])
      }else{
        NA
      }
    }) %>% remove_NA() %>% dplyr::bind_rows() %>% dplyr::mutate("StepwiseID"=j)
    if(nrow(statDFList[[j]])>0){
      statSummaryDFList[[j]] <- statDFList[[j]] %>%
        dplyr::select(StepwiseID, Algorithm, RemovedFeature, Accuracy) %>%
        dplyr::group_by(StepwiseID, Algorithm, RemovedFeature) %>%
        dplyr::summarise(MeanAccuracy=suppressWarnings(mean(Accuracy)))
      vars <- statSummaryDFList[[j]] %>%
        dplyr::filter(Algorithm=="Stack")
      vars <- DescTools::Sort(x=vars, ord="MeanAccuracy", decreasing=T)[-1,]
      vars <- as.character(vars[["RemovedFeature"]])
    }else{
      break
    }
  }
  statDF <- dplyr::bind_rows(statDFList)
  statSummaryDF <- dplyr::bind_rows(statSummaryDFList)

  # Baseline
  statDF.baseline <- try(machineLearning.Comparison.Bootstrap(
    ml.data=dplyr::select(ml.data, dplyr::matches(paste0(c(outcomeLabelName, featureSet), collapse="|"))),
    outcomeLabelName=outcomeLabelName,
    ...
  ), silent=T)
  if(class(statDF.baseline)!="try-error"){
    statDF.baseline <- statDF.baseline$"MLStatDF" %>% dplyr::mutate("RemovedFeature"=NA, "StepwiseID"=0)
    statSummaryDF.baseline <- statDF.baseline %>%
      dplyr::select(StepwiseID, Algorithm, RemovedFeature, Accuracy) %>%
      dplyr::group_by(StepwiseID, Algorithm, RemovedFeature) %>%
      dplyr::summarise(MeanAccuracy=suppressWarnings(mean(Accuracy)))
    statDF <- dplyr::bind_rows(statDF.baseline, statDF)
    statSummaryDF <- dplyr::bind_rows(statSummaryDF.baseline, statSummaryDF)
  }

  # Extract the most important variables
  bestStepID <- statDF %>%
    dplyr::filter(Algorithm=="Stack") %>%
    DescTools::Sort(ord=c("Accuracy", "StepwiseID"), decreasing=T) %>%
    (function(d){d[["StepwiseID"]][1]})+1
  bestVars <- dplyr::filter(statSummaryDF, Algorithm=="Stack", StepwiseID==bestStepID)[["RemovedFeature"]] %>% unique()

  return(list("StepwiseStatDF"=statDF, "StepwiseStatSummaryDF"=statSummaryDF, "MostImportantFeatures"=bestVars))
}
