# Predict HAND status using pre-trained classifiers.
#' @title Predict HAND status using pre-trained classifiers.
#' @param X Either a dataframe or a set of sequences.
#' @param ... Parameters passed to an appropriate working function.
#' @param ml.data.predict A machine-learning-ready dataframe,
#' @param sequences A set of amino acid sequences of HIV C2V3C3 region.
#' @param gapPositions Alignment gap positions.
#' @param analysisType An indicator for the type of  analysis. Currently, "boolian" and "AAIndex" are implemented.
#' @param aaIndexType Currently, "all" and "physicochemical" are implemented. If set to be "all", a total of 544 AAIndices included in the BioSeqClass package are used. If set to be "physicochemical", pre-selected 78 AAIndices are used. See the original paper for further details.
#' @param normalized Logical. Whether AAIndex values are rescaled to [0,1] range.
#' @param ml.metadata A metadata.
#' @param colNames.from Names of the columns in the metadata to be incorporated in the subsequent analysis.
#' @param colNames.to (Optional) Column names of the included metadata. Default values are the same as colNames.from.
#' @param outcomeLabelName A column name representing the neurological outcomes.
#' @param sourceLabelName A column name representing the sample sources.
#' @param pp An indicator for preprocessing. If set "internal", variables are first screened using the zero-variance and near-zero-variance methods implemented in the caret package, and then further screened using the correlation method implemented in the caret package. If set "external", the matrix is pre-processed using externally provided functions. Omit preprocessing by setting "none".
#' @param pp.by.sequence (Optional) Users can provide their own pre-processing function (per sequence).
#' @param pp.by.patient (Optional) Users can provide their own pre-processing function (per patient).
#' @param modelList A list of pre-trained models.
#' @param model.stack A pre-trained stacking model.
#' @export
#' @rdname machineLearning.Prediction
machineLearning.Stacking.Predict <- function(X, ...){
  if(is.data.frame(X)) return(machineLearning.Stacking.Predict_From_PreprocessedDataset(ml.data.predict=X, ...))
  if(is.character(unlist(X))) return(machineLearning.Stacking.Predict_From_Metadataset(sequences=X, ...))
}

#' @export
#' @rdname machineLearning.Prediction
machineLearning.Stacking.Predict_From_PreprocessedDataset <- function(
  ml.data.predict, modelList, model.stack
){
  algorithmLabelList <- names(modelList)
  predProbList <- list(1:length(algorithmLabelList))
  for(i in 1:length(algorithmLabelList)){
    predProbList[[i]] <- predict(modelList[[i]], newdata=ml.data.predict, type="prob")
  }
  names(predProbList) <- algorithmLabelList
  predProbDF <- as.data.frame(predProbList)
  predOutcomeDF <- as.data.frame(predict(model.stack, predProbDF, type="prob"))
  colnames(predOutcomeDF) <- paste0("Stack.", colnames(predOutcomeDF))
  predOutcomeDF <- data.frame("PredictedOutcome"=predict(model.stack, predProbDF, type="raw"),
                              predProbDF,
                              predOutcomeDF)
  return(predOutcomeDF)
}

#' @export
#' @rdname machineLearning.Prediction
machineLearning.Stacking.Predict_From_Metadataset <- function(
  sequences, gapPositions=c(125),
  analysisType="AAIndex", aaIndexType="physicochemical", normalized=T,
  ml.metadata,
  colNames.from, colNames.to,
  outcomeLabelName, sourceLabelName,
  pp="external", pp.by.sequence, pp.by.patient,
  modelList, model.stack
){
  ML_Matrix_Pred <- machineLearning.Matrix(
    sequences=sequences, gapPositions=gapPositions,
    analysisType=analysisType, aaIndexType=aaIndexType, normalized=normalized
  )
  ML_Data_Pred <- machineLearning.DataFormat(
    ml.metadata=ml.metadata, ml.matrix=ML_Matrix_Pred,
    colNames.from=colNames.from, colNames.to=colNames.to,
    outcomeLabelName=outcomeLabelName, sourceLabelName=sourceLabelName,
    pp=pp, pp.by.sequence=pp.by.sequence, pp.by.patient=pp.by.patient
  )
  predOutcomeDF <- machineLearning.Stacking.Predict_From_PreprocessedDataset(
    ml.data.predict=ML_Data_Pred$"MLDataFrame",
    modelList=modelList, model.stack=model.stack
  )
  predOutcomeDF <- data.frame(ML_Data_Pred$"MLMetadataDF", predOutcomeDF)
  return(list("MLDataFrame"=ML_Data_Pred$"MLDataFrame", "MLPredDF"=predOutcomeDF))
}
