# Compare multiple machine learning algorithms.
#' @title Compare multiple machine learning algorithms.
#' @param trainingData A training dataset.
#' @param testingData A testing dataset.
#' @param ml.data A machine-learning-ready dataframe returned by the machineLearning.DataFormat function.
#' @param ml.metadata A metadata.
#' @param holdoutTrainingRatio A ratio of data to be used for model training.
#' @param outcomeLabelName A column name for neurological outcomes.
#' @param patientIDLabelName A column name for patient IDs.
#' @param algorithmList Machine learning algorithms to be compared.
#' @param algorithmLabelList User-defined labels for the machine learning algorithms to be compared.
#' @param trainControlOptions A trainControl object for the train function in caret.
#' @param stackingAlgorithm An algorithm used for model stacking.
#' @param tuneL An integer denoting the amount of granularity in the grid search of tuning parameters in the train function.
#' @param outputFileName.TXT (Optional) Print and save the comparison results in a TXT format.
#' @param seed A random seed.
#' @importFrom dplyr select_
#' @importFrom dplyr bind_cols
#' @importFrom caret preProcess
#' @importFrom caret train
#' @importFrom caret confusionMatrix
#' @export
#' @rdname machineLearning.Comparison
machineLearning.Comparison.default <- function(trainingData, testingData, outcomeLabelName,
                                               algorithmList, algorithmLabelList, trainControlOptions,
                                               stackingAlgorithm="xgbTree", tuneL=3){
  # Remove numeric predictors with no/little variation from the training data
  trainingData.pp <- caret::preProcess(dplyr::select_(trainingData, paste0("-",outcomeLabelName)), method=c("zv", "nzv"))
  trainingData.pp <- predict(trainingData.pp, newdata=trainingData)

  # Construct various classifiers
  modelList <- list(1:length(algorithmList))
  f <- as.formula(paste0(outcomeLabelName, " ~ ."))
  for(i in 1:length(algorithmList)){
    modelList[[i]] <- invisible(caret::train(f, data=trainingData.pp, method=algorithmList[[i]], trControl=trainControlOptions, tuneLength=tuneL))
  }
  names(modelList) <- algorithmLabelList

  # Compare the performances of the best-tuned classifiers
  predList <- list(1:length(algorithmList))
  for(i in 1:length(algorithmList)){
    predList[[i]] <- predict(modelList[[i]], newdata=testingData)
  }
  names(predList) <- algorithmLabelList
  cmList <- list(1:length(algorithmList))
  for(i in 1:length(algorithmList)){
    cmList[[i]] <- caret::confusionMatrix(predList[[i]], testingData[[outcomeLabelName]])
  }
  names(cmList) <- algorithmLabelList

  # Stack all the classifiers
  predProbList <- list(1:length(algorithmList))
  for(i in 1:length(algorithmList)){
    predProbList[[i]] <- predict(modelList[[i]], newdata=trainingData, type="prob")
  }
  names(predProbList) <- algorithmLabelList
  predProbDF <- as.data.frame(predProbList)
  predProbDF[[outcomeLabelName]] <- trainingData[[outcomeLabelName]]
  model.stack <- caret::train(f, data=predProbDF, method=stackingAlgorithm, trControl=trainControlOptions, tuneLength=tuneL)

  # Evaluate the performance of a stacked classifier
  predProbList.test <- list(1:length(algorithmList))
  for(i in 1:length(algorithmList)){
    predProbList.test[[i]] <- predict(modelList[[i]], newdata=testingData, type="prob")
  }
  names(predProbList.test) <- algorithmLabelList
  predProbDF.test <- as.data.frame(predProbList.test)
  predProbDF.test[[outcomeLabelName]] <- testingData[[outcomeLabelName]]
  pred.test.stack <- predict(model.stack, predProbDF.test)
  cm.stack <- caret::confusionMatrix(pred.test.stack, testingData[[outcomeLabelName]])
  modelList[["Stack"]] <- model.stack
  cmList[["Stack"]] <- cm.stack

  # Sensitivity, specificity etc...
  snList <- sapply(cmList, function(cm){cm[["byClass"]]['Sensitivity']})
  spList <- sapply(cmList, function(cm){cm[["byClass"]]['Specificity']})
  accList <- sapply(cmList, function(cm){cm[["overall"]]['Accuracy']})
  accLoList <- sapply(cmList, function(cm){cm[["overall"]]['AccuracyLower']})
  accUpList <- sapply(cmList, function(cm){cm[["overall"]]['AccuracyUpper']})
  statDF <- data.frame("Algorithm"=unlist(c(algorithmLabelList,"Stack")),
                       "Sensitivity"=snList, "Specificity"=spList, "Accuracy"=accList, "AccuracyLower"=accLoList, "AccuracyUpper"=accUpList)
  statDF[["Algorithm"]] <- factor(statDF[["Algorithm"]], levels=c(algorithmLabelList,"Stack"))
  rownames(statDF) <- NULL

  # Outputs
  return(list("modelList"=modelList, "cmList"=cmList, "statDF"=statDF))
}

#' @importFrom dplyr %>%
#' @importFrom caret createDataPartition
#' @importFrom dplyr mutate
#' @importFrom dplyr mutate_
#' @importFrom dplyr filter
#' @export
#' @rdname machineLearning.Comparison
machineLearning.Comparison <- function(ml.data, ml.metadata, holdoutTrainingRatio=0.8,
                                       outcomeLabelName, patientIDLabelName,
                                       algorithmList, algorithmLabelList, trainControlOptions,
                                       stackingAlgorithm="xgbTree", tuneL=3,
                                       outputFileName.TXT=NULL, seed=123456789){
  set.seed(seed)

  # Prepare datasets for model training and model validation
  holdoutIDs.DF <- dplyr::select_(ml.metadata, outcomeLabelName, patientIDLabelName) %>% unique.data.frame()
  holdoutIDs <- caret::createDataPartition(holdoutIDs.DF[[outcomeLabelName]], p=holdoutTrainingRatio, list=F)
  patientIDs_training <- holdoutIDs.DF[holdoutIDs,][[patientIDLabelName]]
  patientIDs_testing <- holdoutIDs.DF[-holdoutIDs,][[patientIDLabelName]]
  rowIDs.DF_training <- ml.metadata %>%
    dplyr::mutate(RowID=row.names(.)) %>%
    dplyr::mutate_(PatID=patientIDLabelName) %>%
    dplyr::filter(PatID %in% patientIDs_training)
  rowIDs.DF_testing <- ml.metadata %>%
    dplyr::mutate(RowID=row.names(.)) %>%
    dplyr::mutate_(PatID=patientIDLabelName) %>%
    dplyr::filter(PatID %in% patientIDs_testing)
  ML_training <- ml.data[as.numeric(rowIDs.DF_training$"RowID"),]
  ML_testing <- ml.data[as.numeric(rowIDs.DF_testing$"RowID"),]
  ML_training_Meta <- data.frame("Data"="Train", ml.metadata[as.numeric(rowIDs.DF_training$"RowID"),], ML_training[-grep(outcomeLabelName, colnames(ML_training), value=F)])
  ML_testing_Meta <- data.frame("Data"="Test", ml.metadata[as.numeric(rowIDs.DF_testing$"RowID"),], ML_testing[-grep(outcomeLabelName, colnames(ML_testing), value=F)])

  # Perform machine learning
  ML_Comparison_List <- machineLearning.Comparison.default(
    trainingData=ML_training, testingData=ML_testing, outcomeLabelName=outcomeLabelName,
    algorithmList=algorithmList, algorithmLabelList=algorithmLabelList,
    trainControlOptions=trainControlOptions,
    stackingAlgorithm=stackingAlgorithm, tuneL=tuneL
    )

  # Compare the performances of various classifiers
  if(!is.null(outputFileName.TXT)){
    sink(outputFileName.TXT)
    print(ML_Comparison_List)
    sink()
  }

  # Outputs
  return(list("MLDataFrame_Training"=ML_training_Meta, "MLDataFrame_Testing"=ML_testing_Meta,
              "MLModelList"=ML_Comparison_List$"modelList", "MLConfusionMatrixList"=ML_Comparison_List$"cmList",
              "MLStatDF"=ML_Comparison_List$"statDF"))

}


#' @param ... Parameters for the machineLearning.Comparison function.
#' @param seedList A set of random seeds.
#' @importFrom dplyr %>%
#' @importFrom dplyr bind_rows
#' @importFrom tidyr gather
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom Rmisc CI
#' @importFrom dplyr filter
#' @importFrom ggsci scale_color_npg
#' @export
#' @rdname machineLearning.Comparison
machineLearning.Comparison.Bootstrap <- function(..., seedList=1:10){
  # Machine learning with different random seeds
  remove_NA <- function(l){l[!is.na(l)]}
  statDFList <- list()
  for(i in 1:length(seedList)){
    statDFList[[i]] <- try(machineLearning.Comparison(..., outputFileName.TXT=NULL, seed=seedList[[i]]), silent=T)
    if(class(statDFList[[i]])!="try-error"){
      statDFList[[i]] <- statDFList[[i]]$"MLStatDF"
      statDFList[[i]][["RandomSeed"]] <- seedList[[i]]
    }else{
      NA
    }
  } %>% remove_NA()
  statDF <- dplyr::bind_rows(statDFList)
  statDF.summary <- statDF %>% tidyr::gather(Stat, Value, -Algorithm, -RandomSeed)
  statDF.summary[["Algorithm"]] <- factor(statDF.summary[["Algorithm"]], levels=unique(statDF[["Algorithm"]]))
  statDF.summary <- statDF.summary %>%
    dplyr::group_by(Algorithm, Stat) %>%
    dplyr::summarise(mean=mean(Value), ci.95.lo=Rmisc::CI(Value)["lower"], ci.95.up=Rmisc::CI(Value)["upper"])

  # The mean accuracy of stacked classifiers
  ave.acc <- statDF.summary %>%
    dplyr::filter(Algorithm=="Stack", Stat=="Accuracy")
  ave.acc <- ave.acc$"mean" %>% as.matrix() %>% as.vector()

  # Accuracy comparison plot
  acc.plot <- ggplot(statDF, aes(x=Algorithm, y=Accuracy, color=Algorithm)) +
    geom_jitter(shape=19, width=0.2) +
    stat_summary(fun.y="median", geom="crossbar", mapping=aes(ymin=..y.., ymax=..y..), width=0.5) +
    ylim(0,1) +
    ggsci::scale_color_npg(guide=F) +
    theme_Publication()

  # Output
  return(list("MLStatDF"=statDF, "MLStatDF.Summary"=statDF.summary, "StackMeanAccuracy"=ave.acc, "AccuracyPlot"=acc.plot))
}
