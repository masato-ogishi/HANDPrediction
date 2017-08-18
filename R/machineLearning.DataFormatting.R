# Convert the positions of the input alignment to the positions in HXB2 env reference sequence.
#' @title Convert the positions of the input alignment to the positions in HXB2 env reference sequence.
#' @param startPosition A position from which the C2 region starts; a default is 197 in the HXB2 sequence.
#' @param endPosition A position at which the C3 region ends; a default is 384 in the HXB2 sequence.
#' @param gapPositions Gap positions; a default is c(125), which is consistent with the alignment used in the original study. Users should carefully examine their alignments before setting this parameter.
#' @importFrom purrr is_empty
#' @export
#' @rdname HXB2_C2V3C3_AlignmentPositionConverter
HXB2_C2V3C3_AlignmentPositionConverter <- function(startPosition=197, endPosition=384, gapPositions=c(125)){
  vec <- paste0("Pos", formatC(startPosition:endPosition, width=3, flag="0"))
  if(is.na(gapPositions)) return(vec)
  if(is.null(gapPositions)) return(vec)
  if(purrr::is_empty(gapPositions)) return(vec)
  for(i in 1:length(gapPositions)){
    vec <- append(vec, values=paste0("Pos",(startPosition-1+gapPositions[[i]]),"G"),
                  after=(gapPositions[[i]]+i-1)) #"G" for an alignment gap
  }
  return(vec)
}

# Convert an amino acid sequence to a numerical matrix.
#' @title Convert an amino acid sequence to a numerical matrix.
#' @param sequence An amino acid sequence of HIV C2V3C3 region.
#' @param sequences A set of amino acid sequences.
#' @param gapPositions Alignment gap positions.
#' @param analysisType An indicator for the type of  analysis. Currently, "boolian" and "AAIndex" are implemented.
#' @param aaIndexType An indicator for the type of the set of AAIndex to be used. Currently, "all" and "physicochemical" are implemented. If set to be "all", a total of 544 AAIndices included in the BioSeqClass package are used. If set to be "physicochemical", pre-selected 78 AAIndices are used. See the original paper for further details.
#' @param normalized Logical. Whether AAIndex values are rescaled to [0,1] range.
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @export
#' @rdname machineLearning.Matrix
aaseqToBoolian <- function(sequence){
  seq <- str_replace_all(sequence, "\\r", "")
  seq.aaList <- unlist(str_split(seq, ""))
  seq.aa.bool <- numeric(21*nchar(seq))
  for(i in 1:nchar(seq)){seq.aa.bool[aaLetters[seq.aaList[[i]]]+(i-1)*21] <- 1}
  return(seq.aa.bool)
}

#' @importFrom stringr str_replace_all
#' @importFrom stringr str_replace
#' @export
#' @rdname machineLearning.Matrix
aaseqToBoolian.batch <- function(sequences, gapPositions=c(125)){
  seqs <- str_replace_all(sequences, "\\r", "")
  len <- max(nchar(seqs))
  colHeaders <- list()
  for(i in 1:len){colHeaders[[i]] <- paste0(HXB2_C2V3C3_AlignmentPositionConverter(gapPositions=gapPositions)[[i]], "_AA_", names(aaLetters))}
  colHeaders <- unlist(colHeaders)
  d <- as.data.frame(t(as.data.frame(sapply(seqs, aaseqToBoolian))))
  rownames(d) <- NULL
  colnames(d) <- colHeaders
  return(d)
}

#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @export
#' @rdname machineLearning.Matrix
aaseqToAAIndex <- function(sequence, aaIndexType="all", normalized=T){
  seq <- str_replace_all(sequence, "\\r", "")
  seq.aaList <- unlist(str_split(seq, ""))
  if(normalized==T){
    val <- switch(aaIndexType,
                  "all"=as.vector(unlist(aaIndexData.all.norm[seq.aaList])),
                  "physicochemical"=as.vector(unlist(aaIndexData.pc.norm[seq.aaList])))
  }
  if(normalized==F){
    val <- switch(aaIndexType,
                  "all"=as.vector(unlist(aaIndexData.all[seq.aaList])),
                  "physicochemical"=as.vector(unlist(aaIndexData.pc[seq.aaList])))
  }
  return(val)
}

#' @importFrom stringr str_replace_all
#' @importFrom stringr str_split
#' @export
#' @rdname machineLearning.Matrix
aaseqToAAIndex.batch <- function(sequences, gapPositions=c(125), aaIndexType="all", normalized=T){
  seqs <- str_replace_all(sequences, "\\r", "")
  len <- max(nchar(seqs))
  inds <- switch(aaIndexType,
                 "all"=aaIndexData.all$"aaIndexName",
                 "physicochemical"=aaIndexData.pc$"aaIndexName")
  colHeaders <- list()
  for(i in 1:len){colHeaders[[i]] <- paste0(HXB2_C2V3C3_AlignmentPositionConverter(gapPositions=gapPositions)[[i]], "_AAIndex_", inds)}
  colHeaders <- unlist(colHeaders)
  d <- as.data.frame(t(as.data.frame(lapply(seqs, function(s){aaseqToAAIndex(sequence=s, aaIndexType=aaIndexType, normalized=normalized)}))))
  rownames(d) <- NULL
  colnames(d) <- colHeaders
  return(d)
}

#' @importFrom stringr str_replace_all
#' @export
#' @rdname machineLearning.Matrix
machineLearning.Matrix <- function(sequences, gapPositions=c(125),
                                   analysisType="AAIndex", aaIndexType="physicochemical", normalized=T){
  # Analysis flow chart
  type <- "aaindex.physicochemical"
  if(analysisType=="boolian") type <- "boolian"
  if(aaIndexType=="all") type <- "aaindex.all"

  # A numerical matrix
  seqs <- stringr::str_replace_all(sequences, "\\r", "")
  mat <- switch(type,
              "boolian"=aaseqToBoolian.batch(seqs, gapPositions=gapPositions),
              "aaindex.all"=aaseqToAAIndex.batch(seqs, gapPositions=gapPositions, aaIndexType="all", normalized=normalized),
              "aaindex.physicochemical"=aaseqToAAIndex.batch(seqs, gapPositions=gapPositions, aaIndexType="physicochemical", normalized=normalized)
  )
  rownames(mat) <- paste0("Seq_", formatC(1:length(seqs), width=nchar(as.character(length(seqs))), flag="0"))
  return(mat)
}

# Convert the sequence-derived numerical matrix into a machine-learning-ready dataframe.
#' @title Convert the sequence-derived numerical matrix into a machine-learning-ready dataframe.
#' @param ml.matrix A numerical matrix returned by the machineLearning.Matrix function.
#' @param ml.metadata A metadata.
#' @param colNames.from Names of the columns in the metadata to be incorporated in the subsequent analysis.
#' @param colNames.to (Optional) Column names of the included metadata. Default values are the same as colNames.from.
#' @param outcomeLabelName A column name representing the neurological outcomes.
#' @param sourceLabelName A column name representing the sample sources.
#' @param pp An indicator for preprocessing. If set "internal", variables are first screened using the zero-variance and near-zero-variance methods implemented in the caret package, and then further screened using the correlation method implemented in the caret package. If set "external", the matrix is pre-processed using externally provided functions. Omit preprocessing by setting "none".
#' @param pp.by.sequence (Optional) Users can provide their own pre-processing function (per sequence).
#' @param pp.by.patient (Optional) Users can provide their own pre-processing function (per patient).
#' @importFrom caret preProcess
#' @importFrom caret dummyVars
#' @importFrom dplyr select_at
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise_at
#' @importFrom dplyr vars
#' @importFrom dplyr funs
#' @importFrom dplyr contains
#' @importFrom dplyr ungroup
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @export
#' @rdname machineLearning.DataFormat
machineLearning.DataFormat <- function(ml.matrix, ml.metadata,
                                       colNames.from, colNames.to=colNames.from,
                                       outcomeLabelName=NULL, sourceLabelName,
                                       pp="internal", pp.by.sequence=NULL, pp.by.patient=NULL){
  d <- as.data.frame(ml.matrix)

  # Filtering & pre-processing (sequence-level)
  preProc.by.sequence <- NULL
  if(pp=="internal"){
    pp.zv <- caret::preProcess(d, c("zv", "nzv"))
    d <- predict(pp.zv, d)
    pp.cor <- caret::preProcess(d, c("corr"))
    d <- predict(pp.cor, d)
    preProc.by.sequence <- colnames(d)
  }
  if(pp=="external"){
    d <- dplyr::select(d, dplyr::one_of(pp.by.sequence))
  }
  if(pp=="none"){
    d <- d
  }

  # Combine metadata and values
  for(i in 1:length(colNames.from)){
    d[[colNames.to[[i]]]] <- ml.metadata[[colNames.from[[i]]]]
  }

  # Summarize values by groups
  d_grouped <- d %>%
    dplyr::group_by_(.dots=colNames.to) %>%
    dplyr::summarise_at(dplyr::vars(dplyr::contains("_AA")), dplyr::funs(mean, median, min, max, sd), na.rm=T) %>%
    dplyr::ungroup()
  d_grouped_metadata <- d_grouped %>% dplyr::select(dplyr::one_of(colNames.to))
  d_grouped <- d_grouped %>% dplyr::select(-dplyr::one_of(colNames.to))

  # Filtering & pre-processing (patient-level)
  preProc.by.patient <- NULL
  if(pp=="internal"){
    pp.zv <- caret::preProcess(d_grouped, c("zv", "nzv"))
    d_grouped <- predict(pp.zv, d_grouped)
    pp.cor <- caret::preProcess(d_grouped, c("corr"))
    d_grouped <- predict(pp.cor, d_grouped)
    pp.cs <- caret::preProcess(d_grouped, c("center", "scale"))
    d_grouped <- predict(pp.cs, d_grouped)
    preProc.by.patient <- list("Columns"=colnames(d_grouped), "Model"=pp.cs)
  }
  if(pp=="external"){
    d_grouped <- d_grouped %>% dplyr::select(dplyr::one_of(pp.by.patient$"Columns"))
    d_grouped <- predict(pp.by.patient$"Model", d_grouped)
  }
  if(pp=="none"){
    d_grouped <- d_grouped
  }

  # Create dummy variables for source tissues
  f <- as.formula(paste0("~", sourceLabelName))
  d_grouped_source <- predict(caret::dummyVars(f, data=d_grouped_metadata), d_grouped_metadata) %>% as.data.frame()
  d_grouped_wide <- dplyr::bind_cols(d_grouped_source, d_grouped)

  # Append outcome (Optional)
  if(!is.null(outcomeLabelName)){
    d_out <- try(dplyr::select_at(d_grouped_metadata, outcomeLabelName), silent=T)
    if(class(d_out)!="try-error") d_grouped_wide <- dplyr::bind_cols(d_out, d_grouped_wide)
  }

  # Remove NAs
  rows.noNA <- complete.cases(d_grouped_wide)
  d_grouped_wide <- d_grouped_wide[rows.noNA,]
  d_grouped_metadata <- d_grouped_metadata[rows.noNA,] # NAs in the metadata should be tolerated

  # Remove unused factor levels
  d_grouped_wide <- droplevels(d_grouped_wide)
  d_grouped_metadata <- droplevels(d_grouped_metadata)

  # Results
  return(list("MLDataFrame"=d_grouped_wide, "MLMetadataDF"=d_grouped_metadata,
              "MLDataPreProcessing_BySequence"=preProc.by.sequence,
              "MLDataPreProcessing_ByPatient"=preProc.by.patient))
}
