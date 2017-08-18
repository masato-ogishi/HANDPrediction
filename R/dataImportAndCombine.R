# Load & combine a sequence dataset and a metadata into a single dataframe
#' @title Load & combine a sequence dataset and a metadata into a single dataframe
#' @param fileName.Alignment.FASTA A file name of the HIV env sequence alignment in a FASTA format.
#' @param fileName.HIVNeuroMetadata.CSV A file name of the metadata in a CSV format.
#' @param fileName.SampleTissueCategoryDesignSheet.CSV A file name of the tissue category design sheet in a CSV format.
#' @param colName.SequenceID The column name representing sequene ID in the metadata.
#' @param colName.Sample.Tissue The column name representing sample source tissues in the metadata.
#' @param colName.Clinical.Status (Optional) The column name representing neurological outcomess in the metadata.
#' @importFrom stringr str_replace_all
#' @importFrom stringr str_detect
#' @importFrom stringr str_count
#' @importFrom stringr str_length
#' @importFrom seqinr read.alignment
#' @importFrom stringr str_replace
#' @importFrom stringr fixed
#' @importFrom dplyr filter
#' @importFrom dplyr select
#' @importFrom stringdist stringdistmatrix
#' @export
dataImportAndCombine <- function(
  fileName.Alignment.FASTA,
  fileName.HIVNeuroMetadata.CSV,
  fileName.SampleTissueCategoryDesignSheet.CSV,
  colName.SequenceID,
  colName.Sample.Tissue,
  colName.Clinical.Status=NULL
){
  # Load & filter sequences
  filteredSeqAlignment <- function(alignmentFileName, sequenceIDColName="SequenceID", threshold=0.1){
    SeqFilter <- function(seq){
      s <- stringr::str_replace_all(seq, "-", "X") # Alignment gaps will be regarded as ambiguities (=X).
      s.nogap <- stringr::str_replace_all(seq, "-", "")
      removeQ <- "Pass"
      if(stringr::str_detect(s, stringr::fixed("*"))) removeQ <- "Fail"  # Sequences containing a stop codon will be discarded.
      if(stringr::str_count(s, stringr::fixed("X")) > stringr::str_length(s)*threshold) removeQ <- "Fail"  # Sequences containing too many X will be discarded.
      list("Sequence.AA.Unaligned"=s.nogap, "Sequence.AA"=s, "RemoveQ"=removeQ)
    }
    s <- seqinr::read.alignment(alignmentFileName, format="fasta", forceToLower=F)
    d <- as.data.frame(sapply(s$seq, SeqFilter))
    d <- as.data.frame(t(d))
    rownames(d) <- 1:s$nb
    d[[sequenceIDColName]] <- stringr::str_replace(s$nam, stringr::fixed("\r"), "")
    d$"Sequence.AA.Unaligned" <- toupper(unlist(d$"Sequence.AA.Unaligned"))
    d$"Sequence.AA" <- toupper(unlist(d$"Sequence.AA"))
    d <- dplyr::filter(d, RemoveQ=="Pass")
    d <- dplyr::select(d, -RemoveQ)
    return(d)
  }
  df <- filteredSeqAlignment(fileName.Alignment.FASTA)

  # Input validation
  HXB2.C2V3C3 <- "NTSVITQACPKVSFEPIPIHYCAPAGFAILKCNNKTFNGTGPCTNVSTVQCTHGIRPVVSTQLLLNGSLAEEEVVIRSVNFTDNAKTIIVQLNTSVEINCTRPNNNTRKRIRIQRGPGRAFVTIGKIGNMRQAHCNISRAKWNNTLKQIASKLREQFGNNKTIIFKQSSGGDPEIVTHSFNCGGEFFY"
  editDistanceRatios <- as.vector(stringdist::stringdistmatrix(HXB2.C2V3C3, df$"Sequence.AA", method="lv"))/nchar(HXB2.C2V3C3)
  if(max(editDistanceRatios)>0.5){
    print("Note: Extremely large deviation from the HXB2 reference sequence was detected. Check your input sequences.")
  }
  df[["Sequence.AA.VarianceFromHXB2"]] <- editDistanceRatios

  # Combine metadata
  df.meta <- read.csv(fileName.HIVNeuroMetadata.CSV)
  colnames(df.meta)[grep(colName.SequenceID, colnames(df.meta), value=F)] <- "SequenceID"
  df <- merge(df, df.meta, by="SequenceID", all.x=T, sort=F)

  # Combine sample tissue annotation data
  colnames(df)[grep(colName.Sample.Tissue, colnames(df), value=F)] <- "Sample.Tissue"
  df.sampleTissueCategories <- read.csv(fileName.SampleTissueCategoryDesignSheet.CSV)
  df <- merge(df, df.sampleTissueCategories, all.x=T, sort=F, by="Sample.Tissue")
  df$"Sample.Tissue.Category" <- as.factor(df$"Sample.Tissue.Category")

  # Combine neurological outcome annotation data (optional)
  if(!is.null(colName.Clinical.Status)){
    colnames(df)[grep(colName.Clinical.Status, colnames(df), value=F)] <- "Clinical.Status"
    df$"Clinical.Status" <- as.factor(df$"Clinical.Status")
  }

  # Assign unique sequence IDs
  seq.unique <- unique(df$"Sequence.AA.Unaligned")
  seq.unique.names <- paste0("Unique_", formatC(1:length(seq.unique), width=ceiling(log10(length(seq.unique))), flag="0"))
  df <- merge(df, data.frame("UniqueID"=seq.unique.names, "Sequence.AA.Unaligned"=seq.unique), all.x=T, sort=F, by="Sequence.AA.Unaligned")

  # Reorder columns
  if(!is.null(colName.Clinical.Status)){
    cols <- setdiff(colnames(df),
                    c("SequenceID", "Clinical.Status", "Sample.Tissue", "Sample.Tissue.Category",
                      "UniqueID", "Sequence.AA", "Sequence.AA.Unaligned", "Sequence.AA.VarianceFromHXB2"))
    df <- df[c("SequenceID", "Clinical.Status", "Sample.Tissue", "Sample.Tissue.Category",
               cols,
               "UniqueID", "Sequence.AA.Unaligned", "Sequence.AA", "Sequence.AA.VarianceFromHXB2")]
  } else {
    cols <- setdiff(colnames(df),
                    c("SequenceID", "Sample.Tissue", "Sample.Tissue.Category",
                      "UniqueID", "Sequence.AA", "Sequence.AA.Unaligned", "Sequence.AA.VarianceFromHXB2"))
    df <- df[c("SequenceID", "Sample.Tissue", "Sample.Tissue.Category",
               cols,
               "UniqueID", "Sequence.AA.Unaligned", "Sequence.AA", "Sequence.AA.VarianceFromHXB2")]
  }

  # A combined metadataset
  return(df)
}
