library(tidyverse)
library(stringr)
library(seqinr)
library(Biostrings)
library(BioSeqClass)
library(scales)
aaLetters <- 1:21
names(aaLetters) <- c(names(AMINO_ACID_CODE)[1:20],"X") # "X" for an ambiguity.
ambiguityValue <- 0 # A surrogate AAIndex value for "X"
data(aa.index)
aaIndexNames <- names(aa.index)
aaIndexDescriptions <- c()
for(i in 1:length(aaIndexNames)){aaIndexDescriptions[[i]] <- aa.index[[i]][["D"]]}
aaIndexData <- data.frame("aaIndexName"=aaIndexNames, "aaIndexDescription"=aaIndexDescriptions)
aaIndexValues <- list()
for(i in 1:length(aaIndexNames)){aaIndexValues[[aaIndexNames[[i]]]] <- aa.index[[i]][["I"]]}
aaIndexValues <- as.data.frame(t(as.data.frame(aaIndexValues)))
aaIndexValues[["aaIndexName"]] <- rownames(aaIndexValues)
aaIndexData.all <- merge(x=aaIndexData, y=aaIndexValues, by="aaIndexName", all=T, sort=F) %>%
  dplyr::mutate(X=ambiguityValue) %>%
  (function(d){d[complete.cases(d),]})
aaIndex.pc <-
  c(grep("Hydro", aaIndexDescriptions, ignore.case=T, value=T),
    grep("Charge", aaIndexDescriptions, ignore.case=T, value=T),
    grep("Polar", aaIndexDescriptions, ignore.case=T, value=T),
    grep("Distribution", aaIndexDescriptions, ignore.case=T, value=T),
    grep("Flexi", aaIndexDescriptions, ignore.case=T, value=T))
aaIndexData.pc <- aaIndexData.all %>% dplyr::filter(aaIndexDescription %in% aaIndex.pc)
aaIndexData.all.norm <- aaIndexData.all %>%
  dplyr::select(-aaIndexName, -aaIndexDescription, -X) %>%
  t() %>% as.data.frame() %>%
  (function(dat){dplyr::bind_rows(lapply(dat, rescale))}) %>%
  t() %>% as.data.frame() %>%
  dplyr::mutate(X=ambiguityValue) %>%
  magrittr::set_colnames(names(aaLetters))
aaIndexData.all.norm <- dplyr::bind_cols(dplyr::select(aaIndexData.all, aaIndexName, aaIndexDescription), aaIndexData.all.norm)
aaIndexData.pc.norm <- aaIndexData.all.norm %>% dplyr::filter(aaIndexDescription %in% aaIndex.pc)

## devtools::use_data(aaLetters, aaIndexData.all, aaIndexData.pc, aaIndexData.all.norm, aaIndexData.pc.norm, internal = TRUE, overwrite = TRUE)
