# HANDPrediction: toward sequence-based diagnosis of HIV-associated neurocognitive disorder (HAND)

Introduction
---------------------------
This is the R package for the prediction of the probability of HAND for given patients from sets of HIV envelope C2V3C3 sequences.
- Users should prepare following files; (1) env C2V3C3 alignment in FASTA format; (2) metadata; (3) sample source category design sheet.
- Users are recommended to prepare their alignment using [HIVAlign tool](https://www.hiv.lanl.gov/content/sequence/VIRALIGN/viralign.html) with the HMM-align checkbox checked.
- The metadata should contain columns named as 'SequenceID', 'Clinical.Status' and 'Sample.Tissue'.
- The sample shource design sheet should contain columns named as 'Sample.Tissue' and 'Sample.Tissue.Category'. Currently, only supported Sample.Tissue.Category values are 'CNS', 'Blood', 'Lymph', and 'Others'.

Installation
---------------------------
1. Install the latest version of [R](https://cran.r-project.org/bin/windows/base/) (>=3.4.1) and [RStudio](https://www.rstudio.com/products/rstudio/download2/).
2. Upgrade [Bioconductor](https://www.bioconductor.org/install/) to the latest version (>=3.5).
3. Install the package <i>HANDPrediction</i>.  
``` r
if(!require(devtools)) install.packages("devtools")
devtools::install_github("masato-ogishi/HANDPrediction")

# Alternatively, please download this repository as a Zip file, unzip it to the directory you want, and run the following command.
devtools::install_local("path/to/the/unzipped/folder")
```

Reference
---------------------------
- Masato Ogishi and Hiroshi Yotsuyanagi. Prediction of HIV-associated neurocognitive disorder (HAND) from three genetic features of envelope gp120 glycoprotein. Retrovirology. 2018 Jan 27;15(1):12. doi: 10.1186/s12977-018-0401-x.
