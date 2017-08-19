# Screen important features.
#' @title Screen important features.
#' @param ml.data A machine-learning-ready dataframe returned by the machineLearning.DataFormat function.
#' @param outcomeLabelName A column name for neurological outcomes.
#' @param modelList A list of caret models.
#' @param topN An integer indicating how many the most important features are retained.
#' @param pThreshold A threshold of significance. Features with significantly distinct value distributions between neurological outcomes are considered important.
#' @param removedVars (Optional) Feature names that should be removed manually.
#' @param outputFileName.DOCX (Optional) A file name of a summary table in a Word format.
#' @param xLabel (Optional) An x-axis label for the violin plot.
#' @param yLabel (Optional) A y-axis label for the violin plot.
#' @param colorSet The color set to be used.
#' @importFrom dplyr %>%
#' @importFrom DescTools Sort
#' @importFrom caret varImp
#' @importFrom dplyr slice
#' @importFrom stringr str_split
#' @importFrom purrr is_empty
#' @importFrom optimbase transpose
#' @importFrom dplyr bind_rows
#' @importFrom ReporteRs borderProperties
#' @importFrom ReporteRs FlexTable
#' @importFrom ReporteRs spanFlexTableRows
#' @importFrom ReporteRs setFlexTableBorders
#' @importFrom ReporteRs docx
#' @importFrom ReporteRs addFlexTable
#' @importFrom ReporteRs writeDoc
#' @importFrom dplyr group_by_
#' @importFrom dplyr summarise
#' @importFrom dplyr filter
#' @importFrom data.table transpose
#' @importFrom dplyr select
#' @importFrom dplyr one_of
#' @importFrom dplyr ungroup
#' @importFrom tidyr gather
#' @importFrom stringr str_replace
#' @import ggplot2
#' @export
machineLearning.FeatureImportanceAnalysis <- function(
  ml.data, outcomeLabelName,
  modelList, topN=20, pThreshold=0.05, removedVars=c(),
  outputFileName.DOCX=NULL,
  xLabel=NULL, yLabel=NULL, colorSet
){
  # Screen important features in each model
  algorithmLabelList <- names(modelList)
  var.list <- list()
  pos.list <- list()
  imp.list <- list()
  for(i in 1:length(algorithmLabelList)){
    df.imp <- suppressWarnings(try(caret::varImp(modelList[[algorithmLabelList[[i]]]]), silent=T))
    if(!class(df.imp)=="try-error"){
      df.imp <- DescTools::Sort(df.imp[["importance"]], decreasing=T)
      df.imp[["Feature"]] <- rownames(df.imp)
      df.imp <- dplyr::slice(df.imp, grep("AAIndex_", rownames(df.imp), value=F))
      var.list[[i]] <- df.imp[["Feature"]][1:topN]
      pos.list[[i]] <- sort(unique(stringr::str_split(df.imp[["Feature"]][1:topN], "_", simplify=T)[,1]))
      imp.list[[i]] <- df.imp[[1]][1:topN]
    }
  }
  var.all <- sort(unique(unlist(var.list)))
  pos.all <- sort(unique(unlist(pos.list)))
  d.list.pos <- list()
  for(i in 1:length(pos.all)){
    d.list.alg <- list()
    for(j in 1:length(algorithmLabelList)){
      v <- grep(pos.all[[i]], var.list[[j]], value=T)
      d <- data.frame()
      if(!purrr::is_empty(v)){
        d <- stringr::str_split(v, "_", simplify=T)[,c(3,4)]
        if(is.matrix(d)) d <- as.data.frame(d) else d <- as.data.frame(optimbase::transpose(as.matrix(d)))
        colnames(d) <- c("AAIndex", "Stat")
        d <- data.frame("Position"=pos.all[[i]], "Algorithm"=algorithmLabelList[[j]], d)
        imp <- imp.list[[j]][grep(pos.all[[i]], var.list[[j]], value=F)]
        d <- data.frame(d, "Importance"=imp)
      }
      d.list.alg[[j]] <- d
    }
    d.list.pos[[i]] <- suppressWarnings(dplyr::bind_rows(d.list.alg))
  }
  d <- suppressWarnings(dplyr::bind_rows(d.list.pos))

  # Export a formatted table to Word
  if(!is.null(outputFileName.DOCX)){
    no_border=borderProperties(width=0)
    big_border=borderProperties(width=2)
    std_border=borderProperties(width=1)
    table <- FlexTable(d,
                       header.cell.props=cellProperties(background.color="#003366"),
                       header.text.props=textBold(color="white"),
                       add.rownames=T)
    table <- spanFlexTableRows(table, j="Position", runs=as.character(d$"Position"))
    table <- spanFlexTableRows(table, j="Algorithm", runs=as.character(paste0(d$"Position",".",d$"Algorithm")))
    table <- setFlexTableBorders(table, footer=F,
                                 inner.vertical=no_border, inner.horizontal=std_border,
                                 outer.vertical=no_border, outer.horizontal=big_border)
    doc <- docx( )
    doc <- addFlexTable(doc, table)
    writeDoc(doc, file=outputFileName.DOCX)
  }

  # Collect features regarded as important by two or more algorithms
  d.imp <- d %>%
    dplyr::group_by_("Position","AAIndex","Stat") %>%
    dplyr::summarise(Count=n()) %>%
    dplyr::filter(Count>=2)
  var.imp <- sapply(as.list(data.table::transpose(dplyr::select(d.imp, -Count))),
                    function(v){
                      v[[2]] <- paste0("AAIndex_", v[[2]])
                      paste0(v, collapse="_")
                    }) %>% as.character()
  d.pvalues <- ml.data
  colnames(d.pvalues)[grep(outcomeLabelName, colnames(d.pvalues))] <- "Outcome"
  d.pvalues <- d.pvalues %>%
    dplyr::select(one_of(c("Outcome",var.imp))) %>%
    dplyr::group_by_("Outcome") %>%
    dplyr::filter(n()>1)  %>%
    dplyr::ungroup()
  pvalue.list <- sapply(var.imp,
                        function(colName){
                          t.test(as.formula(paste0(colName, "~Outcome")), data=d.pvalues, var.equal=F)$p.value
                        }) ## Welch's t test
  pvalue.fdradjusted.list <- p.adjust(pvalue.list, method="BH") ## Adjusted by FDR
  d.pvalues <- data.frame("Feature"=var.imp, "P.Value"=pvalue.list, "Adj.P.Value"=pvalue.fdradjusted.list,
                          "P.Value.Label"=paste0("P = ", format.pval(pvalue.list, digits=1, eps=0.001)),
                          "Adj.P.Value.Label"=paste0("P = ", format.pval(pvalue.fdradjusted.list, digits=1, eps=0.001))) %>%
    dplyr::filter(Adj.P.Value < pThreshold) %>%
    dplyr::filter(!Feature %in% removedVars)
  var.imp <- unique(as.character(d.pvalues$"Feature"))

  # Violin plot for important features
  d.plot <- ml.data
  colnames(d.plot)[grep(outcomeLabelName, colnames(d.plot))] <- "Outcome"
  d.plot[["Outcome"]] <- factor(d.plot[["Outcome"]], levels=rev(levels(d.plot[["Outcome"]]))) # rev for coord_flip
  d.plot <- d.plot %>%
    dplyr::select(one_of(c("Outcome", var.imp))) %>%
    tidyr::gather(Feature, Value, -Outcome)
  d.pvalues.plot <- data.frame(d.pvalues, "Outcome"=1.5, "Value"=-2.2)
  var.label.FUN <- function(strings){ str_replace(strings, "AAIndex_", "") }
  violine.plot <- suppressWarnings(
    ggplot(data=d.plot, aes_string(x="Outcome", y="Value", color="Outcome")) +
      geom_violin(scale="width", width=0.9, trim=F, fill="#999999", linetype="blank", alpha=I(1/3)) +
      stat_summary(geom="pointrange",fun.y=mean, fun.ymin=function(x){mean(x)-sd(x)}, fun.ymax=function(x){mean(x)+sd(x)}, size=1, alpha=.5) +
      geom_text(data=d.pvalues.plot, aes_string(label="Adj.P.Value.Label"), color="black", size=4) +
      facet_wrap(~Feature, labeller=labeller(Feature=var.label.FUN)) +
      ylim(-3,3) +
      xlab(xLabel) + ylab(yLabel) +
      scale_color_manual(values=colorSet) + guides(color="none") +
      coord_flip() +
      theme_Publication()
  )

  # Outputs
  return(list("FeatureImportanceDF"=d, "FeaturePValueDF"=d.pvalues, "MostImportantFeatures"=var.imp, "ViolinPlot"=violine.plot))
}

# Residue frequencies and AAIndex values at specified important positions.
#' @title Residue frequencies and AAIndex values at specified important positions.
#' @param featureSet An important feature set.
#' @param sequences Sequences to be analyzed.
#' @param seqWeights Weights per sequence. Patient-to-patient differences in sequencing depth should be normalized. Can be provided either as fractions (<1) or sequence counts (>1).
#' @param outcomes A character vector containing neurological outcomes.
#' @param xLabel An x-axis label.
#' @param yLabel1 A y1-axis label [Residue frequency].
#' @param yLabel2 A y2-axis label [AAIndex value].
#' @param legendLabel A legend label.
#' @param colorSet The color set to be used.
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom dplyr select
#' @importFrom dplyr select_at
#' @importFrom dplyr bind_rows
#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom stringr str_split
#' @importFrom magrittr set_colnames
#' @importFrom tidyr gather
#' @importFrom optimbase transpose
#' @importFrom scales rescale
#' @import ggplot2
#' @export
machineLearning.ResidueAAIndexDualPlot <- function(
  featureSet, sequences, seqWeights=NULL, outcomes,
  xLabel=NULL, yLabel1=NULL, yLabel2=NULL, legendLabel=NULL, colorSet
){
  # Normalize sequencing depths
  if(is.null(seqWeights)){
    sequences_norm <- sequences
    outcomes_norm <- droplevels(outcomes)
  } else {
    if(min(seqWeights)>=1){
      w <- 1/seqWeights
    } else {
      w <- seqWeights
    } # w must be fractions.
    sequences_norm <- mapply(function(x,y){rep(x, times=round(1000*y))}, sequences, w) %>% unlist()
    outcomes_norm <- mapply(function(x,y){rep(x, times=round(1000*y))}, outcomes, w) %>% unlist() %>% droplevels()
  }

  # AAIndex values
  d_aaindex <- data.frame(
    "Position"=as.vector(unlist(as.data.frame(stringr::str_split(featureSet,"_"))[1,])),
    "AAIndex"=as.vector(unlist(as.data.frame(stringr::str_split(featureSet,"_"))[3,]))
  ) %>% dplyr::mutate(Label1=paste0(Position,",",AAIndex))
  d_aaindex <- merge(
    aaIndexData.all.norm %>%
      dplyr::select(-aaIndexName, -aaIndexDescription) %>%
      t() %>% as.data.frame() %>%
      magrittr::set_colnames(aaIndexData.all.norm$"aaIndexName") %>%
      dplyr::select_at(as.character(d_aaindex$"AAIndex")) %>%
      dplyr::mutate(Residue=names(aaLetters)) %>%
      tidyr::gather(AAIndex, Value, -Residue) %>%
      dplyr::mutate(Residue=factor(Residue, levels=names(aaLetters))),
    d_aaindex,
    all=T, by="AAIndex", sort=F
  ) %>% dplyr::mutate(Label2=paste0(Position,",",Residue))

  # AA residue frequencies
  d_freq <- as.data.frame(stringr::str_split(sequences_norm, "")) %>%
    as.matrix() %>% optimbase::transpose() %>% as.data.frame()
  colnames(d_freq) <- HXB2_C2V3C3_AlignmentPositionConverter()
  d_freq_sub <- list()
  p <- unique(as.character(d_aaindex$"Position"))
  for(i in 1:length(p)){
    d_freq_sub[[i]] <- as.data.frame(table(outcomes_norm, d_freq[[p[i]]]))
    colnames(d_freq_sub[[i]]) <- c("Outcome", "Residue", "Count")
    d_freq_sub[[i]][["Position"]] <- p[i]
  }
  d_freq <- suppressWarnings(
    dplyr::bind_rows(d_freq_sub) %>%
      dplyr::mutate(Residue=factor(Residue, levels=names(aaLetters))) %>%
      na.omit() %>%
      dplyr::mutate(Label2=paste0(Position,",",Residue)) %>%
      dplyr::select(-Residue, -Position)
  )
  d <- merge(d_aaindex, d_freq, all=T, by="Label2", sort=F) %>% dplyr::mutate(Label3=paste0(Label1,",",Outcome))
  d <- merge(d, d %>% group_by(Label3) %>% dplyr::summarise(Total.Count=sum(Count)), all.x=T, by="Label3", sort=F)
  d[["Freq"]] <- d$"Count"/d$"Total.Count"
  d[["Residue.Num"]] <- as.numeric(d[["Residue"]])
  d <- na.omit(d)

  # Combined ggplot
  scale_to_value1 <- function(values) scales::rescale(values, to=range(na.omit(d)$Freq))
  scale_to_value2 <- function(values) scales::rescale(values, to=range(d$Value))
  ggplot(NULL) +
    geom_bar(data=na.omit(d), aes(x=Residue.Num, y=Freq, fill=Outcome), stat="identity", position="dodge") +
    geom_line(data=d, aes(x=Residue.Num, y=scale_to_value1(Value)), stat="identity", position="identity") +
    geom_point(data=d, aes(x=Residue.Num, y=scale_to_value1(Value))) +
    facet_wrap(~Label1, scales="free_y") +
    scale_x_continuous(xLabel, breaks=1:20, labels=names(aaLetters)[1:20]) +
    scale_y_continuous(yLabel1, sec.axis=sec_axis(~ scale_to_value2(.), name=yLabel2)) +
    scale_fill_manual(values=colorSet, name=legendLabel) +
    theme_Publication() +
    theme(legend.position="top", legend.direction="horizontal")
}

