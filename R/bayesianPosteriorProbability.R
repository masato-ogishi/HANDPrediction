# Visualize the distribution of Bayesian posterior probabilities in ggplot
#' @title Visualize the distribution of Bayesian posterior probabilities in ggplot
#' @param predictedOutcomes The prediction result. Either "HAND", "NonHAND", or both.
#' @param priorProbHAND Prior probability of HAND. It is the user's responsibility to provide a rational estimate. The default is set to be 0.33, which is derived from The Multicenter AIDS Cohort Study [PMID: 26718568].
#' @param confMat The caret confusion matrix of the stacked classifier.
#' @param sensitivity (Optional) The sensitivity of the stacked classifier.
#' @param specificity (Optional) The specificity of the stacked classifier.
#' @importFrom dplyr %>%
#' @importFrom dplyr mutate
#' @importFrom tidyr gather
#' @import ggplot2
#' @export
#' @rdname bayesianPosteriorProbability
bayesianPosteriorProbability <- function(predictedOutcomes=c("HAND","NonHAND"), priorProbHAND=0.33,
                                        confMat=NULL, sensitivity=NULL, specificity=NULL){
  postProbs <- function(p.prior, sn, sp, out){
    odds.prior <- p.prior/(1-p.prior)
    likelihood <- switch(out,
                         "HAND"=sn/(1-sp),
                         "NonHAND"=(1-sn)/sp,
                         1)
    odds.posterior <- odds.prior*likelihood
    return(odds.posterior/(1+odds.posterior))
  }
  if(is.null(confMat)){
    res <- sapply(predictedOutcomes, function(o){postProbs(priorProbHAND, sensitivity, specificity, o)})
  }else{
    res <- sapply(predictedOutcomes, function(o){
      postProbs(priorProbHAND, confMat$byClass["Sensitivity"], confMat$byClass["Specificity"], o)
    })
  }
  names(res) <- predictedOutcomes
  return(res)
}

#' @export
#' @rdname bayesianPosteriorProbability
bayesianPosteriorProbabilityPlot <- function(predictedOutcomes=c("HAND","NonHAND"),
                                            confMat=NULL, sensitivity=NULL, specificity=NULL){
  priorProbHANDSet <- seq(0.05, 0.95, by=0.05)
  d <- lapply(priorProbHANDSet,
              function(p){bayesianPosteriorProbability(predictedOutcomes, priorProbHAND=p, confMat)}) %>%
    as.data.frame() %>% t() %>% as.data.frame() %>%
    dplyr::mutate(PriorHANDProbability=priorProbHANDSet) %>%
    tidyr::gather(PredictionResult, PosteriorHANDProbability, -PriorHANDProbability)
  ggplot(d, aes(x=PriorHANDProbability, y=PosteriorHANDProbability, color=PredictionResult)) +
    geom_smooth(method='loess') + xlim(0,1) + xlab("Prior probability of HAND") + ylab("Posterior probability of HAND") +
    scale_color_brewer(palette="Set1", name="Predicted HAND status") +
    theme_Publication()
}
