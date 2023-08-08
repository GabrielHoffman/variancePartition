# Gabriel Hoffman
# Nov 1, 2021

#' Apply pre-specified sample weights
#'
#' Apply pre-specified sample weights by scaling existing precision weights
#'
#' @param vobj \code{EList} from \code{voom} or \code{voomWithDreamWeights}.
#' @param weights sample level weights
#'
#' @details Apply pre-specified sample-level weights to the existing precision weights estimated from the data.  While the \code{limma::voomWithQualityWeights} function of Lui et al. (2015) estimates the sample-level weights from \code{voom} fit, here the weights are fixed beforehand.
#'
#' @references{
#'   \insertRef{liu2015weight}{variancePartition}
#' }
#'
#' @export
#' @seealso \code{limma::voomWithQualityWeights}
applyQualityWeights <- function(vobj, weights) {
  # apply weights like in voomWithQualityWeights
  vobj$weights <- t(weights * t(vobj$weights))

  vobj$targets$sample.weights <- weights

  vobj
}
