#' Extract variance statistics
#'
#' Extract variance statistics from list of models fit with \code{lm()} or \code{lmer()}
#'
#' @param modelList list of \code{lmer()} model fits
#' @param ... other arguments
#'
#' @return
#' \code{data.frame} of fraction of variance explained by each variable, after correcting for all others.
# @details
#' @examples
#' # library(variancePartition)
#'
#' library(BiocParallel)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1 | Individual) + (1 | Tissue)
#'
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel(geneExpr, form, info)
#'
#' # violin plot of contribution of each variable to total variance
#' plotVarPart(sortCols(varPart))
#'
#' # Advanced:
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results <- fitVarPartModel(geneExpr, form, info)
#'
#' # Step 2: extract variance fractions
#' varPart <- extractVarPart(results)
#'
#' @export
extractVarPart <- function(modelList, ...) {
  # get results from first model to enumerate all variables present
  singleResult <- calcVarPart(modelList[[1]], ...)

  # for each model fit, get R^2 values
  entry <- 1
  varPart <- lapply(modelList, function(entry) {
    calcVarPart(entry, ...)
  })

  varPartMat <- data.frame(matrix(unlist(varPart), nrow = length(varPart), byrow = TRUE))
  colnames(varPartMat) <- names(varPart[[1]])
  rownames(varPartMat) <- names(modelList)

  modelType <- ifelse(is(modelList[[1]], "lm"), "anova", "linear mixed model")

  new("varPartResults", varPartMat, type = modelType, method = "Variance explained (%)")
}
