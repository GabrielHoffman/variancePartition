#' Fit linear (mixed) model, report variance fractions
#'
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Report fraction of variance attributable to each variable
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.   See Details.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param control control settings for \code{lmer()}
#' @param hideErrorsInBackend default FALSE.  If TRUE, hide errors in \code{attr(.,"errors")} and \code{attr(.,"error.initial")}
#' @param showWarnings default TRUE. Indicate model failures
#' @param BPPARAM parameters for parallel evaluation
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#'
#' @return
#' list() of where each entry is a model fit produced by \code{lmer()} or \code{lm()}
#'
#' @details
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' If there are no random effects, so formula is ~ a + b + c, a 'standard' linear model is used:
#'
#' \code{fit <- lm( exprObj[j,] ~ a + b + c, data=data)}
#'
#' In both cases, \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run in parallel.  Parallel processing must be enabled before calling this function.  See below.
#'
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to \code{lm}/\code{lmer}.
#'
#' \code{REML=FALSE} uses maximum likelihood to estimate variance fractions.  This approach produced unbiased estimates, while \code{REML=TRUE} can show substantial bias.  See Vignette "3) Theory and practice of random effects and REML"
#'
#' @examples
#'
#' # load library
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
#' # Step 1: fit linear mixed model on gene expression
#' # If categorical variables are specified, a linear mixed model is used
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
#' # Note: fitExtractVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package = "Biobase")
#'
#' # ExpressionSet example
#' form <- ~ (1 | sex) + (1 | type) + score
#' info2 <- Biobase::pData(sample.ExpressionSet)
#' varPart2 <- fitExtractVarPartModel(sample.ExpressionSet, form, info2)
#'
# # Parallel processing using multiple cores with reduced memory usage
# param = SnowParam(4, progressbar=TRUE)
# varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2, BPPARAM = param)
#'
#' @export
#' @docType methods
#' @rdname fitExtractVarPartModel-method
#' @importFrom BiocParallel SerialParam bpiterate bplapply bpok
setGeneric("fitExtractVarPartModel",
  signature = "exprObj",
  function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {
    standardGeneric("fitExtractVarPartModel")
  }
)


# matrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,matrix-method
setMethod(
  "fitExtractVarPartModel", "matrix",
  function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {

    weights <- matrix(1, nrow(exprObj), ncol(exprObj))
    exprObj <- new("EList", list(E = exprObj, weights = weights))

    .fitExtractVarPartModel(exprObj, formula, data,
      REML = REML,
      useWeights = useWeights,
      control = control,
      hideErrorsInBackend = hideErrorsInBackend,
      showWarnings = showWarnings,
      BPPARAM = BPPARAM,
      ...
    )
  }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod(
  "fitExtractVarPartModel", "data.frame",
  function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {

    fitExtractVarPartModel(as.matrix(exprObj), formula, data,
      REML = REML,
      useWeights = useWeights,
      control = control,
      hideErrorsInBackend = hideErrorsInBackend,
      showWarnings = showWarnings,
      BPPARAM = BPPARAM,
      ...
    )
  }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod(
  "fitExtractVarPartModel", "EList",
  function(exprObj, formula, data, REML = FALSE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {

    .fitExtractVarPartModel(exprObj, formula, data,
      REML = REML,
      useWeights = useWeights,
      control = control,
      hideErrorsInBackend = hideErrorsInBackend,
      showWarnings = showWarnings,
      BPPARAM = BPPARAM,
      ...
    )
  }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod(
  "fitExtractVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {
    
    exprObj <- as.matrix(exprs(exprObj))

    fitExtractVarPartModel(exprObj, formula, data,
      REML = REML,
      useWeights = useWeights,
      control = control,
      hideErrorsInBackend = hideErrorsInBackend,
      showWarnings = showWarnings,
      BPPARAM = BPPARAM,
      ...
    )
  }
)

# sparseMatrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,sparseMatrix-method
setMethod(
  "fitExtractVarPartModel", "sparseMatrix",
  function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, hideErrorsInBackend = FALSE, showWarnings = TRUE, BPPARAM = SerialParam(), ...) {
    .fitExtractVarPartModel(exprObj, formula, data,
      REML = REML,
      useWeights = useWeights,
      control = control,
      hideErrorsInBackend = hideErrorsInBackend,
      showWarnings = showWarnings,
      BPPARAM = BPPARAM,
      ...
    )
  }
)




# internal driver function
.fitExtractVarPartModel <- function(exprObj,
                                    formula,
                                    data,
                                    REML = FALSE,
                                    useWeights = TRUE,
                                    control = vpcontrol,
                                    hideErrorsInBackend = FALSE,
                                    showWarnings = TRUE,
                                    BPPARAM = SerialParam(),
                                    ...) {
  # filter and check input data
  objFlt <- filterInputData(exprObj, formula, data, useWeights = useWeights)

  res <- run_lmm(objFlt$exprObj, objFlt$formula, objFlt$data,
    fxn = calcVarPart,
    REML = REML,
    control = control,
    dreamCheck = FALSE,
    BPPARAM = BPPARAM,
    ...
  )

  if (!is.null(res$error.initial) & !hideErrorsInBackend) {
    stop(paste0("Initial model failed:\n", res$error.initial))
  }

  # Get results from successful models
  res.unlist <- unlist(res$succeeded, recursive = FALSE)

  if (length(res.unlist) == 0 & !hideErrorsInBackend) {
    txt <- paste("All models failed.  The first model fails with:\n", res$errors[1])
     stop(txt)
  }

  # Aggregate successful results
  df <- as.data.frame(Reduce(rbind, res.unlist))
  rownames(df) <- names(res.unlist)

  modelType <- ifelse(.isMixedModelFormula(formula),
    "linear mixed model", "anova"
  )

  res2 <- new("varPartResults", df,
    type = modelType,
    method = "Variance explained (%)"
  )

  attr(res2, "errors") <- res$errors
  attr(res2, "error.initial") <- res$error.initial

  if (!is.null(attr(res2, "errors")) & !hideErrorsInBackend) {
    txt <- paste("Model failed for", length(attr(res2, "errors")), "responses.\n  See errors with attr(., 'errors')")
    warning(txt)
  }

  res2
}
