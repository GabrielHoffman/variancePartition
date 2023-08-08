#' Class VarParCIList
#'
#' Class \code{VarParCIList}
#'
#' @name VarParCIList-class
#' @rdname VarParCIList-class
#' @exportClass VarParCIList
setClass("VarParCIList", representation(method = "character"), contains = "list")


#' Linear mixed model confidence intervals
#'
#' Fit linear mixed model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Then perform parametric bootstrap sampling to get a 95\% confidence intervals for each variable for each gene.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  Strongly discourage against changing this option, but here for compatibility.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList} from \code{voom()} or \code{weightsMatrix} is specified
#' @param control control settings for \code{lmer()}
#' @param nsim number of bootstrap datasets
#' @param ... Additional arguments for \code{lmer()} or l\code{m()}
#'
#' @return
#' \code{list()} of where each entry is the result for a gene.  Each entry is a matrix of the 95\% confidence interval of the variance fraction for each variable
#'
#' @details
#' A linear mixed model is fit for each gene, and \code{bootMer()} is used to generate parametric bootstrap confidence intervals.  \code{use.u=TRUE} is used so that the \eqn{\hat{u}} values from the random effects are used as estimated and are not re-sampled.  This gives confidence intervals as if additional data were generated from these same current samples.  Conversely, \code{use.u=FALSE} assumes that  this dataset is a sample from a larger population.   Thus it simulates \eqn{\hat{u}} based on the estimated variance parameter.  This approach gives confidence intervals as if additional data were collected from the larger population from which this dataset is sampled.  Overall, \code{use.u=TRUE} gives smaller confidence intervals that are appropriate in this case.
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
#' # Compute bootstrap confidence intervals for each variable for each gene
#' resCI <- varPartConfInf(geneExpr[1:5, ], form, info, nsim = 100)
#'
#' @importFrom stats as.formula
#' @importFrom lme4 bootMer
#' @export
varPartConfInf <- function(exprObj, formula, data, REML = FALSE, useWeights = TRUE, control = vpcontrol, nsim = 1000, ...) {
  if (!is.numeric(nsim) || nsim <= 0) {
    stop("Must specify nsim as positive number")
  }

  # need derivs in boostrapping
  control$calc.derivs <- TRUE

  formula <- as.formula(formula)

  # define bootstrap function
  bootStrapFxn <- local(function(fit) {
    if (is(fit, "lmerMod")) {
      # use bootMer from lme4 to sample to bootstraps and refit the model
      # use.u=TRUE says that the \hat{u} values from the random effects are used
      bootRes <- bootMer(fit, calcVarPart, use.u = TRUE, nsim = nsim, parallel = "no", ncpus = 1)
      apply(bootRes$t, 2, function(x) quantile(x, c(0.025, 0.975), na.rm = TRUE))
    } else if (is(fit, "lm")) {
      stop("Bootstrap of fixed effect ANOVA model not currently supported")
    }
  })

  # fit the model and run the bootstrap function for each gene
  res <- fitVarPartModel(exprObj = exprObj, formula = formula, data = data, REML = REML, useWeights = useWeights, fxn = bootStrapFxn, control = control, ...)

  # res@method = "bootstrap"

  return(res)
}




# form = ~ Age + Individual + Tissue

#  fit = lm(geneExpr[1,]~ Age + Individual + Tissue, info, weights=1:nrow(info))



# bs2 <- function(formula, data, indices) {
# 	d <- data[indices,]
# 	fit2 <- lm(formula, data=d, weights=fit$weights[indices])
# 	return( calcVarPart(fit2) )
# }

# bootRes <- boot(data=data.frame(info, gene=geneExpr[1,]), statistic=bs2, R=1000, formula= paste("gene ~", as.character(form))



# apply( bootRes$t, 2, function(x) quantile(x, c(0.025, 0.975)))

# calcVarPart( fit )
