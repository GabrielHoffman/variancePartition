# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai


#' Class varParFrac
#'
#' Class \code{varParFrac}
#'
#' @name varParFrac-class
#' @rdname varParFrac-class
#' @exportClass varParFrac
setClass("varParFrac")

#' Class VarParFitList
#'
#' Class \code{VarParFitList}
#'
#' @name VarParFitList-class
#' @rdname VarParFitList-class
#' @exportClass VarParFitList
setClass("VarParFitList", representation(method = "character"), contains = "list")


#' Class varPartResults
#'
#' Class \code{varPartResults}
#'
#' @name varPartResults-class
#' @rdname varPartResults-class
#' @exportClass varPartResults
setClass("varPartResults", representation(type = "character", method = "character"), contains = "data.frame")

# # @export
# setMethod("print", "varPartResults",
#   function(x, ...) {
#   	print( x )
#   }
# )
# setMethod("print", "varPartResults",
#   function(x, ...) {
#   	print( x )
#   }
# )

#' Convert to matrix
#'
#' Convert varPartResults to matrix
#'
#' @param x varPartResults
#' @param ... other arguments.
#'
#' @return
#' matrix
#' @examples
#' # load library
#' # library(variancePartition)
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
#' # Fit model
#' varPart <- fitExtractVarPartModel(geneExpr[1:5, ], form, info)
#'
#' # convert to matrix
#' as.matrix(varPart)
#'
#' @export
# @docType methods
#' @aliases as.matrix
setMethod(
  "as.matrix", "varPartResults",
  function(x, ...) {
    df <- as.data.frame(x@.Data)
    colnames(df) <- names(x)
    rownames(df) <- x@row.names
    # colnames(df) = colnames(x)
    # rownames(df) = rownames(x)

    return(as.matrix(df))
  }
)


#' Convert to data.frame
#'
#' Convert varPartResults to data.frame
#'
#' @param x varPartResults
#' @param row.names pass thru to generic
#' @param optional pass thru to generic
#' @param ... other arguments.
#'
#' @return
#' data.frame
#' @examples
#' # load library
#' # library(variancePartition)
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
#' # Fit model
#' varPart <- fitExtractVarPartModel(geneExpr[1:5, ], form, info)
#'
#' # convert to matrix
#' as.data.frame(varPart)
#'
#' @export
#' @rawNamespace S3method("as.data.frame", 'varPartResults')
# @rdname as.data.frame
# @method as.data.frame varPartResults
# S3/S4 combo for varPartResults
as.data.frame.varPartResults <- function(x, row.names = NULL, optional = FALSE, ...) {
  df <- as.data.frame(x@.Data)
  colnames(df) <- names(x)
  rownames(df) <- x@row.names
  df
}

# @rdname as.data.frame-methods
# @aliases as.data.frame,varPartResults,varPartResults-method
# setMethod("as.data.frame", "varPartResults", as.data.frame.varPartResults)




#' Default colors for ggplot
#'
#' Return an array of n colors the same as the default used by ggplot2
#'
#' @param n number of colors
#'
#' @return
#' array of colors of length n
#'
#' @examples
#' ggColorHue(4)
#' @export
ggColorHue <- function(n) {
  hues <- seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#' Residuals from model fit
#'
#' Extract residuals for each gene from model fit with fitVarPartModel()
#'
#' @param object object produced by fitVarPartModel()
#' @param ... other arguments.
#'
#' @return
#' Residuals extracted from model fits stored in object
#'
#' @details
#' If model is fit with missing data, residuals returns NA for entries that were
#' missing in the original data
#'
#' @examples
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
#' # Fit model
#' modelFit <- fitVarPartModel(geneExpr, form, info)
#'
#' # Extract residuals of model fit
#' res <- residuals(modelFit)
#'
#' @export
# @docType methods
#' @aliases residuals,VarParFitList-method
setMethod(
  "residuals", "VarParFitList",
  function(object, ...) {
    # extract residuals
    res <- lapply(object, function(fit) {
      residuals(fit)
    })

    resMatrix <- matrix(NA, nrow = length(res), ncol = length(res[[1]]))
    colnames(resMatrix) <- names(res[[1]])

    # fill non-missing entries with residuals
    for (j in 1:nrow(resMatrix)) {
      resMatrix[j, ] <- res[[j]]
    }

    # get total number of samples
    n_samples <- sapply(res, length)

    if (max(n_samples) - min(n_samples) > 0) {
      stop("Samples were dropped from model fit.  Either expression data or metadata contained NA values")
    }

    # identify which samples were omitted
    # excludeList <- sapply( object, function(fit)
    # 	getOmitted( fit )
    # )

    # # get total number of samples
    # n_samples = sapply(res, length)# + sapply(excludeList, length)

    # if( max(n_samples) - min(n_samples) > 0 ){
    # 	stop("Samples were dropped from model fit.  Either expression data or metadata contained NA values")
    # }

    # # create matrix of total size, including omitted samples
    # resMatrix <- matrix(NA, nrow=length(res), ncol=n_samples[1])

    # # fill non-missing entries with residuals
    # for( j in 1:nrow(resMatrix)){
    # 	excl = excludeList[[j]]

    # 	if( is.null(excl) ){
    # 		resMatrix[j,] = res[[j]]
    # 	}else{
    # 		resMatrix[j,-excl] = res[[j]]
    # 	}
    # }

    # # gene names along rows
    # rownames(resMatrix) <- names(object)

    # # get columns names, filling NA values from excludeList
    # sampleNames = rep(NA, ncol(resMatrix))
    # excl = excludeList[[1]]

    # if( is.null(excl) ){
    # 	sampleNames = names(fitted.values(object[[1]]))
    # }else{
    # 	sampleNames[-excl] = names(fitted.values(object[[1]]))
    # 	sampleNames[excl] = names(excl)
    # }
    # colnames(resMatrix) = sampleNames

    return(as.matrix(resMatrix))
  }
)

setGeneric("getOmitted",
  signature = "fit",
  function(fit, ...) {
    standardGeneric("getOmitted")
  }
)

setMethod(
  "getOmitted", "lm",
  function(fit, ...) {
    fit$na.action
  }
)

setMethod(
  "getOmitted", "lmerMod",
  function(fit, ...) {
    attr(fit@frame, "na.action")
  }
)




#' Collinearity score
#'
#' Collinearity score for a regression model indicating if variables are too highly correlated to give meaningful results
#'
#' @param fit regression model fit from lm() or lmer()
#'
#' @return
#' Returns the collinearity score between 0 and 1, where a score > 0.999 means the degree of collinearity is too high.  This function reports the correlation matrix between coefficient estimates for fixed effects.  The collinearity score is the maximum absolute correlation value of this matrix. Note that the values are the correlation between the parameter estimates, and not between the variables themselves.
#' @export
#' @examples
#'
# load data
#' # load library
#' # library(variancePartition)
#'
#' # load simulated data:
#' data(varPartData)
#' #
#' form <- ~ Age + (1 | Individual) + (1 | Tissue)
#'
#' res <- fitVarPartModel(geneExpr[1:10, ], form, info)
#'
#' # evaluate the collinearity score on the first model fit
#' # this reports the correlation matrix between coefficients estimates
#' # for fixed effects
#' # the collinearity score is the maximum absolute correlation value
#' # If the collinearity score > .999 then the variance partition
#' # estimates may be problematic
#' # In that case, a least one variable should be omitted
#' colinearityScore(res[[1]])
#'
colinearityScore <- function(fit) {
  # get correlation matrix
  V <- vcov(fit)
  if (any(is.na(vcov(fit))) || nrow(V) == 0) {
    C <- NA
  } else {
    C <- cov2cor(as.matrix(V))
  }

  score <- 0
  if (any(is.na(vcov(fit)))) {
    score <- 1
  } else if (length(C) > 1 && nrow(C) > 1) {
    # return largest off-diagonal absolute correlation
    score <- max(abs(C[lower.tri(C)]))
  }

  attr(score, "vcor") <- C
  return(score)
}




#' Check if model contains a random effect
#'
#' Check if model contains a random effect
#'
#' @param formula model formula
#'
#' @importFrom stats as.formula
# @importFrom lme4 lmerControl
#' @importFrom reformulas findbars
#' @keywords internal
.isMixedModelFormula <- function(formula) {
  !is.null(findbars(as.formula(formula)))
}






# return array indicating if variable in formula has missing data
hasMissingData <- function(form, info) {
  # get variable names
  variableNames <- all.vars(form)

  # return indicator of any NA's
  sapply(variableNames, function(v) {
    any(is.na(info[[v]]))
  })
}





#' Adapted from lme4:::reOnly
#'
#' Adapted from lme4:::reOnly
#'
#' @param f formula
#' @param response (FALSE) is there a response in the formula
#'
#' @importFrom reformulas findbars
#' @importFrom stats reformulate
reOnly <- function(f, response = FALSE) {
  reformulate(paste0(
    "(", vapply(findbars(f), deparse1, ""),
    ")"
  ), response = if (response && length(f) == 3L) {
    f[[2]]
  })
}
