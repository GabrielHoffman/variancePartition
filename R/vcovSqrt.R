# Gabriel Hoffman
# May 11, 2023
#
# Evaluate the sqrt of the variance-covarince
# matrix on result of dream() to



# TODO
# contrasts
# coefs: cant subset vcov matrix, must do shur complement
# in decorrelate, set eclairs argument to sample size = sqrt(n)








#' Sqrt of co-variance matrix for \code{dream()} fit
#'
#' Define generic \code{vcovSqrt()} for result of \code{lmFit()} and \code{dream()}
#
#' @param object \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} object returned by \code{voom()}
#' @param coef name of coefficient to be extracted
#' @param approx use fast approximation
#'
#' @return Computes factor of covariance matrix so that \code{vcov(object)} is the same as \code{crossprod(vcovSqrt(object))}
#'
#' @examples
#' # load simulated data:
#' # geneExpr: matrix of *normalized* gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' form <- ~Batch
#'
#' fit <- dream(geneExpr[1:2, ], form, info)
#' fit <- eBayes(fit)
#'
#' # Compute covariance directly
#' Sigma <- vcov(fit, geneExpr[1:2, ])
#'
#' # Compute factor of covariance
#' S <- crossprod(vcovSqrt(fit, geneExpr[1:2, ]))
#'
#' @importFrom stats coefficients
#' @docType methods
#' @rdname vcovSqrt-method
#' @export
setGeneric("vcovSqrt",
  signature = "object",
  function(object, vobj, coef, approx = TRUE) {
    standardGeneric("vcovSqrt")
  }
)




#' @rdname vcovSqrt-method
#' @aliases mvTest,MArrayLM-method
#' @export
setMethod("vcovSqrt", c("MArrayLM"), function(object, vobj, coef, approx = TRUE) {
  if (missing(vobj)) {
    stop("Must include response object as argument")
  }

  # Check that model fit and vobj have the same number of responses
  if (nrow(object) != nrow(vobj)) {
    stop("Model fit and data must have the same number of responses")
  }

  # Check that the responses have the same name
  if (!identical(rownames(object), rownames(vobj))) {
    stop("Model fit and data must have the same responses")
  }

  if (object$method != "ls") {
    stop("Only valid for models fit with lmFit()")
  }

  if (is(vobj, "EList")) {
    weights <- t(vobj$weights)
    colnames(weights) <- rownames(vobj)
    rownames(weights) <- colnames(vobj)
  } else {
    weights <- matrix(1, ncol(vobj), nrow(vobj))
    colnames(weights) <- rownames(vobj)
    rownames(weights) <- colnames(vobj)
  }

  # check that coef is valid
  if (!missing(coef)) {
    i <- match(coef, colnames(coefficients(object)))
    if (any(is.na(i))) {
      txt <- paste("Coefficients not valid:", paste(coef[is.na(i)], collapse = ", "))
      stop(txt)
    }
  } else {
    coef <- NULL
  }

  # subsetting MArrayLM objects keep all residuals
  # so subset manually here
  features <- rownames(coefficients(object))

  idx <- match(features, rownames(residuals(object)))

  resids <- t(residuals(object)[idx, , drop = FALSE])

  # use exact calculation for linear model
  # vcov = crossprod(P)
  P <- eval_vcov_sqrt(
    resids = resids,
    X = object$design,
    W = weights[rownames(resids), , drop = FALSE],
    rdf = object$df.residual[1],
    coef = coef,
    contrasts = object$contrasts,
    approx = approx
  )

  P
})




#' @rdname vcovSqrt-method
#' @aliases mvTest,MArrayLM2-method
#' @export
setMethod("vcovSqrt", c("MArrayLM2"), function(object, vobj, coef, approx = TRUE) {
  if (missing(vobj)) {
    stop("Must include response object as argument")
  }

  # Check that model fit and vobj have the same number of responses
  if (nrow(object) != nrow(vobj)) {
    stop("Model fit and data must have the same number of responses")
  }

  # Check that the responses have the same name
  if (!identical(rownames(object), rownames(vobj))) {
    stop("Model fit and data must have the same responses")
  }

  if (is(vobj, "EList")) {
    weights <- t(vobj$weights)
    colnames(weights) <- rownames(vobj)
    rownames(weights) <- colnames(vobj)
  } else {
    weights <- matrix(1, ncol(vobj), nrow(vobj))
    colnames(weights) <- rownames(vobj)
    rownames(weights) <- colnames(vobj)
  }

  # check that coef is valid
  if (!missing(coef)) {
    i <- match(coef, colnames(object$cov.coefficients.list[[1]]))
    if (any(is.na(i))) {
      txt <- paste("Coefficients not valid:", paste(coef[is.na(i)], collapse = ", "))
      stop(txt)
    }
  } else {
    coef <- NULL
  }

  # subsetting MArrayLM objects keep all residuals
  # so subset manually here
  features <- rownames(coefficients(object))

  idx <- match(features, rownames(residuals(object)))

  resids <- t(residuals(object)[idx, , drop = FALSE])

  # use approximate calculation for linear mixed model
  eval_vcov_sqrt_approx(
    resids = resids,
    W = weights[rownames(resids), , drop = FALSE],
    ccl = object$cov.coefficients.list,
    X = object$design,
    coef = coef,
    contrasts = object$contrasts
  )
})


eval_vcov_sqrt <- function(resids, X, W, rdf, coef, contrasts, approx) {
  # which coefficients to include
  if (!is.null(contrasts)) {
    if (is.null(coef)) coef <- colnames(contrasts)
  } else {
    if (is.null(coef)) coef <- colnames(X)
  }

  # scale weights to have mean 1 for each column
  W <- sweep(W, 2, colMeans(W), "/")

  # pre-compute square root of W
  sqrtW <- sqrt(W)

  # weighted residuals
  R <- resids * sqrtW

  # number of features
  m <- ncol(resids)

  res <- lapply(seq(m), function(i) {
    X_i <- X * sqrtW[, i]
    a <- R[, i, drop = FALSE] / sqrt(rdf)
    if (approx) {
      b <- matrExp(solve(crossprod(X_i)), 0.5)
    } else {
      b <- t(solve(crossprod(X_i), t(X_i)))
    }
    kronecker(a, b, make.dimnames = TRUE)
  })
  P <- do.call(cbind, res)

  if (!is.null(contrasts)) {
    # use contrasts

    # extract single contrast
    L <- contrasts[, coef, drop = FALSE]

    # expand L to multiple responses
    D <- bdiag(lapply(seq(m), function(i) L))
    D <- as.matrix(D)

    # assign names
    colnames(D) <- c(outer(coef, colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # apply linear contrasts
    P <- P %*% D
  } else {
    # use coef
    # subsect to selected coefs
    i <- match(coef, colnames(X))

    # names of coefficients to retain
    keep <- c(outer(colnames(X)[i], colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # subset covariance matrix
    P <- P[, keep, drop = FALSE]
  }
  P
}


eval_vcov_sqrt_approx <- function(resids, W, ccl, X, coef, contrasts) {
  if (identical(colnames(contrasts), rownames(contrasts))) {
    contrasts <- NULL
  }

  # which coefficients to include
  if (!is.null(contrasts)) {
    if (is.null(coef)) coef <- colnames(contrasts)
  } else {
    if (is.null(coef)) coef <- colnames(X)
  }

  # scale weights to have mean 1
  W <- sweep(W, 2, colMeans(W), "/")

  # store dimensions of data
  k <- ncol(ccl[[1]])
  m <- ncol(resids)

  # residual correlation
  scale_res <- scale(resids * sqrt(W)) / sqrt(nrow(resids) - 1)

  res <- lapply(seq(m), function(i) {
    a <- scale_res[, i, drop = FALSE]
    b <- matrExp(ccl[[i]], 0.5)
    kronecker(a, b, make.dimnames = TRUE)
  })
  P <- do.call(cbind, res)

  if (!is.null(contrasts)) {
    coef <- coef[coef %in% colnames(contrasts)]

    # names of coefficients to retain
    keep <- c(outer(coef, colnames(resids), function(a, b) paste(b, a, sep = ":")))

    P <- P[, keep, drop = FALSE]
  } else {
    # use coef
    # subsect to selected coefs
    i <- match(coef, colnames(X))

    # names of coefficients to retain
    keep <- c(outer(colnames(X)[i], colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # subset covariance matrix
    P <- P[, keep, drop = FALSE]
  }
  P
}
