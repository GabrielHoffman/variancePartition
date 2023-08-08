# Gabriel Hoffman
# October 1, 2022
#
# Evaluate vcov() on result of dream() to
# compute covariance between estimated coefficients

# TODO:
# Write vcov for two lm() or lmer() fits


#' Co-variance matrix for \code{dream()} fit
#'
#' Define generic \code{vcov()} for result of \code{lmFit()} and \code{dream()}
#
#' @param object \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} object returned by \code{voom()}
#' @param coef name of coefficient to be extracted
#'
#' @return variance-covariance matrix
#' @importFrom stats coefficients
#' @export
setMethod("vcov", c("MArrayLM"), function(object, vobj, coef) {
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
  eval_vcov(
    resids = resids,
    X = object$design,
    W = weights[rownames(resids), , drop = FALSE],
    rdf = object$df.residual[1],
    coef = coef,
    contrasts = object$contrasts
  )
})



#' Co-variance matrix for \code{dream()} fit
#'
#' Define generic \code{vcov()} for result of \code{lmFit()} and \code{dream()}
#
#' @param object \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} object returned by \code{voom()}
#' @param coef name of coefficient to be extracted
#'
#' @return variance-covariance matrix
#' @importFrom stats coefficients
#' @export
setMethod("vcov", c("MArrayLM2"), function(object, vobj, coef) {
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
  eval_vcov_approx(
    resids = resids,
    W = weights[rownames(resids), , drop = FALSE],
    ccl = object$cov.coefficients.list,
    X = object$design,
    coef = coef,
    contrasts = object$contrasts
  )
})



# Evaluate variance-covariance matrix in multivariate regression
#
# This method is exact for linear regression, even when each response has its own weight vector
#
# @param resids matrix of residuals from regression
# @param X design matrix
# @param W matrix of precision weights
# @param rdf residual degrees of freedom
# @param coef name of coefficient to be extracted
#
#' @importFrom Matrix bdiag
eval_vcov <- function(resids, X, W, rdf, coef, contrasts) {
  # With no weights:
  # kronecker(crossprod(res), solve(crossprod(X)), make.dimnames=TRUE) / rdf

  # which coefficients to include
  if (!is.null(contrasts)) {
    if (is.null(coef)) coef <- colnames(contrasts)
    coef <- unique(coef)
  } else {
    if (is.null(coef)) coef <- colnames(X)
  }

  # scale weights to have mean 1 for each column
  W <- sweep(W, 2, colMeans(W), "/")

  # pre-compute square root of W
  sqrtW <- sqrt(W)

  # all pairs of responses
  Sigma <- crossprod(resids * sqrtW)

  # store dimensions of data
  k <- ncol(X)
  m <- ncol(resids)

  # matrix to store results
  Sigma_vcov <- matrix(0, m * k, m * k)

  # Equivalent to kronecker product when weights are shared
  # outer loop
  for (i in seq(1, m)) {
    # define positions in output matrix
    idx1 <- seq((i - 1) * k + 1, i * k)

    # scale X by weights
    X_i <- X * sqrtW[, i]

    # evaluate (X^TX)^-1 X^T for i
    A <- solve(crossprod(X_i), t(X_i))

    # inner loop
    for (j in seq(i, m)) {
      idx2 <- seq((j - 1) * k + 1, j * k)
      X_j <- X * sqrtW[, j]
      B <- solve(crossprod(X_j), t(X_j))

      # standard method using observed covariates
      value <- (Sigma[i, j] / rdf) * tcrossprod(A, B)

      Sigma_vcov[idx1, idx2] <- value
      Sigma_vcov[idx2, idx1] <- value
    }
  }

  # assign names
  colnames(Sigma_vcov) <- c(outer(colnames(X), colnames(resids), function(a, b) paste(b, a, sep = ":")))
  rownames(Sigma_vcov) <- colnames(Sigma_vcov)

  if (!is.null(contrasts)) {
    # use contrasts

    # extract single contrast
    L <- contrasts[, coef, drop = FALSE]

    # expand L to multiple responses
    D <- bdiag(lapply(seq(m), function(i) L))

    # assign names
    colnames(D) <- c(outer(coef, colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # apply linear contrasts
    Sigma_vcov <- crossprod(D, Sigma_vcov) %*% D
    Sigma_vcov <- as.matrix(Sigma_vcov)
  } else {
    # use coef
    # subsect to selected coefs
    i <- match(coef, colnames(X))

    # names of coefficients to retain
    keep <- c(outer(colnames(X)[i], colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # subset covariance matrix
    Sigma_vcov <- Sigma_vcov[keep, keep]
  }
  Sigma_vcov
}


# Evaluate variance-covariance matrix in multivariate regression
#
# This method is approximate since part of the calculations assume equal weights.  This is useful for the linear mixed model where the exact calculation is very challanging
#
# @param resids matrix of residuals from regression
# @param W matrix of precision weights
# @param ccl list of vcov matrices for each response
# @param coef name of coefficient to be extracted
#
eval_vcov_approx <- function(resids, W, ccl, X, coef, contrasts) {
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
  Sigma <- crossprod(scale_res)

  # matrix to store results
  Sigma_vcov <- matrix(0, m * k, m * k)

  # Equivalent to kronecker product when weights are shared
  # outer loop
  for (i in seq(1, m)) {
    # define positions in output matrix
    idx1 <- seq((i - 1) * k + 1, i * k)

    # Use eigen-decomp and transforming eigen-values
    # since using contrasts can make this matrix singular
    sqrt_inv_i <- matrExp(ccl[[i]], -0.5)

    # inner loop
    for (j in seq(i, m)) {
      idx2 <- seq((j - 1) * k + 1, j * k)

      sqrt_inv_j <- matrExp(ccl[[j]], -0.5)

      value <- Sigma[i, j] * ccl[[i]] %*% crossprod(sqrt_inv_i, sqrt_inv_j) %*% ccl[[j]]

      Sigma_vcov[idx1, idx2] <- value
      Sigma_vcov[idx2, idx1] <- value
    }
  }

  # assign names
  colnames(Sigma_vcov) <- c(outer(colnames(ccl[[1]]), colnames(resids), function(a, b) paste(b, a, sep = ":")))
  rownames(Sigma_vcov) <- colnames(Sigma_vcov)

  # select coefficients
  if (!is.null(contrasts)) {
    coef <- coef[coef %in% colnames(contrasts)]

    # names of coefficients to retain
    keep <- c(outer(coef, colnames(resids), function(a, b) paste(b, a, sep = ":")))

    Sigma_vcov <- Sigma_vcov[keep, keep, drop = FALSE]
  } else {
    # use coef
    # subsect to selected coefs
    i <- match(coef, colnames(X))

    # names of coefficients to retain
    keep <- c(outer(colnames(X)[i], colnames(resids), function(a, b) paste(b, a, sep = ":")))

    # subset covariance matrix
    Sigma_vcov <- Sigma_vcov[keep, keep, drop = FALSE]
  }

  Sigma_vcov
}



# if( !is.null(contrasts)){
# 	# use contrasts
# 	# subsect to selected coefs
# 	coef = coef[coef %in% colnames(contrasts)]

# 	# names of coefficients to retain
# 	keep = c(outer(coef, colnames(resids), function(a,b) paste(b,a, sep=':')))

# 	# subset covariance matrix
# 	Sigma_vcov = Sigma_vcov[keep,keep]

# }else{
# 	# use coef from design matrix
# 	# names of coefficients to retain
# 	keep = c(outer(coef, colnames(resids), function(a,b) paste(b,a, sep=':')))

# 	# subset covariance matrix
# 	Sigma_vcov = Sigma_vcov[keep,keep]
# }



# Like standard sign function,
# except sign(x) giving 0 is reset to give 1
sign0 <- function(x) {
  # use standard sign function
  res <- sign(x)

  # get entries that equal 0
  # and set them to 1
  i <- which(res == 0)
  if (length(i) > 0) {
    res[i] <- 1
  }

  res
}

# Raise eigen-values of a matrix to exponent alpha
matrExp <- function(S, alpha, symmetric = TRUE, tol = sqrt(.Machine$double.eps)) {
  # pass R check
  vectors <- values <- NULL

  # eigen decomposition
  dcmp <- eigen(S, symmetric = symmetric)

  # Modify sign of vectors, so diagonal is always positive
  # This removes an issue of sensitivity to small numerical changes
  values <- sign0(diag(dcmp$vectors))
  dcmp$vectors <- sweep(dcmp$vectors, 2, values, "*")

  # identify positive eigen-values
  idx <- dcmp$values > tol

  # a matrix square root is U %*% diag(lambda^alpha) %*% U^T, alpha = 0.5
  # make sure to return to original axes
  res <- with(dcmp, vectors[, idx, drop = FALSE] %*% (values[idx]^alpha * t(vectors[, idx, drop = FALSE])))

  rownames(res) <- rownames(S)
  colnames(res) <- colnames(S)

  res
}
