# Gabriel Hoffman
# December 13, 2020
#
# Compute the residual degrees of freedom for a linear mixed model

#' Residual degrees of freedom
#'
#' Residual degrees of freedom
#'
#' @param fit model fit from \code{lm()}, \code{glm()}, \code{lmer()}
#'
#' @examples
#' library(lme4)
#'
#' fit <- lm(Reaction ~ Days, sleepstudy)
#' rdf(fit)
#'
#' @export
#' @seealso \code{rdf.merMod}
rdf <- function(fit) {
  if (is(fit, "lm")) {
    return(fit$df.residual)
  }
  if (is(fit, "glm")) {
    return(fit$df.residual)
  }
  if (is(fit, "merMod")) {
    return(rdf.merMod(fit))
  }
  NA
}


# trace only gives correct rdf if evals are 1 and 0
# trace is given in the book with H is idempotent.
# I = diag(1, n)
# tr((I-H) %*% (I-H))
# When H is not idempotent, this actually gives the degrees of freedom
# of the best chi-square fit to the distribution
# if x ~ sum(lambda * rchisq(n)), then best chi-square approximation
# is a 1 momoent match to the mean of sum(lambda)
# When no strinkage, rdf.trace is exact.
# The approximation improves as sum(lambda) increases
# For smaller rdf, the chi-square is shifted slighly left
# rdf.naive = function(fit){
#   tr = function(A) sum(diag(A))
#    # get hat matrix from linear mixed model
#   H = lme4:::hatvalues.merMod(fit, fullHatMatrix=TRUE)

#   # number of samples
#   n = nrow(H)

#   # I = diag(1, n)
#   # tr((I-H) %*% (I-H))
#   n - 2*tr(H) + sum(H*H)
# }

#' Fast approximate residual degrees of freedom
#'
#' Defining \eqn{H = A^TA + B^TB} where \eqn{A} and \eqn{B} are low rank, compute
#' \eqn{n - 2tr(H) + tr(HH)} in \eqn{O(np^2)} instead of \eqn{O(n^2p^2)}.
#'
#' @param A a \code{matrix} or \code{sparseMatrix}
#' @param B a \code{matrix} or \code{sparseMatrix}
#'
# @examples
#
# # sample size
# n = 500
# # number of covariates
# n_a = 20
# n_b = 20
#
# # Simulate low rank matrices
# A = matrix(rnorm(n_a*n), ncol=n)
# B = matrix(rnorm(n_b*n), ncol=n)
#
# # Evaluate RDF
# rdf_from_matrices(A, B)
#
#' @importFrom Matrix drop0
#' @seealso rdf.merMod
#'
rdf_from_matrices <- function(A, B) {
  # desired result, but this simple code is O(n^2)
  # H = crossprod(A) + crossprod(B)
  # n = nrow(H)
  # rdf = n - 2*tr(H) + sum(H*H)

  # Compute spectral decomposition  of sparse matrix using eigen decomposition of crossproduct

  # Pass BiocCheck
  u <- d <- NA

  # A
  #####

  # Even though this is a eigen decomp of a crossprod
  # the eigen-values can be very slightly negative.
  # catch this error and perform SVD instead
  dcmp_A <- tryCatch(
    {
      eigen(tcrossprod(A))
    },
    error = function(e) {
      with(svd(A, nv = 0), list(vectors = u, values = d^2))
    }
  )

  # drop eigen-values less than tol
  tol <- 1e-12
  idx <- which(dcmp_A$values > tol)
  U <- dcmp_A$vectors[, idx, drop = FALSE]
  s_a <- sqrt(dcmp_A$values[idx])
  # VT = crossprod(U, A) / s_a

  # B
  #####
  dcmp_B <- eigen(tcrossprod(B), only.values = TRUE)
  s_b <- sqrt(dcmp_B$values)

  G <- as.matrix(tcrossprod(B, A) %*% U)

  tr_H_H <- sum(s_a^4) +
    2 * sum(G * G) +
    sum(s_b^4)

  n <- ncol(A)

  n - 2 * (sum(s_a^2) + sum(s_b^2)) + tr_H_H
}




#' Approximate residual degrees of freedom
#'
#' Compute the approximate residual degrees of freedom from a linear mixed model.
#'
#' @param model An object of class \code{merMod}
#' @param method Use algorithm that is "linear" (default) or quadratic time in the number of samples
#'
#' @description
#' For a linear model with \eqn{n} samples and \eqn{p} covariates, \eqn{RSS/\sigma^2 \sim \chi^2_{\nu}} where \eqn{\nu = n-p} is the residual degrees of freedom.  In the case of a linear mixed model, the distribution is no longer exactly a chi-square distribution, but can be approximated with a chi-square distribution.
#'
#' Given the hat matrix, \eqn{H}, that maps between observed and fitted responses, the approximate residual degrees of freedom is \eqn{\nu = tr((I-H)^T(I-H))}.  For a linear model, this simplifies to the well known form \eqn{\nu = n - p}. In the more general case, such as a linear mixed model, the original form simplifies only to \eqn{n - 2tr(H) + tr(HH)} and is an approximation rather than being exact.  The third term here is quadratic time in the number of samples, \eqn{n}, and can be computationally expensive to evaluate for larger datasets.  Here we develop a linear time algorithm that takes advantage of the fact that \eqn{H} is low rank.
#'
#' \eqn{H} is computed as \eqn{A^TA + B^TB} for \code{A=CL} and \code{B=CR} defined in the code.  Since \eqn{A} and \eqn{B} are low rank, there is no need to compute \eqn{H} directly.  Instead, the terms \eqn{tr(H)} and \eqn{tr(HH)} can be computed using the eigen decompositions of \eqn{AA^T} and \eqn{BB^T} which is linear time in the number of samples.
#'
#' @return residual degrees of freedom
#' @examples
#' library(lme4)
#'
#' # Fit linear mixed model
#' fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#'
#' # Evaluate the approximate residual degrees of freedom
#' rdf.merMod(fit)
#'
#' @importFrom Matrix solve Diagonal diag crossprod tcrossprod t
#' @importFrom stats weights hatvalues
#' @importFrom lme4 isGLMM getME
#' @seealso rdf_from_matrices
# @method rdf merMod
#' @export
rdf.merMod <- function(model, method = c("linear", "quadratic")) {
  method <- match.arg(method)

  # code adapted from lme4:::hatvalues.merMod
  if (isGLMM(model)) {
    warning("the hat matrix may not make sense for GLMMs")
  }

  if (method == "linear") {
    sqrtW <- Diagonal(x = sqrt(weights(model, type = "prior")))

    # Pass Bioconductor check
    L <- Lambdat <- Zt <- RX <- X <- RZX <- NULL

    rdf <- with(getME(model, c("L", "Lambdat", "Zt", "RX", "X", "RZX")), {
      CL <- solve(L, solve(L, Lambdat %*% Zt %*% sqrtW, system = "P"), system = "L")
      CR <- solve(t(RX), t(X) %*% sqrtW - crossprod(RZX, CL))
      # if (fullHatMatrix)
      #     crossprod(CL) + crossprod(CR)
      # else colSums(CR^2) + colSums(CL^2)

      # H = crossprod(CL) + crossprod(CR)
      # compute residual degrees of freedom as if using H, but use only CL and CR
      rdf_from_matrices(CL, CR)
    })
  } else {
    H <- hatvalues(model, fullHatMatrix = TRUE)

    # number of samples
    n <- nrow(H)

    # Simple code
    # I = diag(1, n)
    # tr((I-H) %*% (I-H))
    rdf <- n - 2 * sum(diag(H)) + sum(H * H)
  }
  rdf
}


#' Shrinkage metric for eBayes
#'
#' Shrinkage metric for eBayes quantifying the amount of shrinkage that is applied to shrink the maximum likelihood residual variance to the empirical Bayes posterior estimate
#'
#' @param sigmaSq maximum likelihood residual variance for every gene
#' @param s2.post empirical Bayes posterior estimate of residual variance for every gene
#'
#' @description Evaluates the coefficient from the linear regression of \code{s2.post ~ sigmaSq}. When there is no shrinkage, this value is 1.  Values less than 1 indicate the amount of shrinkage.
#'
#' @importFrom stats lm coef
#' @export
shrinkageMetric <- function(sigmaSq, s2.post) {
  fit <- lm(s2.post ~ sigmaSq)
  as.numeric(coef(fit)[2])
}



#' Scaled chi-square
#'
#' Scaled chi-square density using a gamma distribution
#'
#' @param x vector of quantiles.
#' @param a scale
#' @param b degrees of freedom
#'
#' @importFrom stats dgamma
dscchisq <- function(x, a, b) {
  dgamma(x, b / 2, 1 / (2 * a))
}


#' Plot Variance Estimates
#'
#' @param fit model fit from \code{dream()}
#' @param fitEB model fit from \code{eBayes()}
#' @param var_true array of true variance values from simulation (optional)
#' @param xmax maximum value on the x-axis
#'
#' @importFrom stats density
#' @import ggplot2
#' @export
plotVarianceEstimates <- function(fit, fitEB, var_true = NULL, xmax = quantile(fit$sigma^2, .999)) {
  # Pass R CMD check
  Method <- y <- NA

  # x values where to evaluate the scaled chi-square density
  x <- seq(1e-3, xmax, length.out = 1000)

  # MLE
  d_mle <- density(fit$sigma^2, from = 0, to = xmax)
  df_combine <- data.frame(Method = "MLE", x = d_mle$x, y = d_mle$y)

  # EB posterior
  if (var(fitEB$s2.post) > 0) {
    d_posterior <- density(fitEB$s2.post, from = 0, to = xmax)
    df_combine <- rbind(
      df_combine,
      data.frame(Method = "EB posterior", x = d_posterior$x, y = d_posterior$y)
    )
  } else {
    ymax <- max(df_combine$y)
    df_combine <- rbind(df_combine, data.frame(
      Method = "EB posterior",
      x = c(fitEB$s2.post[1] - .01, fitEB$s2.post[1], fitEB$s2.post[1] + 0.01),
      y = c(0, ymax, 0)
    ))
  }


  # define colors
  col <- c("MLE" = "red", "EB prior" = "orange", "EB posterior" = "blue")

  # True variance
  if (!is.null(var_true)) {
    col <- c("True variance" = "green", col)
    if (var(var_true) > 0) {
      d_true <- density(var_true, from = 0, to = xmax)
      df_combine <- rbind(
        df_combine,
        data.frame(Method = "True variance", x = d_true$x, y = d_true$y)
      )
    } else {
      ymax <- max(df_combine$y)
      df_combine <- rbind(df_combine, data.frame(
        Method = "True variance",
        x = c(var_true[1] - .01, var_true[1], var_true[1] + 0.01),
        y = c(0, ymax, 0)
      ))
    }
  }

  # compute prior density even when eBayes() uses trend=TRUE
  # so that the s2.prior is a vector

  # these can be scalars or vectors (depending on robust and trend)
  # if scalar, use rep() to create a vector
  df.prior <- fitEB$df.prior
  s2.prior <- fitEB$s2.prior
  n_genes <- length(fitEB$s2.post)

  if (length(df.prior) == 1) {
    df.prior <- rep(df.prior, n_genes)
  }
  if (length(s2.prior) == 1) {
    s2.prior <- rep(s2.prior, n_genes)
  }

  dst <- sapply(seq_len(n_genes), function(i) {
    dscchisq(x, (s2.prior[i] / df.prior[i]), df.prior[i])
  })
  scale_chiSq_density <- rowSums(dst) / ncol(dst)

  # Empirical Bayes prior density
  # if df.prior is finite
  if (all(is.finite(fitEB$df.prior))) {
    df_combine <- rbind(df_combine, data.frame(
      Method = "EB prior", x = x,
      y = scale_chiSq_density
    ))
  } else {
    # if df.prior is infinite, this is a point mass
    ymax <- max(df_combine$y)
    df_combine <- rbind(df_combine, data.frame(
      Method = "EB prior",
      x = c(fitEB$s2.prior - .01, fitEB$s2.prior, fitEB$s2.prior + 0.01),
      y = c(0, ymax, 0)
    ))
  }

  # order methods
  df_combine$Method <- factor(df_combine$Method, c("True variance", "MLE", "EB prior", "EB posterior"))
  df_combine$Method <- droplevels(df_combine$Method)

  # plot
  ymax <- max(df_combine$y) * 1.05
  ggplot(df_combine, aes(x, y, color = Method)) +
    geom_line() +
    theme_bw(16) +
    theme(legend.position = "right", aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    scale_color_manual(values = col[levels(df_combine$Method)]) +
    xlab(bquote(hat(sigma)^2)) +
    ylab("Density") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, xmax)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, ymax))
}
