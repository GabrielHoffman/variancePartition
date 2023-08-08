# Gabriel Hoffman
#
# October 13, 2020
# Simplify calculations of variance fractions and
# add compatability with glm's


#' Compute variance statistics
#'
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#'
#' @param fit model fit from lm() or lmer()
#' @param returnFractions default: TRUE.  If TRUE return fractions that sum to 1.  Else return unscaled variance components.
#' @param ... additional arguments (not currently used)
#'
#' @return
#' fraction of variance explained / ICC for each variable in the regression model
#'
#' @details
#' For linear model, variance fractions are computed based on the sum of squares explained by each component.  For the linear mixed model, the variance fractions are computed by variance component estimates for random effects and sum of squares for fixed effects.
#'
#' For a generalized linear model, the variance fraction also includes the contribution of the link function so that fractions are reported on the linear (i.e. link) scale rather than the observed (i.e. response) scale. For linear regression with an identity link, fractions are the same on both scales.  But for logit or probit links, the fractions are not well defined on the observed scale due to the transformation imposed by the link function.
#'
#' The variance implied by the link function is the variance of the corresponding distribution:
#'
#'  logit -> logistic distribution -> variance is pi^2/3
#'
#'  probit -> standard normal distribution -> variance is 1
#'
#' For the Poisson distribution with rate \eqn{\lambda}, the variance is \eqn{log(1 + 1/\lambda)}.
#'
#' For the negative binomial distribution with rate \eqn{\lambda} and shape \eqn{\theta}, the variance is \eqn{log(1 + 1/\lambda + 1/\theta)}.
#'
#' Variance decomposition is reviewed by Nakagawa and Schielzeth (2012), and expanded to other GLMs by Nakagawa, Johnson and Schielzeth (2017).  See McKelvey and Zavoina (1975) for early work on applying to GLMs.  Also see DeMaris (2002)
#'
#' We note that Nagelkerke's pseudo R^2 evaluates the variance explained by the full model.  Instead, a variance partitioning approach evaluates the variance explained by each term in the model, so that the sum of each systematic plus random term sums to 1 (Hoffman and Schadt, 2016; Nakagawa and Schielzeth, 2012).
#'
#' @references{
#'   \insertRef{nakagawa2017coefficient}{variancePartition}
#'
#'   \insertRef{nakagawa2013general}{variancePartition}
#'
#'   \insertRef{mckelvey1975statistical}{variancePartition}
#'
#'   \insertRef{demaris2002explained}{variancePartition}
#'
#'   \insertRef{hoffman2016variancepartition}{variancePartition}
#' }
#'
#' @examples
#' library(lme4)
#' data(varPartData)
#'
#' # Linear mixed model
#' fit <- lmer(geneExpr[1, ] ~ (1 | Tissue) + Age, info)
#' calcVarPart(fit)
#'
#' # Linear model
#' # Note that the two models produce slightly different results
#' # This is expected: they are different statistical estimates
#' # of the same underlying value
#' fit <- lm(geneExpr[1, ] ~ Tissue + Age, info)
#' calcVarPart(fit)
#'
#' @export
#' @docType methods
#' @rdname calcVarPart-method
setGeneric("calcVarPart",
  signature = "fit",
  function(fit, returnFractions = TRUE, ...) {
    standardGeneric("calcVarPart")
  }
)



# New version on March 30, 2021
#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lm-method
setMethod(
  "calcVarPart", "lm",
  function(fit, returnFractions = TRUE, ...) {
    # check validity of model fit
    checkModelStatus(fit, ...)

    # create design matrix
    dsgn <- model.matrix(fit$terms, fit$model)

    # loop through all groupings of variables:
    # continuous variables are 1 column,
    # categorical variables depend on the number of levels
    # i=0 indicates the intercept, but skip this
    # since it doesn't constribute to variance
    fxeff <- sapply(seq_len(max(fit$assign)), function(i) {
      idx <- which(fit$assign == i)
      dsgn[, idx, drop = FALSE] %*% fit$coefficients[idx]
    })
    colnames(fxeff) <- attr(fit$terms, "term.labels")

    # get weights
    w <- weights(fit)
    if (is.null(w)) {
      w <- rep(1, nrow(fit$model))
    }

    # get sum of squares explained by each variable
    SS <- apply(fxeff, 2, function(x) {
      weighted.var(x, w) * (length(x) - 1)
    },
    simplify = FALSE
    )

    # Compute residual sum of squares
    SS["Residuals"] <- sigma(fit)^2 * rdf(fit)

    SS <- unlist(SS)

    if (returnFractions) {
      # get variance fractions by dividing each SS by total sum of squares
      res <- SS / sum(SS)
    } else {
      res <- SS / nrow(fit$model)
    }

    res
  }
)




# from modi::weighted.var()
#' @importFrom stats weighted.mean
weighted.var <- function(x, w, na.rm = FALSE) {
  if (missing(w)) {
    w <- rep.int(1, length(x))
  } else if (length(w) != length(x)) {
    stop("x and w must have the same length")
  }
  if (min(w) < 0) {
    stop("there are negative weights")
  }
  if (is.integer(w)) {
    w <- as.numeric(w)
  }
  if (na.rm) {
    w <- w[obs.ind <- !is.na(x)]
    x <- x[obs.ind]
  }
  w <- w * length(w) / sum(w)
  return(sum(w * (x - weighted.mean(x, w))^2) / (sum(w) - 1))
}


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lmerMod-method
setMethod(
  "calcVarPart", "lmerMod",
  function(fit, returnFractions = TRUE, ...) {
    # check validity of model fit
    checkModelStatus(fit, ...)

    # extract variance components
    vc <- unlist(getVarianceComponents(fit))

    if (returnFractions) {
      # create fractions
      res <- vc / sum(vc)
    } else {
      res <- vc
    }

    # remove ".(Intercept)" string
    names(res) <- gsub("\\.\\(Intercept\\)", "", names(res))

    res
  }
)


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,glm-method
setMethod(
  "calcVarPart", "glm",
  function(fit, returnFractions = TRUE, ...) {
    checkModelStatus(fit, ...)

    cvp_glm(fit, returnFractions = returnFractions, ...)
  }
)


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,negbin-method
#' @importFrom aod negbin
setMethod(
  "calcVarPart", "negbin",
  function(fit, returnFractions = TRUE, ...) {
    checkModelStatus(fit, ...)

    cvp_glm(fit, returnFractions = returnFractions, ...)
  }
)

# Compute distribution variances for GLMs described Nakagawa, 2017
#' @importFrom stats family
getDistrVar <- function(fit) {
  # Pass BiocCheck
  link <- NA

  # compute residual term for each link
  famLink <- with(family(fit), paste(gsub("\\(.*", "", family), link))

  distVar <- switch(famLink,
    "binomial logit" = (pi^2) / 3,
    "binomial probit" = 1,
    "gaussian identity" = sigma(fit)^2 * rdf(fit),
    "poisson log" = {
      beta_0 <- coef(summary(fit))["(Intercept)", "Estimate"]
      log(1 + 1 / exp(beta_0))
    },
    "Negative Binomial log" = {
      beta_0 <- coef(summary(fit))["(Intercept)", "Estimate"]

      if (is(fit, "negbin")) theta <- fit$theta
      if (is(fit, "glmerMod")) theta <- getME(fit, "glmer.nb.theta")
      log(1 + 1 / exp(beta_0) + 1 / theta)
    }
  )

  if (is.null(distVar)) {
    stop("glm family/link not supported: ", famLink)
  }
  distVar
}

# evaluate GLM's
cvp_glm <- function(fit, returnFractions = TRUE, ...) {
  # get weights
  w <- weights(fit)
  if (is.null(w)) {
    w <- rep(1, nrow(fit$model))
  }

  # Compute eta for each term
  # predicted value in linear space for each term
  Eta <- predict(fit, type = "terms")

  # residual variance based on link function
  distVar <- getDistrVar(fit)

  # variance on linear scale
  # get variance of each term
  # append with variance due to link function
  var_term <- apply(Eta, 2, function(x) {
    weighted.var(x, w) * (length(x) - 1)
  })
  var_term <- c(var_term, Residuals = distVar)

  names(var_term) <- c(colnames(Eta), "Residuals")

  if (returnFractions) {
    # compute fraction
    res <- var_term / sum(var_term)
  } else {
    res <- var_term
  }

  res
}



#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,glmer-method
setMethod(
  "calcVarPart", "glmerMod",
  function(fit, returnFractions = TRUE, ...) {
    checkModelStatus(fit, ...)

    cvp_glmm(fit, returnFractions = returnFractions, ...)
  }
)


# evaluate GLMM's
cvp_glmm <- function(fit, returnFractions = TRUE, ...) {
  # Extract variance components
  vc <- getVarianceComponents(fit)

  # extract distribution-specific variance
  vc$Residuals <- getDistrVar(fit)

  vc <- unlist(vc)

  # remove ".(Intercept)" string
  names(vc) <- gsub("\\.\\(Intercept\\)", "", names(vc))

  if (returnFractions) {
    # create fractions
    res <- vc / sum(vc)
  } else {
    res <- vc
  }

  res
}



#' @importFrom lme4 VarCorr fixef
getVarianceComponents <- function(fit) {
  # get weights
  w <- weights(fit)
  if (is.null(w)) {
    w <- rep(1, nrow(fit$model))
  }

  # get random effects estimates
  varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)

  # order variables by name
  # essential so that all models are ordered the same
  varComp <- varComp[order(names(varComp))]

  # extract predictor for each fixed effect
  idx <- which(colnames(fit@pp$X) != "(Intercept)")

  # if there are fixed effects
  if (length(idx) > 0) {
    # this part is now in 1.5.2
    # better estimates of fixed effects
    fxeff <- sapply(idx, function(i) {
      fit@pp$X[, i] * lme4::fixef(fit)[i]
    })
    colnames(fxeff) <- colnames(fit@pp$X)[idx]

    # compute variance of each fixed effect
    # variance of each term
    N <- nrow(fxeff)

    fixedVar <- apply(fxeff, 2, function(x) {
      weighted.var(x, w) * (N - 1) / N
    })

    varFixedTotal <- weighted.var(rowSums(fxeff), w) * (N - 1) / N

    for (key in names(fixedVar)) {
      varComp[[key]] <- as.array(fixedVar[[key]] / sum(fixedVar) *
        varFixedTotal)
    }
  }

  # get residuals
  varComp$Residuals <- sigma(fit)^2

  return(varComp)
}
