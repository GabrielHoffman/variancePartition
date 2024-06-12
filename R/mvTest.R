# Gabriel Hoffman
# October 4, 2022
#
# Perform multivariate tests on results from dream() using vcov()






#' Multivariate tests on results from \code{dream()}
#'
#' Evaluate multivariate tests on results from \code{dream()} using \code{vcov()} to compute the covariance between estimated regression coefficients across multiple responses.  A joint test to see if the coefficients are jointly different from zero is performed using meta-analysis methods that account for the covariance.
#'
#' @param fit \code{MArrayLM} or \code{MArrayLM2} returned by \code{dream()}
#' @param vobj matrix or \code{EList} object returned by \code{voom()}
#' @param features a) indeces or names of features to perform multivariate test on, b) list of indeces or names.  If missing, perform joint test on all features.
#' @param coef name of coefficient or contrast to be tested
#' @param method statistical method used to perform multivariate test.  See details.  \code{'FE'} is a fixed effect test that models the covariance between coefficients.  \code{'FE.empirical'} use compute empirical p-values by sampling from the null distribution and fitting with a gamma. \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'hotelling'} performs the Hotelling T2 test. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#' @param shrink.cov shrink the covariance matrix between coefficients using the Schafer-Strimmer method
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other arugments
#'
#' @details See package \code{remaCor} for details about the \code{remaCor::RE2C()} test, and see \code{remaCor::LS()} for details about the fixed effect test.  When only 1 feature is selected, the original p-value is returned and the test statistic is set to \code{NA}.
#'
#' For the \code{"RE2C"} test, the final test statistic is the sum of a test statistic for the mean effect (\code{stat.FE}) and heterogeneity across effects (\code{stat.het}).  \code{mvTest()} returns 0 if \code{stat.het} is negative in extremely rare cases.
#'
#' @return
#' Returns a \code{data.frame} with the statistics from each test, the \code{pvalue} from the test, \code{n_features},  \code{method}, and \code{lambda} from the Schafer-Strimmer method to shrink the estimated covariance.  When \code{shrink.cov=FALSE}, \code{lambda = 0}.
#'
#' @examples
#' # library(variancePartition)
#' library(edgeR)
#' library(BiocParallel)
#'
#' data(varPartDEdata)
#'
#' # normalize RNA-seq counts
#' dge <- DGEList(counts = countMatrix)
#' dge <- calcNormFactors(dge)
#'
#' # specify formula with random effect for Individual
#' form <- ~ Disease + (1 | Individual)
#'
#' # compute observation weights
#' vobj <- voomWithDreamWeights(dge[1:20, ], form, metadata)
#'
#' # fit dream model
#' fit <- dream(vobj, form, metadata)
#' fit <- eBayes(fit)
#'
#' # Multivariate test of features 1 and 2
#' mvTest(fit, vobj, 1:2, coef = "Disease1")
#'
#' # Test multiple sets of features
#' lst <- list(a = 1:2, b = 3:4)
#' mvTest(fit, vobj, lst, coef = "Disease1", BPPARAM = SnowParam(2))
#' @export
#' @docType methods
#' @rdname mvTest-method
setGeneric("mvTest",
  signature = c("fit", "vobj", "features"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, BPPARAM = SerialParam(), ...) {
    standardGeneric("mvTest")
  }
)


#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,integer-method
#' @importFrom remaCor RE2C LS LS.empirical hotelling
#' @importFrom stats coefficients pchisq cov2cor
#' @importFrom corpcor estimate.lambda
#' @export
setMethod(
  "mvTest", c("MArrayLM", "EList", "vector"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, BPPARAM = SerialParam(), ...) {
    method <- match.arg(method)

    i <- match(coef, colnames(coefficients(fit)))
    if (any(is.na(i))) {
      txt <- paste("Coefficients not valid:", paste(coef[is.na(i)], collapse = ", "))
      stop(txt)
    }

    if (missing(features)) {
      features <- rownames(fit)
    } else {
      # only test features that are characters
      feat.test <- features[is.character(features)]
      i <- match(feat.test, rownames(fit))
      if (any(is.na(i))) {
        txt <- paste("Features not found:", paste(feat.test[is.na(i)], collapse = ", "))
        stop(txt)
      }
    }

    # extract coefficients from features
    tab <- topTable(fit[features, ], coef = coef, sort.by = "none", number = Inf)
    beta <- tab$logFC

    n_features <- length(features)

    # Residual degrees of freedom
    # since mixed model estimates rdf per-gene
    # use minmum to be conservative
    if (is(fit, "MArrayLM2")) {
      # mixed model
      nu <- min(fit[features, ]$rdf)
    } else {
      # fixed effect model
      nu <- min(fit$df.residual)
    }

    # extract covariance
    #-------------------
    # get Sigma directly
    # Sigma = vcov(fit[features,], vobj[features,], coef)

    # Instead, get sqrt of covariance so that Sigma is crossprod(P)
    # Then estimate shrinkage intensity
    # and shrink covariance
    # this is important when p approaches n
    # Note that the test below does not model uncertainy in lambda
    P <- vcovSqrt(fit[features, ], vobj[features, ], coef, approx = TRUE)
    Sigma <- crossprod(P)

    if( shrink.cov ){
      lambda = 0.01

      Sigma <- (1 - lambda) * Sigma + lambda * diag(diag(Sigma), ncol(Sigma))
    }else{
      lambda = 0
    }

    # Meta-analyis method
    #####################

    if (method == "FE.empirical") {
      if (n_features == 1) {
        # for one test, return estimated t-stat as stat
        df <- data.frame(
          beta = beta,
          se = Sigma[1,1],
          stat = tab$t,
          pvalue = tab$P.Value,
          n_features = 1,
          lambda = lambda,
          method = method
        )
      } else {
        # Ensures that sampling from Wishart is defined
        nu <- max(nu, n_features)

        res <- LS.empirical(beta, sqrt(diag(Sigma)), cov2cor(Sigma), nu, ...)

        df <- data.frame(          
          beta = res$beta,
          se = res$se,
          stat = res$beta / res$se,
          pvalue = res$p,
          n_features = n_features,
          lambda = lambda,
          method = method
        )
      }
    } else if (method == "FE") {
      if (n_features == 1) {
        # for one test, return estimated t-stat as stat
        df <- data.frame(
          beta = beta,
          se = Sigma[1,1],
          stat = tab$t,
          pvalue = tab$P.Value,
          n_features = 1,
          lambda = lambda,
          method = method
        )
      } else {
        # Use asymptotic normal null distribution
        res <- LS(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

        df <- data.frame(
          beta = res$beta,
          se = res$se,
          stat = res$beta / res$se,
          pvalue = res$p,
          n_features = n_features,
          lambda = lambda,
          method = method
        )
      }
    } else if (method == "RE2C") {
      if (n_features == 1) {
        # for one test, heterogeneity is zero
        df <- data.frame(
          stat.FE = tab$t^2,
          stat.het = 0,
          pvalue = tab$P.Value,
          n_features = 1,
          lambda = lambda,
          method = method
        )
      } else {
        res <- RE2C(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

        df <- data.frame(
          stat.FE = res$stat1,
          stat.het = max(0, res$stat2),
          pvalue = res$RE2Cp,
          n_features = n_features,
          lambda = lambda,
          method = method
        )
      }
    } else if (method == "tstat") {
      if (n_features == 1) {
        df <- data.frame(
          stat = tab$t,
          pvalue = tab$P.Value,
          n_features = n_features,
          lambda = lambda,
          method = method
        )
      } else {
        Sigma_corr <- cov2cor(Sigma)
        # tstat = crossprod(tab$t, solve(Sigma_corr, tab$t))

        # in case Sigma_corr is not invertable
        tstat <- tryCatch(crossprod(tab$t, solve(Sigma_corr, tab$t)), error = function(e) {
          warning("Covariance matrix is not invertable. Returning NA values.")
          NA
        })

        pv <- pchisq(tstat, length(tab$t), lower.tail = FALSE)

        df <- data.frame(
          stat = tstat,
          pvalue = pv,
          n_features = n_features,
          lambda = lambda,
          method = method
        )
      }
    } else if (method == "hotelling") {
      # in case Sigma_corr is not invertable
      res <- tryCatch(hotelling(beta, Sigma, nu),
        error = function(e) {
          warning("Covariance matrix is not invertable. Returning NA values.")
          NA
        }
      )

      df <- data.frame(
        stat = res$stat,
        pvalue = res$p.value,
        n_features = n_features,
        lambda = lambda,
        method = method
      )
    } else if (method == "sidak") {
      pv <- 1 - (1 - min(tab$P.Value))^nrow(tab)

      df <- data.frame(
        stat = NA,
        pvalue = pv,
        n_features = n_features,
        lambda = lambda,
        method = method
      )
    } else if (method == "fisher") {
      stat <- -2 * sum(log(tab$P.Value))
      k <- nrow(tab)
      pv <- pchisq(stat, 2 * k, lower.tail = FALSE)

      df <- data.frame(
        stat = stat,
        pvalue = pv,
        n_features = n_features,
        lambda = lambda,
        method = method
      )
    }

    df
  }
)

# if no features are specified, test all features
#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,missing-method
#' @export
setMethod(
  "mvTest", c("MArrayLM", "EList", "missing"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, BPPARAM = SerialParam(), ...) {
    mvTest(fit, vobj, features = seq(nrow(vobj)), coef, method, shrink.cov = shrink.cov, BPPARAM, ...)
  }
)



#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,list-method
#' @importFrom stats runif
#' @export
setMethod(
  "mvTest", c("MArrayLM", "EList", "list"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, BPPARAM = SerialParam(), ...) {
    if (is.null(names(features))) {
      stop("features list must have non-null names(features)")
    }

    it <- iterRowsSplit(vobj$E, vobj$weights, fit, splitList = features)

    # BPPARAM$exportglobals <- FALSE
    res <- bpiterate(it, mvTest,
      coef = coef,
      method = method,
      shrink.cov = shrink.cov,
      BPPARAM = BPPARAM
    )

    do.call(rbind, res)
  }
)

#' Class mvTest_input
#'
#' Class \code{mvTest_input} work is with \code{iterRowsSplit()}
#'
#' @name mvTest_input-class
#' @rdname mvTest_input-class
#' @exportClass mvTest_input
setClass("mvTest_input", representation(vobj = "EList", fit = "MArrayLM", ID = "character"))

#' @rdname mvTest-method
#' @aliases mvTest,mvTest_input,method
#' @export
setMethod(
  "mvTest", c("mvTest_input"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, ...) {
    res <- mvTest(fit@fit, fit@vobj,
      coef = coef,
      method = method,
      shrink.cov = shrink.cov,
      ...
    )

    data.frame(ID = fit@ID, res)
  }
)



#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,matrix-method
#' @export
setMethod(
  "mvTest", c("MArrayLM", "matrix"),
  function(fit, vobj, features, coef, method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"), shrink.cov = TRUE, ...) {
    method <- match.arg(method)

    # create weights
    W <- matrix(1, nrow(vobj), ncol(vobj))

    vobj <- new("EList", list(E = vobj, weights = W))

    mvTest(fit, vobj, features, coef, method, shrink.cov = shrink.cov, , ...)
  }
)
