#' Differential expression with linear mixed model
#'
#' Fit linear mixed model for differential expression and preform hypothesis test on fixed effects as specified in the contrast matrix \code{L}
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.
#' @param data data.frame with columns corresponding to formula
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough approximation for most datasets. "adaptive" (Default) uses KR for <= 20 samples.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param control control settings for \code{lmer()}
#' @param hideErrorsInBackend default FALSE.  If TRUE, hide errors in \code{attr(.,"errors")} and \code{attr(.,"error.initial")}
#' @param BPPARAM parameters for parallel evaluation
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is TRUE.  See Details.
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#'
#' @return
#' MArrayLM2 object (just like MArrayLM from limma), and the directly estimated p-value (without eBayes)
#'
#' @details
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression (Hoffman and Roussos, 2021).  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is
#'
#' \code{fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)}
#'
#' \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run code in parallel.  Parallel processing must be enabled before calling this function.  See below.
#'
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lmer.
#'
#' Hypothesis tests and degrees of freedom are producted by \code{lmerTest} and \code{pbkrtest} pacakges
#'
#' While \code{REML=TRUE} is required by \code{lmerTest} when ddf='Kenward-Roger', ddf='Satterthwaite' can be used with \code{REML} as \code{TRUE} or \code{FALSE}.  Since the Kenward-Roger method gave the best power with an accurate control of false positive rate in our simulations, and since the Satterthwaite method with REML=TRUE gives p-values that are slightly closer to the Kenward-Roger p-values, \code{REML=TRUE} is the default.  See Vignette "3) Theory and practice of random effects and REML"
#'
#' @references{
#'   \insertRef{hoffman2021dream}{variancePartition}
#' }
#' @examples
#' # library(variancePartition)
#'
#' # load simulated data:
#' # geneExpr: matrix of *normalized* gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#'
#' form <- ~ Batch + (1 | Individual) + (1 | Tissue)
#'
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' # NOTE: dream() runs on *normalized* data
#' fit <- dream(geneExpr[1:10, ], form, info)
#' fit <- eBayes(fit)
#'
#' # view top genes
#' topTable(fit, coef = "Batch2", number = 3)
#'
#' # get contrast matrix testing if the coefficient for Batch3 is
#' # different from coefficient for Batch2
#' # Name this comparison as 'compare_3_2'
#' # The variable of interest must be a fixed effect
#' L <- makeContrastsDream(form, info, contrasts = c(compare_3_2 = "Batch3 - Batch2"))
#'
#' # plot contrasts
#' plotContrasts(L)
#'
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit2 <- dream(geneExpr[1:10, ], form, info, L)
#' fit2 <- eBayes(fit2)
#'
#' # view top genes for this contrast
#' topTable(fit2, coef = "compare_3_2", number = 3)
#'
#' # Parallel processing using multiple cores with reduced memory usage
#' param <- SnowParam(4, "SOCK", progressbar = TRUE)
#' fit3 <- dream(geneExpr[1:10, ], form, info, L, BPPARAM = param)
#' fit3 <- eBayes(fit3)
#'
#' # Fit fixed effect model for each gene
#' # Use lmFit in the backend
#' form <- ~Batch
#' fit4 <- dream(geneExpr[1:10, ], form, info, L)
#' fit4 <- eBayes(fit4)
#'
#' # view top genes
#' topTable(fit4, coef = "compare_3_2", number = 3)
#'
#' # Compute residuals using dream
#' residuals(fit4)[1:4, 1:4]
#'
#' @export
# @docType methods
#' @rdname dream-method
#' @importFrom pbkrtest get_SigmaG
#' @importFrom BiocParallel bpiterate bpok SerialParam
#' @importFrom lme4 VarCorr
#' @importFrom stats hatvalues as.formula
#' @importFrom methods as
dream <- function(exprObj,
                  formula,
                  data,
                  L,
                  ddf = c("adaptive", "Satterthwaite", "Kenward-Roger"),
                  useWeights = TRUE,
                  control = vpcontrol,
                  hideErrorsInBackend = FALSE,
                  REML = TRUE,
                  BPPARAM = SerialParam(),
                  ...) {
  ddf <- match.arg(ddf)

  # if ddf is "adaptive" and # samples < 20, use "Satterthwaite"
  if( ddf == "adaptive" ){
    ddf = ifelse( ncol(exprObj) < 20, "Kenward-Roger", "Satterthwaite")
  }

  # filter and check input data
  objFlt <- filterInputData(exprObj, formula, data, useWeights = useWeights)

  # process contrasts to create L matrix
  ctst <- createContrastL(objFlt$formula, objFlt$data, L)

  if (.isMixedModelFormula(formula)) {
    # run model
    res <- run_lmm(objFlt$exprObj, objFlt$formula, objFlt$data,
      fxn = create_eval_dream(L = ctst$L, ddf = ddf, ctst$ univariateContrasts),
      REML = REML,
      control = control,
      BPPARAM = BPPARAM,
      dreamCheck = TRUE,
      useInitialFit = FALSE,
      # rescaleWeights = FALSE,
      ...
    )

    if (!is.null(res$error.initial) & !hideErrorsInBackend) {
      stop(paste0("Initial model failed:\n", res$error.initial))
    }

    res.unlist <- unlist(res$succeeded, recursive = FALSE)

    if (length(res.unlist) == 0 & !hideErrorsInBackend) {
      txt <- paste("All models failed.  The first model fails with:\n", res$errors[1])
      stop(txt)
    }

    res2 <- combineResults(objFlt$exprObj, ctst$L, res.unlist, ctst$univariateContrasts)

    attr(res2, "error.initial") <- res$error.initial
    attr(res2, "errors") <- res$errors
  } else {
    # run limma fixed effects models
    design <- model.matrix(objFlt$formula, objFlt$data)

    if (kappa(design) > 1e7) {
      stop("Design matrix is singular, covariates are very correlated")
    }

    res2 <- lmFit(objFlt$exprObj, design)
    res2$residuals <- residuals(res2, objFlt$exprObj)
    res2$logLik = logLik(res2, objFlt$exprObj)
    res2$BIC = BIC(res2, objFlt$exprObj)

    # apply contrasts
    if (!ctst$univariateContrasts) {
      res2 <- contrasts.fit(res2, ctst$L)
    }

    # compute hat values
    res2$hatvalues <- hatvalues(res2, objFlt$exprObj)
  }

  res2$formula <- objFlt$formula
  res2$data <- objFlt$data

  if (!is.null(attr(res2, "errors")) & !hideErrorsInBackend) {
    txt <- paste("Model failed for", length(attr(res2, "errors")), "responses.\n  See errors with attr(., 'errors')")
    warning(txt)
  }

  res2
}


#' @importFrom lme4 VarCorr refit
#' @importFrom stats sigma logLik
create_eval_dream <- function(L, ddf, univariateContrasts) {
  function(x) {
    # convert to result of lmerTest::lmer()
    fit <- as_lmerModLmerTest2(x)

    # check L
    if (!identical(rownames(L), names(fixef(fit)))) {
      stop("Names of entries in L must match fixed effects")
    }

    # check that fixed effects and names of L match
    a <- names(fixef(fit))
    b <- rownames(L)

    if (!identical(a, b)) {
      stop("Terms in contrast matrix L do not match model:\n  Model: ", paste(a, collapse = ","), "\n  L: ", paste(b, collapse = ","), "\nNote thhat order must be the same")
    }

    mod <- eval_contrasts(fit, L, ddf)

    ret <- list(
      coefficients = mod$beta,
      design = fit@pp$X,
      # approximate df of residuals
      rdf = rdf.merMod(fit),
      # df of test statistics
      df.residual = mod$df,
      Amean = mean(fit@frame[, 1], na.rm = TRUE),
      method = "lmer",
      sigma = mod$sigma,
      stdev.unscaled = mod$SE / mod$sigma,
      pValue = mod$pValue,
      residuals = residuals(fit)
    )

    # get variance terms for random effects
    varComp <- lapply(VarCorr(fit), function(x) {
      attr(x, "stddev")^2
    })
    varComp[["resid"]] <- sigma(fit)^2

    if (univariateContrasts) {
      V <- mod$vcov
    } else {
      # do expensive evaluation only when L is defined by the user
      V <- crossprod(chol(mod$vcov) %*% L)
    }

    # compute hatvalues
    h <- hatvalues(fit)

    list(
      ret = new("MArrayLM", ret),
      varComp = varComp,
      logLik = as.numeric(logLik(fit)),
      # effective degrees of freedom as sum of diagonals of hat matrix
      edf = sum(h),
      hatvalues = h,
      vcov = V
    )
  }
}

#' @importFrom lmerTest as_lmerModLmerTest contest
#' @importFrom pbkrtest vcovAdj.lmerMod
eval_contrasts <- function(fit, L, ddf, kappa.tol = 1e6, pd.tol = 1e-8) {
  # convert to lmerModLmerTest object
  if (!is(fit, "lmerModLmerTest")) {
    fit <- as_lmerModLmerTest2(fit)
  }

  # determine ddf method
  # if KR method is specified but gives invalid covariance matrix
  # fall back on Satterthwaite
  # compute variance-covariance matrix
  if (ddf == "Kenward-Roger") {
    # get vcov using KR adjusted method
    V <- as.matrix(vcovAdj.lmerMod(fit, 0))

    # if poor condition number, or not positive definite
    # fall back on Satterthwaite
    if (kappa(V) > kappa.tol | !isPositiveDefinite(V, pd.tol)) {
      ddf <- "Satterthwaite"
    }
  }

  if (ddf == "Satterthwaite") {
    # standard
    V <- as.matrix(vcov(fit))
  }

  cons <- contest(fit, t(L), ddf = ddf, joint = FALSE, confint = FALSE)

  # extract results
  df <- cons$df
  pValue <- as.numeric(cons[, "Pr(>|t|)"])

  beta <- matrix(cons$Estimate, ncol = 1)
  colnames(beta) <- "logFC"
  rownames(beta) <- colnames(L)

  SE <- matrix(cons[["Std. Error"]], ncol = 1)
  colnames(SE) <- "logFC"
  rownames(SE) <- colnames(L)

  list(
    cons = cons,
    df = df,
    sigma = sigma(fit),
    beta = beta,
    SE = SE,
    pValue = pValue,
    vcov = V
  )
}



#' @importFrom stats getCall formula
as_lmerModLmerTest2 <- function(model, tol = 1e-08) {
  if (!inherits(model, "lmerMod")) {
    stop("model not of class 'lmerMod': cannot coerce to class 'lmerModLmerTest")
  }
  mc <- getCall(model)
  args <- c(as.list(mc), devFunOnly = TRUE)
  if (!"control" %in% names(as.list(mc))) {
    args$control <- lme4::lmerControl(check.rankX = "silent.drop.cols")
  }
  Call <- as.call(c(list(quote(lme4::lmer)), args[-1]))
  ff <- environment(formula(model))
  pf <- parent.frame(n = 2) # GEH: n=2 instead of 1
  sf <- sys.frames()[[1]]
  ff2 <- environment(model)
  devfun <- tryCatch(eval(Call, envir = pf), error = function(e) {
    tryCatch(eval(Call, envir = ff), error = function(e) {
      tryCatch(eval(Call, envir = ff2), error = function(e) {
        tryCatch(eval(Call, envir = sf), error = function(e) {
          "error"
        })
      })
    })
  })
  if ((is.character(devfun) && devfun == "error") || !is.function(devfun) ||
    names(formals(devfun)[1]) != "theta") {
    stop("Unable to extract deviance function from model fit")
  }
  lmerTest:::as_lmerModLT(model, devfun, tol = tol)
}

createContrastL <- function(formula, data, L) {
  univariateContrasts <- FALSE
  if (missing(L) || is.null(L)) {
    # all univariate contrasts
    L <- .getAllUniContrasts(formula, data)
    univariateContrasts <- TRUE
  } else {
    # format contrasts
    if (is(L, "numeric")) {
      L <- as.matrix(L, ncol = 1)
    } else if (is(L, "data.frame")) {
      L <- as.matrix(L)
    }
    if (is.null(colnames(L))) {
      colnames(L) <- paste0("L", seq_len(ncol(L)))
    }

    # check columns that only have a single 1
    tst <- apply(L, 2, function(x) {
      length(x[x != 0]) == 1
    })
    if (any(tst)) {
      warning("Contrasts with only a single non-zero term are already evaluated by default.")
    }
    # remove univariate contrasts
    # L = L[,!tst,drop=FALSE]

    # add all univariate contrasts
    Luni <- .getAllUniContrasts(formula, data)
    L <- cbind(L, Luni)
  }

  if (ncol(L) == 0) {
    stop("Must include fixed effect in the model for hypothesis testing")
  }

  # check rownames of contrasts
  if (length(unique(colnames(L))) != ncol(L)) {
    stop(paste("Contrast names must be unique: ", paste(colnames(L), collapse = ", ")))
  }

  list(L = L, univariateContrasts = univariateContrasts)
}


#' @importFrom methods as
combineResults <- function(exprObj, L, resList, univariateContrasts) {

  # only keep expression data from genes in resList
  exprObj <- exprObj[names(resList), seq(ncol(exprObj)), drop = FALSE]

  if (is.null(resList) | length(resList) == 0) {
    ret <- list(
      coefficients = NULL,
      design = NULL,
      rdf = NULL,
      df.residual = NULL,
      hatvalues = NULL,
      edf = NULL,
      logLik = NULL,
      Amean = NULL,
      method = "lmer",
      sigma = NULL,
      contrasts = L,
      stdev.unscaled = NULL,
      residuals = NULL
    )
    ret <- new("MArrayLM", ret)
    return(ret)
  }

  coefficients <- t(do.call(cbind, lapply(resList, function(x) x$ret$coefficients)))
  rownames(coefficients) <- names(resList)
  colnames(coefficients) <- colnames(L)

  df.residual <- t(do.call(cbind, lapply(resList, function(x) x$ret$df.residual)))
  colnames(df.residual) <- colnames(L)

  rdf <- c(do.call(cbind, lapply(resList, function(x) x$ret$rdf)))
  names(rdf) <- names(resList)

  pValue <- t(do.call(cbind, lapply(resList, function(x) x$ret$pValue)))
  colnames(pValue) <- colnames(L)

  stdev.unscaled <- t(do.call(cbind, lapply(resList, function(x) x$ret$stdev.unscaled)))
  rownames(stdev.unscaled) <- names(resList)

  residuals <- t(do.call(cbind, lapply(resList, function(x) x$ret$residuals)))

  hatvalues <- t(do.call(cbind, lapply(resList, function(x) x$hatvalues)))

  logLik <- sapply(resList, function(x) x$logLik)

  design <- resList[[1]]$ret$design

  Amean <- sapply(resList, function(x) x$ret$Amean)
  sigma <- sapply(resList, function(x) x$ret$sigma)

  varComp <- lapply(resList, function(x) {
    x <- unlist(x$varComp)
    names(x) <- gsub("\\.\\(Intercept\\)", "", names(x))
    as.data.frame(t(x))
  })
  varComp <- do.call("rbind", varComp)
  rownames(varComp) <- names(resList)

  edf <- sapply(resList, function(x) x$edf)

  ret <- list(
    coefficients = coefficients,
    design = design,
    rdf = c(rdf),
    df.residual = df.residual,
    hatvalues = hatvalues,
    edf = rowSums(hatvalues),
    logLik = logLik,
    Amean = Amean,
    method = "lmer",
    sigma = sigma,
    contrasts = L,
    stdev.unscaled = stdev.unscaled,
    residuals = residuals
  )

  if ("genes" %in% names(exprObj)) {
    ret$genes <- exprObj$genes
  }

  ret <- new("MArrayLM", ret)
  # ret$pValue = pValue

  # set covariance between covariates
  # C = solve(crossprod(ret$design))
  # C = chol2inv(chol(crossprod(ret$design)))
  V <- chol2inv(qr(ret$design)$qr)
  rownames(V) <- colnames(ret$design)
  colnames(V) <- colnames(ret$design)

  if (!univariateContrasts) {
    # do expensive evaluation only when L is defined by the user
    V <- crossprod(chol(V) %*% L)
  }

  # covariance used by limma.
  # Assumes covariance is the same for all genes
  ret$cov.coefficients <- V

  # allows covariance to differ for each gene based on variance components
  ret$cov.coefficients.list <- lapply(resList, function(x) as.matrix(x$vcov))

  ret <- as(ret, "MArrayLM2")
  # add additional information for pinnacle
  # ret = new("MArrayLM2", ret, varComp, sigGStruct)
  attr(ret, "varComp") <- varComp
  # attr(ret, "sigGStruct") = get_SigmaG( fit )$G
  attr(ret, "edf") <- edf

  # compute standard values for t/F/p without running eBayes
  # eBayes can be run afterwards, if wanted
  ret <- .standard_transform(ret)
}



isPositiveDefinite <- function(x, tol = 1e-8) {
  # ensure matrix is square
  stopifnot(nrow(x) == ncol(x))

  # compute eigen values
  ev <- eigen(x, only.values = TRUE)$values

  all(ev > tol)
}
