#' Class MArrayLM2
#'
#' Class \code{MArrayLM2}
#'
#' @name MArrayLM2-class
#' @rdname MArrayLM2-class
#' @exportClass MArrayLM2
setClass(
  "MArrayLM2",
  #  Linear model fit
  representation("MArrayLM") # , varComp="data.frame", sigGStruct='list')
)

setIs("MArrayLM2", "LargeDataObject")
setAs(from = "MArrayLM", to = "MArrayLM2", function(from) {
  structure(from, class = "MArrayLM2")
})


# define S3 version of these functions

#' Residuals for result of dream
#'
#' Residuals for result of dream
#'
#' @param object See \code{?stats::residuals}
#' @param y \code{EList} object used in \code{dream()}
#' @param ... See \code{?stats::residuals}
#' @param type compute either response or pearson residuals
#'
#' @rawNamespace S3method("residuals", MArrayLM2)
#' @export
residuals.MArrayLM2 <- function(object, y, ..., type = c("response", "pearson")) {
  stopifnot(is(object, "MArrayLM2"))

  type <- match.arg(type)

  if (is.null(object$residuals)) {
    stop("Residuals were not computed, must run:\n dream(...,computeResiduals=TRUE)")
  }

  if (!missing(y)) {
    # subset to intersecting features
    featureIds <- intersect(rownames(object$residuals), rownames(y))

    object <- object[featureIds, ]
    object$residuals <- object$residuals[featureIds, , drop = FALSE]
    y <- y[featureIds, , drop = FALSE]

    if (!all.equal(dim(object$residuals), dim(y))) {
      stop("Dimension of object and y must be the same")
    }
    if (!all.equal(rownames(y), rownames(object))) {
      stop("rownames of object and y must be the same")
    }
    if (!all.equal(colnames(y), colnames(object$residuals))) {
      stop("colnames of object and y must be the same")
    }
  }

  residResponse <- object$residuals

  if (type == "response") {
    result <- residResponse
  } else if (type == "pearson") {
    if (missing(y)) {
      stop("Original EList required to fit pearson residuals")
    }
    # from lmer() fit
    # residuals(fitlmer, type="response") * sqrt(w) / sqrt(1-h)
    if (!is.null(y$weights)) {
      w <- y$weights
    } else {
      w <- 1
    }
    # with linear mixed model, hatvalues can be effectively 1
    # so add a small value, sqrt(.Machine$double.eps), to keep this 
    # value positive
    result <- residResponse * sqrt(w) / sqrt(1 - object$hatvalues + sqrt(.Machine$double.eps))
  }
  result
}





# #' @importFrom limma residuals.MArrayLM
# # @export
# residuals.MArrayLM = function( object, ...){

# 	if( is.null(object$residuals) ){
# 		# use residuals computed by limma
# 		res = limma::residuals.MArrayLM( object, ...)
# 	}else{
# 		# use precomputed residuals
# 		res = object$residuals
# 	}
# 	res
# }


# S4 methpds

#' residuals for MArrayLM
#'
#' residuals for MArrayLM
#'
#' @param object MArrayLM object from dream
#' @param y \code{EList} object used in \code{dream()}
#' @param ... other arguments, currently ignored
#' @param type compute either response or pearson residuals
#'
#' @return results of residuals
#' @export
setMethod(
  "residuals", "MArrayLM",
  function(object, y, ..., type = c("response", "pearson")) {
    stopifnot(is(object, "MArrayLM"))

    type <- match.arg(type)

    if (type == "response" & !missing(y)) {
      residResponse <- residuals.MArrayLM(object, y = y, ...)
      return(residResponse)
    }

    if (is.null(object$residuals)) {
      stop("Residuals were not computed, must run:\n dream(...,computeResiduals=TRUE)")
    }

    if (!missing(y)) {
      if (!all.equal(dim(object$residuals), dim(y))) {
        stop("Dimension of object and y must be the same")
      }
      if (!all.equal(rownames(y), rownames(object))) {
        stop("rownames of object and y must be the same")
      }
      if (!all.equal(colnames(y), colnames(object$residuals))) {
        stop("colnames of object and y must be the same")
      }
    }

    residResponse <- object$residuals

    if (type == "response") {
      result <- residResponse
    } else if (type == "pearson") {
      # from glm() fit
      # residuals(fitlm, type="response") * sqrt(fitlm$weights) / sqrt(1-h)
      if (!is.null(y$weights)) {
        w <- y$weights
      } else {
        w <- 1
      }
      result <- residResponse * sqrt(w) / sqrt(1 - object$hatvalues)
    }
    result
  }
)


#' residuals for MArrayLM2
#'
#' residuals for MArrayLM2
#'
#' @param object MArrayLM2 object from dream
#' @param y \code{EList} object used in \code{dream()}
#' @param ... other arguments, currently ignored
#' @param type compute either response or pearson residuals
#'
#' @return results of residuals
#' @export
setMethod(
  "residuals", "MArrayLM2",
  function(object, y, type = c("response", "pearson"), ...) {
    residuals.MArrayLM2(object, y = y, ..., type = type)
  }
)




# Evaluate contrasts for linear mixed model
#
# Evaluate contrasts for linear mixed model
#
# @param fit model fit
# @param L contrast matrix
# @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough approximation for most datasets.
#
# @return
# df, sigma, beta, SE of model
#
# @details If the Kenward-Roger covariance matrix is not positive definite, the Satterthwaite method is used
#
# @export
# @docType methods
# @rdname eval_lmm-method
# @importFrom lmerTest contest
# @importFrom lme4 fixef
# @importFrom stats sigma
# .eval_lmm = function( fit, L, ddf ){

# 	j = 1
# 	# evaluate each contrast
# 	# cons = lmerTest::contest(fit, L, ddf=ddf)
# 	# cons = foreach( j = 1:ncol(L), .combine=rbind) %do% {
# 	# 	lmerTest::contest(fit, L[,j], ddf=ddf)
# 	# }

# 	cons = lapply( seq(ncol(L)), function(j){
# 		contest(fit, L[,j], ddf=ddf)
# 		})
# 	cons = do.call(rbind, cons)

# 	df = as.numeric(cons[,'DenDF'])

# 	if(ddf == "Kenward-Roger"){
# 		# KR
# 		V = pbkrtest::vcovAdj.lmerMod(fit, 0)

# 		# if matrix is not PSD
# 		if( min(diag(as.matrix(V))) < 0){
# 			warning("The adjusted Kenward-Roger covariance matrix is not positive definite.\nUsing Satterthwaite approximation instead")

# 			# Satterthwaite
# 			V = vcov(fit)
# 		}
# 		# df = pbkrtest::get_Lb_ddf(fit, L)
# 	}else{
# 		# Satterthwaite
# 		V = vcov(fit)
# 		# df = as.numeric(contest(fit, L, ddf="Sat")['DenDF'])
# 	}

# 	# sigma = attr(lme4::VarCorr(fit), "sc")

# 	# get contrasts
# 	# beta = as.matrix(sum(L * fixef(fit)), ncol=1)
# 	# colnames(beta) = "logFC"

# 	# beta = foreach( j = 1:ncol(L), .combine=rbind) %do% {
# 	# 	as.matrix(sum(L[,j] * fixef(fit)), ncol=1)
# 	# }
# 	beta = lapply( seq(ncol(L)), function(j){
# 		as.matrix(sum(L[,j] * fixef(fit)), ncol=1)
# 		})
# 	beta = do.call(rbind, beta)

# 	colnames(beta) = "logFC"
# 	rownames(beta) = colnames(L)

# 	# SE = as.matrix(sqrt(sum(L * (V %*% L))), ncol=1)
# 	# colnames(SE) = "logFC"
# 	# SE = foreach( j = 1:ncol(L), .combine=rbind) %do% {
# 	# 	as.matrix(sqrt(sum(L[,j] * (V %*% L[,j]))), ncol=1)
# 	# }
# 	SE = lapply( seq(ncol(L)), function(j){
# 		as.matrix(sqrt(sum(L[,j] * (V %*% L[,j]))), ncol=1)
# 		})
# 	SE = do.call(rbind, SE)

# 	colnames(SE) = "logFC"
# 	rownames(SE) = colnames(L)

# 	# pValue = 2*pt(as.numeric(abs(beta / SE)), df, lower.tail=FALSE)
# 	pValue = as.numeric(cons[,'Pr(>F)'])

# 	list(	cons 	= cons,
# 			df		= df,
# 			sigma	= sigma(fit),
# 			beta	= beta,
# 			SE		= SE,
# 			pValue	= pValue,
# 			vcov 	= V )
# }


#' @importFrom methods is
.checkNA <- function(exprObj) {
  if (is(exprObj, "sparseMatrix") || is(exprObj, "matrix")) {
    countNA <- sum(!is.finite(exprObj))
  } else {
    # is.finite is not defined for data.frames, so convert to matrix first
    countNA <- sum(!is.finite(as.matrix(exprObj)))

    # check if values are NA
    # countNA = sum(!is.finite(exprObj)) # sum(is.nan(exprObj))
  }

  if (countNA > 0) {
    stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
  }

  # check if all genes have variance
  rv <- apply(exprObj, 1, var)
  if (any(rv == 0)) {
    idx <- which(rv == 0)
    stop(paste("Response variable", idx[1], "has a variance of 0"))
  }
}

#' Compute standard post-processing values
#'
#' These values are typically computed by eBayes
#'
#' @param fit result of dream (MArrayLM2)
#' @param sigma vector of standard errors used to compute t-statistic. Can be maximum likelihood estimates, or posterior means
#'
#' @return MArrayLM2 object with values computed
#'
#' @importFrom stats pchisq pf
#' @keywords internal
.standard_transform <- function(fit, sigma = fit$sigma) {
  # If fit$df.prior is not defined, set df.prior to zero
  if (!is.null(fit$df.prior)) {
    fit$df.total <- fit$df.residual + fit$df.prior
  } else {
    fit$df.total <- fit$df.residual
  }

  # t-test
  out <- fit
  out$t <- fit$coefficients / fit$stdev.unscaled / sigma
  out$p.value <- 2 * pt(abs(out$t), df = fit$df.total, lower.tail = FALSE)
  # out$p.value.loge <- log(2) + pt(abs(out$t), df=fit$df.total, lower.tail=FALSE, log.p=TRUE )

  # F-test
  if (!is.null(out$design) && is.fullrank(out$design)) {
    # only evaluate F-stat on real coefficients, not contrasts
    realcoef <- colnames(out)[colnames(out) %in% colnames(out$design)]
    realcoef <- realcoef[realcoef != "(Intercept)"]

    if (is.null(realcoef) || (length(realcoef) == 0)) {
      # this happends when only the intercept term is included
      warning("No testable fixed effects were included in the model.\n  Running topTable() will fail.")
    } else {
      df <- rowMeans(out[, realcoef]$df.total)

      F.stat <- classifyTestsF(out[, realcoef], df = df, fstat.only = TRUE)
      out$F <- as.vector(F.stat)
      df1 <- attr(F.stat, "df1")
      df2 <- attr(F.stat, "df2")
      if (df2[1] > 1e6) { # Work around bug in R 2.1
        out$F.p.value <- pchisq(df1 * out$F, df1, lower.tail = FALSE)
      } else {
        out$F.p.value <- pf(out$F, df1, df2, lower.tail = FALSE)
      }
    }
  }

  # if fit$df.prior does not exist, then remove the df.total term
  if (is.null(fit$df.prior)) {
    out$df.total <- NULL
  }

  out
}




#' Subseting for MArrayLM2
#'
#' Enable subsetting on MArrayLM2 object. Same as for MArrayLM, but apply column subsetting to df.residual and cov.coefficients.list
#'
#' @param object MArrayLM2
#' @param i row
#' @param j col
#'
#' @name [.MArrayLM2
#' @return subset
#' @rawNamespace S3method("[", MArrayLM2)
#' @importFrom stats p.adjust
#' @rdname subset.MArrayLM2-method
#' @aliases subset.MArrayLM2,MArrayLM2-method
#' @keywords internal
#' @export
assign(
  "[.MArrayLM2",
  function(object, i, j) {
    if (nargs() != 3) {
      stop("Two subscripts required", call. = FALSE)
    }

    if( !missing(j) ){
      if( is.logical(j) ) j <- which(j)
    }

    if( !missing(i) ){
      if( is.logical(i) ) i <- which(i)
    }

    # apply standard MArrayLM subsetting
    obj <- as(object, "MArrayLM")

    if (!missing(j)) {
      obj <- obj[, j]
    }
    if (!missing(i)) {
      obj <- obj[i, ]
    }

    # custom code to deal with df.total, df.residual and rdf
    if (is.null(ncol(object$df.total))) {
      if (!missing(i)) {
        obj$df.total <- object$df.total[i]
      } else {
        obj$df.total <- object$df.total
      }
    } else {
      tmp <- object$df.total
      if (!missing(i)) {
        tmp <- object$df.total[i, , drop = FALSE]
      }
      if (!missing(j)) {
        tmp <- tmp[, j, drop = FALSE]
      }
      obj$df.total <- tmp
    }

    if (is.null(ncol(object$df.residual))) {
      if (!missing(i)) {
        obj$df.residual <- object$df.residual[i]
      } else {
        obj$df.residual <- object$df.residual
      }
    } else {
      tmp <- object$df.residual
      if (!missing(i)) {
        tmp <- object$df.residual[i, , drop = FALSE]
      }
      if (!missing(j)) {
        tmp <- tmp[, j, drop = FALSE]
      }
      obj$df.residual <- tmp
    }

    if (!is.null(object$rdf)) {
      if (!missing(i)) {
        names(object$rdf) <- rownames(object)
        obj$rdf <- object$rdf[i]
        names(obj$rdf) <- c()
      } else {
        obj$rdf <- object$rdf
      }
    }

    if (!is.null(object$edf)) {
      if (!missing(i)) {
        names(object$edf) <- rownames(object)
        obj$edf <- object$edf[i]
        # names(obj$edf) <- c()
      } else {
        obj$edf <- object$edf
      }
    }

    if (!is.null(object$logLik)) {
      if (!missing(i)) {
        names(object$logLik) <- rownames(object)
        obj$logLik <- object$logLik[i]
        # names(obj$logLik) <- c()
      } else {
        obj$logLik <- object$logLik
      }
    }

    # obj$pValue = object$pValue[i,j]
    obj$s2.prior <- object$s2.prior
    obj$df.prior <- object$df.prior

    obj <- as(obj, "MArrayLM2")

    #  copy gene-specific covariance, if it exists
    if (!is.null(object$cov.coefficients.list)) {
      if (!missing(i)) {
        if (is.numeric(i)) {
          # extract by index
          obj$cov.coefficients.list <- object$cov.coefficients.list[i]
        } else {
          # extract by matching feature name
          idx <- match(i, rownames(object))
          obj$cov.coefficients.list <- object$cov.coefficients.list[idx]
        }
      } else {
        obj$cov.coefficients.list <- object$cov.coefficients.list
      }
      # name cov.coefficients.list using names of the whole object
      names(obj$cov.coefficients.list) <- rownames(obj)
    }

    if (is.null(obj$df.total)) {
      obj$df.total <- rowMeans(obj$df.residual)
    }

    # subset residuals and hatvalues
    obj$residuals = obj$residuals[rownames(obj),,drop=FALSE]
    obj$hatvalues = obj$hatvalues[rownames(obj),,drop=FALSE]

    # the F-statistic and p-value are evaluated when subsetting is applied
    # so need to apply df2 here
    # If columns have been subsetted, need to re-generate F
    if (!is.null(obj[["F"]]) && !missing(j)) {
      F.stat <- classifyTestsF(obj, df = obj$df.total, fstat.only = TRUE)
      obj$F <- as.vector(F.stat)
      df1 <- attr(F.stat, "df1")
      df2 <- attr(F.stat, "df2")
      if (df2[1] > 1e6) {
        obj$F.p.value <- pchisq(df1 * obj$F, df1, lower.tail = FALSE)
      } else {
        obj$F.p.value <- pf(obj$F, df1, df2, lower.tail = FALSE)
      }
    }
    obj
  }
)


setGeneric("eBayes", function(
    fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4),
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1), legacy = NULL) {
  eBayes(fit, proportion, stdev.coef.lim, trend, robust, winsor.tail.p, legacy)
})


#' eBayes for MArrayLM2
#'
#' eBayes for result of linear mixed model for with \code{dream()} using residual degrees of freedom approximated with \code{rdf.merMod()}
#'
#' @param fit fit
#' @param proportion proportion
#' @param stdev.coef.lim stdev.coef.lim
#' @param trend trend
#' @param robust robust
#' @param winsor.tail.p winsor.tail.p
#' @param legacy legacy
#'
#' @return results of eBayes using approximated residual degrees of freedom
#'
#' @export
#' @rdname eBayes-method
#' @aliases eBayes,MArrayLM2-method
#' @importFrom limma eBayes
#' @seealso \code{dream()}, \code{rdf.merMod()}, \code{limma::eBayes()}
setMethod(
  "eBayes", "MArrayLM2",
  function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4),
           trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1), legacy = NULL) {
    # limma::eBayes() uses df.residual as the residual degrees of freedom,
    # 	while dream() uses rdf.
    # For linear models these values are always equal,
    # but for linear mixed models, rdf must be computed separately
    # save df for test statistics
    df.test <- fit$df.residual
    fit$df.residual <- fit$rdf

    # Use limma::eBayes() but with new df.residual values
    fit_eb <- limma::eBayes(
      fit = fit,
      proportion = proportion,
      stdev.coef.lim = stdev.coef.lim,
      trend = trend,
      robust = robust,
      winsor.tail.p = winsor.tail.p,
      legacy = legacy
    )

    # re-set to the df.residual of the test statistics
    fit_eb$df.residual <- df.test
    fit_eb$rdf <- fit$rdf

    # Calculate p-values using estimated degrees of freedom
    # and posterior variance estimates
    .standard_transform(fit_eb, sigma = sqrt(fit_eb$s2.post))
  }
)




#' Compare p-values from two analyses
#'
#' Plot -log10 p-values from two analyses and color based on donor component from variancePartition analysis
#'
#' @param p1 p-value from first analysis
#' @param p2 p-value from second analysis
#' @param vpDonor donor component for each gene from variancePartition analysis
#' @param dupcorvalue scalar donor component from duplicateCorrelation
#' @param fraction fraction of highest/lowest values to use for best fit lines
#' @param xlabel for x-axis
#' @param ylabel label for y-axis
#'
#' @return ggplot2 plot
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
#' # Perform very simple analysis for demonstration
#'
#' # Analysis 1
#' form <- ~Batch
#' fit <- dream(geneExpr, form, info)
#' fit <- eBayes(fit)
#' res <- topTable(fit, number = Inf, coef = "Batch3")
#'
#' # Analysis 2
#' form <- ~ Batch + (1 | Tissue)
#' fit2 <- dream(geneExpr, form, info)
# fitEB2 = eBayes( fit2 )
#' res2 <- topTable(fit2, number = Inf, coef = "Batch3")
#'
#' # Compare p-values
#' plotCompareP(res$P.Value, res2$P.Value, runif(nrow(res)), .3)
#'
#' @export
# @docType methods
#' @rdname plotCompareP-method
plotCompareP <- function(p1, p2, vpDonor, dupcorvalue, fraction = .2, xlabel = bquote(duplicateCorrelation ~ (-log[10] ~ p)), ylabel = bquote(dream ~ (-log[10] ~ p))) {
  if (length(unique(c(length(p1), length(p2), length(vpDonor)))) != 1) {
    stop("p1, p2 and vpDonor must have the same number of entries")
  }
  if (length(dupcorvalue) != 1) {
    stop("dupcorvalue must be a scalar")
  }

  df2 <- data.frame(p1 = -log10(p1), p2 = -log10(p2), vpDonor = vpDonor, delta = vpDonor - dupcorvalue)

  N <- nrow(df2)

  c1 <- sort(df2$delta)[N * fraction]
  c2 <- sort(df2$delta)[length(df2$delta) - N * fraction]

  l1 <- lm(p2 ~ p1, df2[df2$delta >= c2, ])
  l2 <- lm(p2 ~ p1, df2[df2$delta <= c1, ])

  df_line <- data.frame(rbind(coef(l1), coef(l2)))
  colnames(df_line) <- c("a", "b")
  df_line$type <- c("darkred", "navy")

  lim <- c(0, max(max(df2$p1), max(df2$p2)))

  # xlab("duplicateCorrelation (-log10 p)") + ylab("dream (-log10 p)")
  ggplot(df2, aes(p1, p2, color = vpDonor)) +
    geom_abline() +
    geom_point(size = 2) +
    theme_bw(17) +
    theme(aspect.ratio = 1, plot.title = element_text(hjust = 0.5)) +
    xlim(lim) +
    ylim(lim) +
    xlab(xlabel) +
    ylab(ylabel) +
    geom_abline(intercept = df_line$a, slope = df_line$b, color = df_line$type, linetype = 2) +
    scale_color_gradientn(
      name = "Donor", colours = c("blue", "green", "red"),
      values = rescale(c(0, dupcorvalue, 1)),
      guide = "colorbar", limits = c(0, 1)
    )
}

#' Multiple Testing Genewise Across Contrasts
#'
#'  For each gene, classify a series of related t-statistics as up, down or not significant.
#'
#' @param object numeric matrix of t-statistics or an 'MArrayLM2' object from which the t-statistics may be extracted.
#' @param ... additional arguments
#'
#' @details Works like limma::classifyTestsF, except object can have a list of covariance matrices object$cov.coefficients.list, instead of just one in object$cov.coefficients
#' @seealso \code{limma::classifyTestsF}
# @export
setGeneric("classifyTestsF",
  signature = "object",
  function(object, ...) {
    standardGeneric("classifyTestsF")
  }
)



#' Multiple Testing Genewise Across Contrasts
#'
#'  For each gene, classify a series of related t-statistics as up, down or not significant.
#'
#' @param object numeric matrix of t-statistics or an 'MArrayLM2' object from which the t-statistics may be extracted.
#' @param cor.matrix covariance matrix of each row of t-statistics.  Defaults to the identity matrix.
#' @param df numeric vector giving the degrees of freedom for the t-statistics.  May have length 1 or length equal to the number of rows of tstat.
#' @param p.value numeric value between 0 and 1 giving the desired size of the test
#' @param fstat.only logical, if 'TRUE' then return the overall F-statistic as for 'FStat' instead of classifying the test results
#'
#' @details Works like limma::classifyTestsF, except object can have a list of covariance matrices object$cov.coefficients.list, instead of just one in object$cov.coefficients
#' @seealso \code{limma::classifyTestsF}
#' @importFrom stats qf
# @export
setMethod(
  "classifyTestsF", "MArrayLM2",
  function(object, cor.matrix = NULL, df = Inf, p.value = 0.01, fstat.only = FALSE) {
    # 	Use F-tests to classify vectors of t-test statistics into outcomes
    # 	Gordon Smyth
    # 	20 Mar 2003.  Last revised 6 June 2009.

    # 	Method intended for MArrayLM objects but accept unclassed lists as well
    if (is.list(object)) {
      if (is.null(object$t)) stop("tstat cannot be extracted from object")
      computeCorrMat <- ifelse(is.null(cor.matrix), TRUE, FALSE)
      if (missing(df) && !is.null(object$df.prior) && !is.null(object$df.residual)) {
        df <- object$df.prior + object$df.residual
      }
      tstat <- as.matrix(object$t)
    } else {
      tstat <- as.matrix(object)
    }
    ngenes <- nrow(tstat)
    ntests <- ncol(tstat)
    if (ntests == 1) {
      if (fstat.only) {
        fstat <- tstat^2
        attr(fstat, "df1") <- 1
        attr(fstat, "df2") <- df
        return(fstat)
      } else {
        p <- 2 * pt(abs(tstat), df, lower.tail = FALSE)
        return(new("TestResults", sign(tstat) * (p < p.value)))
      }
    }

    # F test of multiple coefficients
    #-------------------------------

    if (ngenes != length(object$cov.coefficients.list)) {
      msg <- paste0("Number of genes does not equal number of elements in cov.coefficients.list\n", ngenes, " != ", length(object$cov.coefficients.list))
      stop(msg)
    }

    fstat <- rep(NA, ngenes)
    names(fstat) <- rownames(tstat)

    result <- matrix(0, ngenes, ntests, dimnames = dimnames(tstat))

    for (i in seq_len(ngenes)) {
      if (computeCorrMat) {
        if (is.null(object$cov.coefficients.list)) {
          C <- object$cov.coefficients
        } else {
          C <- object$cov.coefficients.list[[i]]
        }
        # subset based on coefficient names
        C <- C[colnames(object), colnames(object)]
        cor.matrix <- cov2cor(C)
      } else {
        cor.matrix <- cov2cor(cor.matrix)
      }

      # cor.matrix is estimated correlation matrix of the coefficients
      # and also the estimated covariance matrix of the t-statistics
      if (is.null(cor.matrix)) {
        r <- ntests
        Q <- diag(r) / sqrt(r)
      } else {
        E <- eigen(cor.matrix, symmetric = TRUE)
        r <- sum(E$values / E$values[1] > 1e-8)
        Q <- limma:::.matvec(E$vectors[, 1:r], 1 / sqrt(E$values[1:r])) / sqrt(r)
      }

      # Return overall moderated F-statistic only
      if (fstat.only) {
        fstat[i] <- drop((tstat[i, , drop = FALSE] %*% Q)^2 %*% array(1, c(r, 1)))
      }
      if (i == 1) {
        attr(fstat, "df1") <- r
        attr(fstat, "df2") <- df[i]
      }

      # Return TestResults matrix
      qF <- qf(p.value, r, df[i], lower.tail = FALSE)
      if (length(qF) == 1) qF <- rep(qF, ngenes)
      x <- tstat[i, ]
      if (any(is.na(x))) {
        result[i, ] <- NA
      } else if (crossprod(crossprod(Q, x)) > qF[i]) {
        ord <- order(abs(x), decreasing = TRUE)
        result[i, ord[1]] <- sign(x[ord[1]])
        for (j in 2:ntests) {
          bigger <- ord[1:(j - 1)]
          x[bigger] <- sign(x[bigger]) * abs(x[ord[j]])
          if (crossprod(crossprod(Q, x)) > qF[i]) {
            result[i, ord[j]] <- sign(x[ord[j]])
          } else {
            break
          }
        }
      }
    }

    if (fstat.only) {
      return(fstat)
    }

    new("TestResults", result)
  }
)
