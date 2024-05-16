#' Table of Top Genes from Linear Model Fit
#'
#' topTable generic
#'
#' @param fit fit
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @return results of toptable
#' @export
#' @importFrom stats qnorm
#' @import limma
#' @rdname toptable-method
#' @aliases toptable,MArrayLM-method
setGeneric("topTable", function(fit, coef = NULL, number = 10, genelist = fit$genes, adjust.method = "BH", sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0, confint = FALSE) {
  standardGeneric("topTable")
})



#' Table of Top Genes from Linear Model Fit
#'
#' topTable generic MArrayLM
#'
#' @param fit fit
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @return results of toptable
#' @export
#' @importFrom stats qnorm
#' @import limma
#' @rdname toptable-method
#' @aliases toptable,MArrayLM-method
setMethod(
  "topTable", "MArrayLM",
  function(fit, coef = NULL, number = 10, genelist = fit$genes,
           adjust.method = "BH", sort.by = "p", resort.by = NULL, p.value = 1,
           lfc = 0, confint = FALSE) {
   
    if( length(coef) > 1 & ! sort.by %in% c("F", "none")){
      sort.by = "F"
    }

    limma::topTable(fit, coef = coef, number = number, genelist = genelist, adjust.method = adjust.method, sort.by = sort.by, resort.by = resort.by, p.value = p.value, lfc = lfc, confint = confint)
  }
)


#' Table of Top Genes from Linear Model Fit
#'
#' topTable generic MArrayLM2
#'
#' @param fit fit
#' @param coef coef
#' @param number number
#' @param genelist genelist
#' @param adjust.method adjust.method
#' @param sort.by sort.by
#' @param resort.by resort.by
#' @param p.value p.value
#' @param lfc lfc
#' @param confint confint
#'
#' @return results of toptable
#' @export
#' @importFrom stats qnorm
#' @import limma
#' @rdname toptable-method
#' @aliases toptable,MArrayLM2-method
setMethod(
  "topTable", "MArrayLM2",
  function(fit, coef = NULL, number = 10, genelist = fit$genes,
           adjust.method = "BH", sort.by = "p", resort.by = NULL, p.value = 1,
           lfc = 0, confint = FALSE) {
    if (!is(fit, "MArrayLM2")) {
      stop("fit must be an MArrayLM2 object")
    }
    if (is.null(fit$t) && is.null(fit$F)) {
      stop("Need to run eBayes or treat first")
    }
    if (is.null(fit$coefficients)) {
      stop("coefficients not found in fit object")
    }
    if (confint && is.null(fit$stdev.unscaled)) {
      stop("stdev.unscaled not found in fit object")
    }
    if (is.null(coef)) {
      if (is.null(fit$treat.lfc)) {
        coef <- 1:ncol(fit)
        cn <- colnames(fit)
        if (!is.null(cn)) {
          i <- which(cn == "(Intercept)")
          if (length(i)) {
            coef <- coef[-i]
            message("Removing intercept from test coefficients")
          }
        }
      } else {
        coef <- ncol(fit)
      }
    }
    if (length(coef) > 1) {
      if (!is.null(fit$treat.lfc)) {
        stop("Treat p-values can only be displayed for single coefficients")
      }
      coef <- unique(coef)
      n.coef <- length(coef)
      if (length(fit$coef[1, coef]) < ncol(fit)) {
        # F-test is performed in subseting code by classifyTestsF()
        fit <- fit[, coef]
        coef <- NULL # since subsetting is performed here, don't need to do it later
      }
      if( ! sort.by %in% c("F", "none")){
        sort.by = "F"
      }
      # tab = topTableF(fit, number = number, genelist = genelist,
      #     adjust.method = adjust.method, sort.by = sort.by,
      #     p.value = p.value, lfc = lfc)

      tab <- limma::topTable(fit,
        coef = coef,
        number = number, genelist = genelist,
        adjust.method = adjust.method, sort.by = sort.by,
        p.value = p.value, lfc = lfc
      )

      if (nrow(tab) > 0) {
        # GEH September 5, 2019
        # convert p-values to F-statistics
        # this corresponds to constant degrees of freedom equals Inf
        tab$P.Value <- pmax(tab$P.Value, 1e-300)
        tab$F.std <- qf(tab$P.Value, df1 = n.coef, df2 = Inf, lower.tail = FALSE)
      }

      return(tab)
    }
    fit <- unclass(fit)
    ebcols <- c("t", "p.value", "lods")
    if (confint) {
      ebcols <- c("s2.post", "df.total", ebcols)
    }

    tab <- .topTableT(
      fit = fit[c("coefficients", "stdev.unscaled")],
      coef = coef, number = number, genelist = genelist, A = fit$Amean,
      eb = fit[ebcols], adjust.method = adjust.method, sort.by = sort.by,
      resort.by = resort.by, p.value = p.value, lfc = lfc,
      confint = confint
    )

    if (nrow(tab) > 0) {
      # GEH September 5, 2019
      # convert p-values to z-scores
      # this corresponds to constant degrees of freedom equals Inf
      tab$P.Value <- pmax(tab$P.Value, 1e-300)
      tab$z.std <- sign(tab$t) * qnorm(tab$P.Value / 2, lower.tail = FALSE)
    }

    tab
  }
)


.topTableT <- function(fit, coef = 1, number = 10, genelist = NULL, A = NULL, eb = NULL, adjust.method = "BH", sort.by = "B", resort.by = NULL, p.value = 1, lfc = 0, confint = FALSE, ...)
                       # 	Summary table of top genes for a single coefficient
                       # 	Gordon Smyth
# 	21 Nov 2002. Forked from toptable() 1 Feb 2018. Last revised 1 Feb 2018.
{
  # 	Check fit
  fit$coefficients <- as.matrix(fit$coefficients)
  rn <- rownames(fit$coefficients)

  # 	Check coef is length 1
  if (length(coef) > 1) {
    coef <- coef[1]
    warning("Treat is for single coefficients: only first value of coef being used")
  }

  # 	Ensure genelist is a data.frame
  if (!is.null(genelist) && is.null(dim(genelist))) genelist <- data.frame(ID = genelist, stringsAsFactors = FALSE)

  # 	Check rownames
  if (is.null(rn)) {
    rn <- 1:nrow(fit$coefficients)
  } else if (anyDuplicated(rn)) {
    if (is.null(genelist)) {
      genelist <- data.frame(ID = rn, stringsAsFactors = FALSE)
    } else if ("ID" %in% names(genelist)) {
      genelist$ID0 <- rn
    } else {
      genelist$ID <- rn
    }
    rn <- 1:nrow(fit$coefficients)
  }

  # 	Check sort.by
  sort.by <- match.arg(sort.by, c("logFC", "M", "A", "Amean", "AveExpr", "P", "p", "T", "t", "B", "none"))
  if (sort.by == "M") sort.by <- "logFC"
  if (sort.by == "A" || sort.by == "Amean") sort.by <- "AveExpr"
  if (sort.by == "T") sort.by <- "t"
  if (sort.by == "p") sort.by <- "P"

  # 	Check resort.by
  if (!is.null(resort.by)) {
    resort.by <- match.arg(resort.by, c("logFC", "M", "A", "Amean", "AveExpr", "P", "p", "T", "t", "B"))
    if (resort.by == "M") resort.by <- "logFC"
    if (resort.by == "A" || resort.by == "Amean") resort.by <- "AveExpr"
    if (resort.by == "p") resort.by <- "P"
    if (resort.by == "T") resort.by <- "t"
  }

  # 	Check A
  if (is.null(A)) {
    if (sort.by == "A") stop("Cannot sort by A-values as these have not been given")
  } else {
    if (NCOL(A) > 1) A <- rowMeans(A, na.rm = TRUE)
  }

  # 	Check for lods component
  if (is.null(eb$lods)) {
    if (sort.by == "B") stop("Trying to sort.by B, but B-statistic (lods) not found in MArrayLM object", call. = FALSE)
    if (!is.null(resort.by)) if (resort.by == "B") stop("Trying to resort.by B, but B-statistic (lods) not found in MArrayLM object", call. = FALSE)
    include.B <- FALSE
  } else {
    include.B <- TRUE
  }

  # 	Extract statistics for table
  M <- fit$coefficients[, coef]
  tstat <- as.matrix(eb$t)[, coef]
  P.Value <- as.matrix(eb$p.value)[, coef]
  if (include.B) B <- as.matrix(eb$lods)[, coef]

  # 	Apply multiple testing adjustment
  adj.P.Value <- p.adjust(P.Value, method = adjust.method)

  # 	Thin out fit by p.value and lfc thresholds
  if (p.value < 1 | lfc > 0) {
    sig <- (adj.P.Value <= p.value) & (abs(M) >= lfc)
    if (any(is.na(sig))) sig[is.na(sig)] <- FALSE
    if (all(!sig)) {
      return(data.frame())
    }
    genelist <- genelist[sig, , drop = FALSE]
    M <- M[sig]
    A <- A[sig]
    tstat <- tstat[sig]
    P.Value <- P.Value[sig]
    adj.P.Value <- adj.P.Value[sig]
    if (include.B) B <- B[sig]
    rn <- rn[sig]
  }

  # 	Are enough rows left?
  if (length(M) < number) number <- length(M)
  if (number < 1) {
    return(data.frame())
  }

  # 	Select top rows
  ord <- switch(sort.by,
    logFC = order(abs(M), decreasing = TRUE),
    AveExpr = order(A, decreasing = TRUE),
    P = order(P.Value, decreasing = FALSE),
    t = order(abs(tstat), decreasing = TRUE),
    B = order(B, decreasing = TRUE),
    none = 1:length(M)
  )
  top <- ord[1:number]

  # 	Assemble output data.frame
  if (is.null(genelist)) {
    tab <- data.frame(logFC = M[top])
  } else {
    tab <- data.frame(genelist[top, , drop = FALSE], logFC = M[top], stringsAsFactors = FALSE)
  }
  if (confint) {
    if (is.numeric(confint)) {
      alpha <- (1 + confint[1]) / 2
    } else {
      alpha <- 0.975
    }

    # t = beta / s
    # so s = beta/t
    se <- fit$coefficients[top, coef] / eb$t[top, coef]
    margin.error <- se * qt(alpha, df = eb$df.total[top])

    # margin.error <- sqrt(eb$s2.post[top])*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])
    # margin.error <- eb$sigma[top]*fit$stdev.unscaled[top,coef]*qt(alpha,df=eb$df.total[top])

    tab$CI.L <- M[top] - margin.error
    tab$CI.R <- M[top] + margin.error
  }
  if (!is.null(A)) tab$AveExpr <- A[top]
  tab <- data.frame(tab, t = tstat[top], P.Value = P.Value[top], adj.P.Val = adj.P.Value[top])
  if (include.B) tab$B <- B[top]
  rownames(tab) <- rn[top]

  # 	Resort table
  if (!is.null(resort.by)) {
    ord <- switch(resort.by,
      logFC = order(tab$logFC, decreasing = TRUE),
      AveExpr = order(tab$AveExpr, decreasing = TRUE),
      P = order(tab$P.Value, decreasing = FALSE),
      t = order(tab$t, decreasing = TRUE),
      B = order(tab$B, decreasing = TRUE)
    )
    tab <- tab[ord, ]
  }

  tab
}
