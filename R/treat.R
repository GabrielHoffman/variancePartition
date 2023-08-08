#' Test if coefficient is different from a specified value
#'
#' Test if coefficient is different from a specified value
#'
#' @param fit fit
#' @param lfc a minimum log2-fold-change below which changes not considered scientifically meaningful
#' @param coef which coefficient to test
#' @param number number of genes to return
#' @param sort.by column to sort by
#'
#' @return results of getTreat
#'
#' @examples
#'
#' data(varPartData)
#'
#' form <- ~ Age + Batch + (1 | Individual) + (1 | Tissue)
#'
#' fit <- dream(geneExpr, form, info)
#' fit <- eBayes(fit)
#'
#' coef <- "Age"
#'
#' # Evaluate treat()/topTreat() in a way that works seamlessly for dream()
#' getTreat(fit, lfc = log2(1.03), coef, sort.by = "none", number = 3)
#'
#' @export
#' @rdname getTreat-method
#' @aliases getTreat,MArrayLM-method
setGeneric("getTreat", function(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p") {
  standardGeneric("getTreat")
})


# Adapted from limma::treat, but returns results table after eBayes has been run
#' @export
#' @rdname getTreat-method
#' @aliases getTreat,MArrayLM-method
setMethod(
  "getTreat", "MArrayLM",
  function(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p") {
    se <- fit$coefficients[, coef] / fit$t[, coef]

    tstat.right <- (abs(fit$coefficients[, coef]) - lfc) / se
    tstat.left <- (abs(fit$coefficients[, coef]) + lfc) / se

    t <- rep(0, length(tstat.right))

    p.value <- pt(tstat.right, df = fit$df.total, lower.tail = FALSE) + pt(tstat.left, df = fit$df.total, lower.tail = FALSE)

    tstat.right <- pmax(tstat.right, 0)
    fc.up <- (fit$coefficients[, coef] > lfc)
    fc.down <- (fit$coefficients[, coef] < -lfc)

    t[fc.up] <- tstat.right[fc.up]
    t[fc.down] <- -tstat.right[fc.down]

    fit.internal <- fit
    fit.internal$t[, coef] <- t
    fit.internal$p.value[, coef] <- p.value

    # topTable(treat(fit, lfc=log2(1.2)), coef=coef, sort.by=sort.by)

    topTable(fit.internal, coef = coef, number = number, sort.by = sort.by)
  }
)


# Adapted from limma::treat, but returns results table after eBayes has been run
#' @export
#' @rdname getTreat-method
#' @aliases getTreat,MArrayLM2-method
setMethod(
  "getTreat", "MArrayLM2",
  function(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p") {
    se <- fit$coefficients[, coef] / fit$t[, coef]

    tstat.right <- (abs(fit$coefficients[, coef]) - lfc) / se
    tstat.left <- (abs(fit$coefficients[, coef]) + lfc) / se

    t <- rep(0, length(tstat.right))

    p.value <- pt(tstat.right, df = fit$df.total[, coef], lower.tail = FALSE) + pt(tstat.left, df = fit$df.total[, coef], lower.tail = FALSE)

    tstat.right <- pmax(tstat.right, 0)
    fc.up <- (fit$coefficients[, coef] > lfc)
    fc.down <- (fit$coefficients[, coef] < -lfc)

    t[fc.up] <- tstat.right[fc.up]
    t[fc.down] <- -tstat.right[fc.down]

    fit.internal <- fit
    fit.internal$t[, coef] <- t
    fit.internal$p.value[, coef] <- p.value

    # topTable(treat(fit, lfc=log2(1.2)), coef=coef, sort.by=sort.by)

    topTable(fit.internal, coef = coef, number = number, sort.by = sort.by)
  }
)
