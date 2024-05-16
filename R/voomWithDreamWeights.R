#' Transform RNA-Seq Data Ready for Linear Mixed Modelling with \code{dream()}
#'
#' Transform count data to log2-counts per million (logCPM), estimate the mean-variance relationship and use this to compute appropriate observation-level weights. The data are then ready for linear mixed modelling with \code{dream()}.   This method is the same as \code{limma::voom()}, except that it allows random effects in the formula
#'
#' @param counts a numeric \code{matrix} containing raw counts, or an \code{ExpressionSet} containing raw counts, or a \code{DGEList} object. Counts must be non-negative and NAs are not permitted.
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work, and \code{lmFit()} followed by contrasts.fit() are run.
#' @param data \code{data.frame} with columns corresponding to formula
#' @param lib.size numeric vector containing total library sizes for each sample.  Defaults to the normalized (effective) library sizes in \code{counts} if \code{counts} is a \code{DGEList} or to the columnwise count totals if \code{counts} is a matrix.
#' @param normalize.method the microarray-style normalization method to be applied to the logCPM values (if any).  Choices are as for the \code{method} argument of \code{normalizeBetweenArrays} when the data is single-channel.  Any normalization factors found in \code{counts} will still be used even if \code{normalize.method="none"}.
#' @param span width of the lowess smoothing window as a proportion. Setting \code{span="auto"} uses \code{fANCOVA::loess.as()} to estimate the tuning parameter from the data
#' @param weights Can be a numeric matrix of individual weights of same dimensions as the \code{counts}, or a numeric vector of sample weights with length equal to \code{ncol(counts)}
#' @param prior.count average count to be added to each observation to avoid taking log of zero. The count applied to each sample is normalized by library size so given equal log CPM for a gene with zero counts across multiple samples
#' @param prior.count.for.weights count added to regularize weights
#' @param plot logical, should a plot of the mean-variance trend be displayed?
#' @param save.plot logical, should the coordinates and line of the plot be saved in the output?
#' @param rescaleWeightsAfter default = FALSE, should the output weights be scaled by the input weights
#' @param scaledByLib if \code{TRUE}, scale pseudocount by \code{lib.size}.  Else to standard constant pseudocount addition 
#' @param priorWeightsAsCounts if \code{weights} is \code{NULL}, set weights to be equal to counts, following delta method for log2 CPM
#' @param BPPARAM parameters for parallel evaluation
#' @param      ... other arguments are passed to \code{lmer}.
#'
#' @return
#' An \code{EList} object just like the result of \code{limma::voom()}
#'
#' @details Adapted from \code{voom()} in \code{limma} v3.40.2
#' @seealso \code{limma::voom()}
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
#' res <- dream(vobj, form, metadata)
#' res <- eBayes(res)
#'
#' # extract results
#' topTable(res, coef = "Disease1", number = 3)
#'
#' @importFrom lme4 VarCorr
#' @importFrom stats approxfun predict as.formula
#' @importFrom limma asMatrixWeights
#' @importFrom stats sigma
#' @importFrom Matrix t
#' @importFrom matrixStats colSums2 rowSums2
#' @importFrom fANCOVA loess.as
#' @export
voomWithDreamWeights <- function(counts, formula, data, lib.size = NULL, normalize.method = "none", span = 0.5, weights = NULL, prior.count = 0.5, prior.count.for.weights = prior.count, plot = FALSE, save.plot = TRUE, rescaleWeightsAfter = FALSE, scaledByLib = FALSE, priorWeightsAsCounts = FALSE, BPPARAM = SerialParam(), ...) {

  objFlt <- filterInputData(counts, formula, data, weights, useWeights = FALSE, isCounts = TRUE)

  counts <- objFlt$exprObj
  formula <- objFlt$formula
  data <- objFlt$data
  weights <- objFlt$weights

  out <- list()

  design <- NULL

  # 	Check counts
  if (is(counts, "DGEList")) {
    out$genes <- counts$genes
    out$targets <- counts$samples
    # if(is.null(design) && diff(range(as.numeric(counts$sample$group)))>0) design <- model.matrix(~group,data=counts$samples)
    if (is.null(lib.size)) lib.size <- counts$samples$lib.size * counts$samples$norm.factors
    counts <- counts$counts
  } else {
    isExpressionSet <- suppressPackageStartupMessages(is(counts, "ExpressionSet"))
    if (isExpressionSet) {
      if (length(Biobase::fData(counts))) out$genes <- Biobase::fData(counts)
      if (length(Biobase::pData(counts))) out$targets <- Biobase::pData(counts)
      counts <- Biobase::exprs(counts)
    } else if (is.null(dim(counts))) {
      stop("counts is type '", class(counts), "' and can't be converted to matrix unambiguously")
    } else {
      counts <- as.matrix(counts)
    }
  }

  n <- nrow(counts)
  if (n < 2L) stop("Need at least two genes to fit a mean-variance trend")

  # Check lib.size
  if (is.null(lib.size)){
    lib.size <- colSums2(counts)
  }

  # Augment observed counts with prior counts
  # scaled so no variance is introduced across samples
  countsAug = augmentPriorCount(counts, lib.size, prior.count, scaledByLib)

  if( priorWeightsAsCounts && (is.null(weights) | all(weights == 1))){
    weights = countsAug + prior.count.for.weights
  }

  # 	Fit linear model to log2-counts-per-million
  y <- t(log2(t(countsAug) / (lib.size + 1) * 1e6))
  y <- normalizeBetweenArrays(y, method = normalize.method)

  if (is.null(weights)) {
    # set sample-level weights to be equal
    weights <- rep(1, ncol(y))
  }

  # if weights is not already a matrix, convert it to a matrix
  if( is.matrix(weights) ){
    stopifnot( dim(weights) == dim(y))
    weightsMatrix = weights
  }else{

    if (length(weights) != ncol(y)) {
      stop("length of weights must equal ncol(counts)")
    }

    # convert sample-level weights array to matrix
    weightsMatrix <- asMatrixWeights(weights, dim(y))

    # Sept 27, 2023
    # Fixes error reported here
    # https://support.bioconductor.org/p/9154670/
    attr(weightsMatrix,"arrayweights") <- NULL
  }
  rownames(weightsMatrix) <- rownames(y)
  
  # rescale input weights to have mean 1 
  # before scaling output weights
  weightsMatrix <- weightsMatrix / rowMeans(weightsMatrix)  

  # put weights into EList
  obj <- new("EList", list(E = y, weights = weightsMatrix))

  # Fit regression model
  #---------------------

  # use pre-specified weights, if available

  if (.isMixedModelFormula(formula)) {
    # fit linear mixed model
    vpList <- fitVarPartModel(obj, formula, data, rescaleWeights = TRUE, ..., fxn = function(fit) {
      # extract
      # 1) sqrt residual variance (i.e. residual standard deviation)
      # 3) fitted values
      list(
        sd = sigma(fit),
        Amean = mean(fit@frame[, 1], na.rm = TRUE),
        fitted.values = predict(fit)
      )
    }, BPPARAM = BPPARAM)

    fit <- list()
    fit$sigma <- sapply(vpList, function(x) x$sd)
    fit$Amean <- sapply(vpList, function(x) x$Amean)
    fit$df.residual <- rep(2, length(fit$sigma)) # check this

    # extract fitted values
    fitted.values <- lapply(vpList, function(x) x$fitted.values)
    fitted.values <- do.call("rbind", fitted.values)
  } else {
    design <- model.matrix(formula, data)

    # fit <- lmFit(obj$E,design,weights=weights,...)
    # Use weights included in EList
    fit <- lmFit(obj, design, ...)

    if (fit$rank < ncol(design)) {
      j <- fit$pivot[1:fit$rank]
      fitted.values <- fit$coef[, j, drop = FALSE] %*% t(fit$design[, j, drop = FALSE])
    } else {
      fitted.values <- fit$coef %*% t(fit$design)
    }
  }

  # only keep y where model has converged
  keepGenes = rownames(fitted.values)
  y <- y[keepGenes, , drop = FALSE]

  if (is.null(fit$Amean)) fit$Amean <- rowMeans(y, na.rm = TRUE)

  # 		If no replication found, set all weight to 1
  NWithReps <- sum(fit$df.residual > 0L)
  if (NWithReps < 2L) {
    if (NWithReps == 0L) warning("The experimental design has no replication. Setting weights to 1.")
    if (NWithReps == 1L) warning("Only one gene with any replication. Setting weights to 1.")
    out$E <- y
    out$weights <- y
    out$weights[] <- 1
    if (!is.null(design)) out$design <- design
    if (is.null(out$targets)) {
      out$targets <- data.frame(lib.size = lib.size)
    } else {
      out$targets$lib.size <- lib.size
    }
    return(new("EList", out))
  }

  # Fit lowess trend to sqrt-standard-deviations by log-count-size
  sx <- fit$Amean + mean(log2(lib.size + 1)) - log2(1e6)

  # get residual standard deviation
  sy <- sqrt(fit$sigma)

  allzero <- rowSums2(counts) == 0
  if (any(allzero)) {
    sx <- sx[!allzero]
    sy <- sy[!allzero]
  }

  if (span == "auto") {
    # fit loess with automatic select of tuning parameter
    fit <- fANCOVA::loess.as(sx, sy)

    # sort by x value
    i <- order(fit$x)
    l <- list(x = fit$x[i], y = fit$fitted[i])
  } else {
    l <- stats::lowess(sx, sy, f = span)
  }

  # 	Make interpolating rule
  suppressWarnings({
    f <- approxfun(l, rule = 2)
  })

  if (plot) {
    plot(sx, sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )", pch = 16, cex = 0.25)
    title("voom: Mean-variance trend")
    lines(l, col = "red")
  }

  # apply interpolation to data
  fitted.cpm <- 2^fitted.values
  fitted.count <- 1e-6 * t(t(fitted.cpm) * (lib.size + 1))
  fitted.logcount <- log2(fitted.count)

  # 	Apply trend to individual observations
  w <- 1 / f(fitted.logcount)^4
  dim(w) <- dim(fitted.logcount)

  # 	Output
  out$E <- y
  out$weights <- w
  if (!is.null(design)) out$design <- design
  if (is.null(out$targets)) {
    out$targets <- data.frame(lib.size = lib.size)
  } else {
    out$targets$lib.size <- lib.size
  }

  # rescale by input weights
  if( rescaleWeightsAfter ){  
    w <- weightsMatrix[keepGenes,,drop=FALSE]

    out$weights <- w * out$weights
    out$targets$sample.weights <- colMeans(weightsMatrix)
  }
  
  # remove rownames to be compatible with voom()
  attr(out$weights, "dimnames") = NULL

  if (save.plot) {
    out$voom.xy <- list(x = sx, y = sy, xlab = "log2( count size + 0.5 )", ylab = "Sqrt( standard deviation )")
    out$voom.line <- l
  }

  # Check max value of precision weights
  maxValue <- max(out$weights)
  if (maxValue > 1e8) {
    txt <- paste0("The maximum precision weight is ", format(maxValue, scientific = TRUE), ", suggesting a poor smoothing fit\non the mean-variance plot for large expression values. Such large weights can\nhave unexpected effects downstream.  Consider examining the mean-variance plot\nand reducing the span parameter.")
    warning(txt)
  }

  out <- new("EList", out)

  if (.isMixedModelFormula(formula)) {
    # attach errors if they exist
    attr(out, "errors") <- attr(vpList, "errors")
  }

  out
}
