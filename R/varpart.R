# July 9, 2026
# Gabriel Hoffman

# variance of fixed effect terms 
#' @importFrom stats cov terms
var_predict_terms <- function( formula, Beta, data, design ){

  if( missing(design) ){
    design <- model.matrix(formula, as.data.frame(data))
  }

  # Assign each coef to a contrast
  asgn <- attr(design, "assign")

  # Names of effects
  nmeffects <- attr(terms(formula), "term.labels")[unique(asgn)]

  if( attr(terms(formula),"intercept") == 1){
    nmeffects <- c("(Intercept)", nmeffects)
  }

  varEta <- sapply(seq(0, max(asgn)), function(i){
    # indeces in this term
    idx <- which(asgn == i)

    # Compute Eta for each response, then eval variance
    # Eta_i = tcrossprod(X[,idx,drop=FALSE], Beta[,idx,drop=FALSE])
    # colVars(Eta_i)

    # Using the fact that var(XB) = diag(B^T cov(X) B)
    Sigma <- cov(design[,idx,drop=FALSE])
    B <- Beta[,idx,drop=FALSE]
    # diag(B %*% Sigma %*% t(B))
    rowSums((B %*% Sigma) * B)
  })  
  varEta <- matrix(varEta, ncol=max(asgn)+1)
  rownames(varEta) <- rownames(Beta)
  colnames(varEta) <- nmeffects

  varEta
}


#' @importFrom fastglmm log_moments_nb_XB
.varpart <- function(formula, design, Beta, theta, offset, phi = 1, method = c("exact", "approximate"), pseudocount = 1, p.tail  = 1e-4, nthreads=1){

  method <- match.arg(method)

  # Mu <- design %*% t(Beta) + offset

  # Compute moments of log(y+c) for the NB distribution
  # faster in C++ code
  res <- log_moments_nb_XB(
    X       = design, 
    Beta    = Beta, 
    offset  = offset, 
    theta   = theta, 
    method  = method,
    c       = pseudocount, 
    p_tail  = p.tail,
    nthreads= nthreads)

  # scale noise variance by the QL dispersion scale
  res$var.noise <- res$var.noise * phi

  # total variance
  var.total <- res$var.signal + res$var.noise

  # total signal
  rho2.signal <- res$var.signal / var.total

  # Exact total noise
  rho2.noise <- res$var.noise / var.total

  # Compute variance fractions of the fixed effects
  eta_var <- var_predict_terms( formula = formula, Beta = Beta, design = design)

  # remove intercept since variance is zero
  eta_var <- eta_var[,colnames(eta_var) != "(Intercept)",drop=FALSE]
  gamma <- eta_var / rowSums(eta_var)

  if( any(is.nan(gamma)) ){
    gamma[is.nan(gamma)] <- 0
  }

  data.frame(
    rho2.signal*gamma, 
    CountNoise = rho2.noise*res$alpha,
    Residuals = rho2.noise*(1-res$alpha)
    )
}

#' Variance Partitioning Analysis
#'
#' Variance partitioning analysis on each gene
#' 
#' @param x regression model fit  
#' @param method select method for count models:  \code{"exact"} or \code{"approximate"} for faster approximation
#' @param pseudocount pseudocount used for \code{"exact"} and \code{"approximate"} methods for count models
#' @param p.tail probability threashold for evaluating expectations for \code{"exact"} methods for count models
#' @param nthreads number of threads used for count models
#' @param ... other arguments
#' 
#' @examples
#' library(ggplot2)
#' 
#' # Simulate counts
#' set.seed(1)
#' countMatrix <- matrix(rnbinom(n=100000, mu=20, size=3), ncol=10)
#' rownames(countMatrix) <- paste0("gene_", seq(nrow(countMatrix)))
#' colnames(countMatrix) <- paste0("sample_", seq(ncol(countMatrix)))
#' 
#' condition <- factor(rep(1:2, each=5))
#' 
#' # DESeq2 model #
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(countMatrix, 
#'   DataFrame(condition), 
#'   ~ condition)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' 
#' # Variance partition analysis
#' vp1 <- varpart(dds)
#' 
#' # Plot contribution of each component
#' plotVarPart(vp1, main="DESeq2") + theme(aspect.ratio=1)
#' 
#' # Plot count noise vs expression magnitude
#' plotTrendVP( dds, vp1, "CountNoise" )
#'
#' # edgeR model #
#' library(edgeR)
#' design <- model.matrix( ~ condition, data.frame(condition))
#' d <- DGEList(countMatrix)
#' d <- normLibSizes(d)
#' d <- estimateDisp(d, design)
#' fit <- glmQLFit(d, design)
#' fit <- glmQLFTest(fit)
#' 
#' vp2 <- varpart(fit, dispObj = d, formula = ~ cond)
#' 
#' # Plot contribution of each component
#' plotVarPart(vp2, main="edgeR") + theme(aspect.ratio=1)
#' 
#' # Plot count noise vs expression magnitude
#' plotTrendVP( dds, vp2, "CountNoise" )
#' 
#' @rdname varpart
#' @export
setGeneric("varpart", function(
  x,    
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  ...)  
  standardGeneric("varpart")
)

#' @rdname varpart
#' @importFrom stats model.matrix coef
#' @importFrom BiocGenerics sizeFactors design
#' @importFrom DESeq2 dispersions
#' @importFrom parallelly availableCores
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors mcols
#' @importFrom SummarizedExperiment colData
#' @export
setMethod("varpart", signature = "DESeqDataSet", 
  function(
  x,    
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  ...){

  method <- match.arg(method)

  # get overdispersion parameter
  theta <- 1 / dispersions(x)

  # create design matrix
  design <- model.matrix(design(x), colData(x))

  # Extract the beta coefficients (log2 fold changes)
  #  scale by log(2) to convert to natural log 
  Beta <- coef(x) * log(2)

  # get library size offset
  os <- log(sizeFactors(x))

  # keep only genes where theta is not NA
  include <- !is.na(theta)

  # if QL dispersion scale was estimated with glmGamPoi
  phi <- mcols(x)$qlDispMAP
  if( is.null(phi) ){
    phi <- 1
  }

  .varpart(
    formula     = design(x),
    design      = design, 
    Beta        = Beta[include,,drop=FALSE],
    theta       = theta[include], 
    offset      = os, 
    phi         = phi,
    method      = method,
    pseudocount = pseudocount, 
    p.tail      = p.tail,
    nthreads    = nthreads
    )
  }
)

#' @param dispObj result of \code{estimateDisp()}
#' @param formula formula used for the design matrix
#'
#' @rdname varpart
#' @export
setMethod("varpart", signature = "DGELRT", 
  function(
  x,    
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  dispObj,
  formula,
  ...){
  
  method <- match.arg(method)

  if( missing(dispObj) ){
    stop("For edgeR analysis, must specify dispObj as the result of estimateDisp()")
  }

  if( missing(formula) ){
    stop("For edgeR analysis, must specify formula used for the design matrix")
  }

  if( is.null(dispObj$tagwise.dispersion) ){
    stop("Dispersions not estimated, use estimateDisp() first")
  }

  .varpart(
    formula     = formula,
    design      = x$design, 
    Beta        = coef(x),
    theta       = 1 / dispObj$tagwise.dispersion, 
    offset      = c(x$offset), 
    phi         = x$s2.post, # QL dispersion scale
    method      = method,
    pseudocount = pseudocount, 
    p.tail      = p.tail,
    nthreads    = nthreads
    )
  }
)



#' Plot Trend of Variance Fractions Vs Count Magnitude
#'
#' Plot trend of variance fractions for a specified component versus count magnitude for each gene and cell cluster
#'
#' @param x object returned by \code{lucida()}
#' @param vp \code{data.frame} from \code{fitVarPart()}
#' @param component variance component to extract from \code{vp}
#' @param ... additional arguments
#'
#' @return Plot of variance fraction vs count magnitude
#'
#' @examples
#' # Simulate counts
#' set.seed(1)
#' countMatrix <- matrix(rnbinom(n=100000, mu=20, size=3), ncol=10)
#' rownames(countMatrix) <- paste0("gene_", seq(nrow(countMatrix)))
#' colnames(countMatrix) <- paste0("sample_", seq(ncol(countMatrix)))
#' 
#' condition <- factor(rep(1:2, each=5))
#' 
#' # DESeq2 model #
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(countMatrix, 
#'   DataFrame(condition), 
#'   ~ condition)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' 
#' # Variance partition analysis
#' vp1 <- varpart(dds)
#' 
#' # Plot count noise vs expression magnitude
#' plotTrendVP( dds, vp1, "CountNoise" )
#'
#'
#' # edgeR model #
#' library(edgeR)
#' design <- model.matrix( ~ condition, data.frame(condition))
#' d <- DGEList(countMatrix)
#' d <- normLibSizes(d)
#' d <- estimateDisp(d, design)
#' fit <- glmQLFit(d, design)
#' fit <- glmQLFTest(fit)
#' 
#' vp2 <- varpart(fit, dispObj = d, formula = ~ cond)
#' 
#' # Plot count noise vs expression magnitude
#' plotTrendVP( dds, vp2, "CountNoise" )
#' 
#' @rdname plotTrendVP-methods
#' @export
setGeneric(
  "plotTrendVP", 
  function(x, vp,component, ...){
    standardGeneric("plotTrendVP")
  }
)

#' @rdname plotTrendVP-methods
#' @importFrom DESeq2 results
#' @importFrom tibble rownames_to_column tibble
#' @importFrom dplyr inner_join `%>%` select filter
#' @importFrom ggplot2 sym ggplot scale_x_log10 geom_smooth theme_classic
#' @export
setMethod(
  "plotTrendVP", c("DESeqDataSet", "data.frame"),
  function(x, vp, component,...){

  ID <- baseMean <- NULL

  if( missing(component) ){
    stop("Must specify component")
  }

  if( ! component %in% colnames(vp) ){
    txt <- paste0("Requested component must be column in vp")
    stop(txt)
  }

  ylab <- paste("Variance explained by", component, "(%)")

  # extract results
  DESeq2::results(x) %>%
    data.frame %>%
    rownames_to_column("ID") %>%
    tibble %>%
    dplyr::select(ID, baseMean) %>%
    filter(baseMean > 0) %>%
    inner_join(vp %>% 
      rownames_to_column("ID"), by=c("ID")) %>%
    ggplot(aes(baseMean, 100*!!sym(component))) +
    geom_point() +
    scale_x_log10() +
    theme_classic() +
    theme(aspect.ratio=1, 
      strip.background = element_rect("grey95")) +
    xlab("Mean log10 counts") +
    scale_y_continuous(limits=c(0,100)) +
    ylab(ylab) +
    geom_smooth( method="nls", formula = y ~ SSlogis(x, Asym, xmid, scal), se=FALSE)
})

#' @rdname plotTrendVP-methods
#' @importFrom edgeR topTags
#' @export
setMethod(
  "plotTrendVP", c("DGELRT", "data.frame"),
  function(x, vp, component,...){

  ID <- logCPM <- NULL

  if( missing(component) ){
    stop("Must specify component")
  }

  if( ! component %in% colnames(vp) ){
    txt <- paste0("Requested component must be column in vp")
    stop(txt)
  }

  ylab <- paste("Variance explained by", component, "(%)")

  # extract results
  topTags(x, n=Inf) %>%
    data.frame %>%
    rownames_to_column("ID") %>%
    tibble %>%
    dplyr::select(ID, logCPM) %>%
    inner_join(vp %>% 
      rownames_to_column("ID"), by=c("ID")) %>%
    ggplot(aes(logCPM, 100*!!sym(component))) +
    geom_point() +
    theme_classic() +
    theme(aspect.ratio=1, 
      strip.background = element_rect("grey95")) +
    xlab("Mean log10 counts") +
    scale_y_continuous(limits=c(0,100)) +
    geom_smooth( method="nls", formula = y ~ SSlogis(x, Asym, xmid, scal), se=FALSE) +
    ylab(ylab)
})
