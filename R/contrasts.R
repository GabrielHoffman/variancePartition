
#' Extract contrast matrix for linear mixed model
#' 
#' Extract contrast matrix, L, testing a single variable.  Contrasts involving more than one variable can be constructed by modifying L directly
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula 
#' @param coefficient the coefficient to use in the hypothesis test
#' 
#' @return
#' Contrast matrix testing one variable
#'
#' @examples
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # get contrast matrix testing if the coefficient for Batch2 is zero 
#' # The variable of interest must be a fixed effect
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L = getContrast( geneExpr, form, info, "Batch3")
#'
#' # get contrast matrix testing if Batch3 - Batch2 = 0
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L = getContrast( geneExpr, form, info, c("Batch3", "Batch2"))
#'
#' # To test against Batch1 use the formula:
#' # 	~ 0 + Batch + (1|Individual) + (1|Tissue) 
#' # to estimate Batch1 directly instead of using it as the baseline
#'
#' @export
# @docType methods
#' @rdname getContrast-method
getContrast = function( exprObj, formula, data, coefficient){ 

	if( length(coefficient) > 2){
		stop("Length of coefficient array limited to 2")
	} else if (length(coefficient) == 0) {
		stop("Need at least one coefficient")
	} else if (any(is.na(coefficient))) {
		stop("Coefficient must not be NA")
	}

	L = .getContrastInit( formula, data)

	# assign coefficient coding
	if( any(!coefficient %in% names(L)) ){
		stop("coefficient is not in the formula.  Valid coef are:\n", paste(names(L), collapse=', '))
	}
	L[coefficient[1]] = 1

	if( length(coefficient) == 2){			
		L[coefficient[2]] = -1
	}
	
	L
}


#' Construct Matrix of Custom Contrasts
#' 
#' Construct the contrast matrix corresponding to specified contrasts of a set of parameters.
#'
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula 
#' @param ... expressions, or character strings which can be parsed to expressions, specifying contrasts
#' @param contrasts character vector specifying contrasts
#' 
#' @return
#' matrix of linear contrasts between regression coefficients
#'
#' @details
#' This function expresses contrasts between a set of parameters as a numeric matrix. The parameters are usually the coefficients from a linear (mixed) model fit, so the matrix specifies which comparisons between the coefficients are to be extracted from the fit. The output from this function is usually used as input to \code{dream()}.
#'
#' This function is inspired by \code{limma::makeContrasts()} but is designed to be compatible with linear mixed models for \code{dream()}
#'
#' Names in ... and contrasts will be used as column names in the returned value.
#'
#' @examples
#' # load library
#' # library(variancePartition)
#' 
#' # Intialize parallel backend with 4 cores
#' library(BiocParallel)
#' register(SnowParam(4))
#' 
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' form <- ~ 0 + Batch + (1|Individual) + (1|Tissue) 
#' 
#' # Define contrasts
#' L = makeContrastsDream( form, info, contrasts = c(Batch1_vs_2 = "Batch1 - Batch2", Batch3_vs_4 = "Batch3 - Batch4", Batch1_vs_34 = "Batch1 - (Batch3 + Batch4)/2"))
#' 
#' # show contrasts matrix
#' L
#' 
#' # Plot to visualize contrasts matrix
#' plotContrasts(L)
#' 
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit = dream( geneExpr[1:10,], form, info, L=L)
#' 
#' # examine contrasts after fitting
#' head(coef(fit))
#' 
#' # show results from first contrast
#' topTable(fit, coef="Batch1_vs_2")
#' 
#' # show results from second contrast
#' topTable(fit, coef="Batch3_vs_4")
#' 
#' # show results from third contrast
#' topTable(fit, coef="Batch1_vs_34")
#'
#' @importFrom rlang new_environment eval_tidy caller_env
#' @export
makeContrastsDream = function(formula, data, ..., contrasts=NULL){
  e <- .getContrastExpressions(..., contrasts = contrasts)
	if( length(e) == 0 ){
		stop("Must specify at least one contrast")
	}
  L_uni <- .getAllUniContrasts(formula, data)
  L_uni_env <- new_environment(
    c(asplit(L_uni, 2)),
    caller_env()
  )
  L <- do.call(cbind, lapply(e, eval_tidy, env = L_uni_env))
  rownames(L) <- rownames(L_uni)
  names(dimnames(L)) <- c("Levels", "Contrasts")

  # detect univariate contrasts
  anyUnivariate = apply(L, 2, function(x) sum(x!=0))

  if( any(anyUnivariate == 1) ){
    txt = paste("All univariate contrasts are already included.\n  Manually specifying them here can cause issues downstream.\n  Terms: ", paste0(names(which(anyUnivariate == 1)), collapse = ', '))
    warning(txt)
  }

  L
}

#' @importFrom rlang enexprs parse_expr
.getContrastExpressions = function(..., contrasts) {
  e <- enexprs(...)
  if (!missing(contrasts) && !is.null(contrasts) && length(contrasts) > 0) {
    if (length(e) > 0) {
      stop("Can't specify both ... and contrasts")
    }
    e <- lapply(as.character(unlist(contrasts)), parse_expr)
    names(e) <- names(unlist(contrasts))
  }
  # Allow contrasts to be specified as string literals in addition to
  # unquoted expressions.
  for (i in seq_along(e)) {
    if (is.character(e[[i]])) {
      e[[i]] <- parse_expr(e[[i]])
    }
  }
  e_text <- vapply(e, deparse1, character(1))
  if (is.null(names(e))) {
    names(e) <- e_text
  } else {
    empty_names <- is.na(names(e)) | names(e) == ""
    names(e)[empty_names] <- e_text[empty_names]
  }
  e
}

#' @importFrom lme4 nobars
.getFixefNames = function( formula, data, ... ) {
  if (!isRunableFormula(, formula, data)) {
    stop("the fixed-effects model matrix is column rank deficient")
  }
  ## See lme4::lFormula
  formula[[length(formula)]] <- nobars(formula[[length(formula)]])
  colnames(model.matrix(object = formula, data = data, ...))
}

#' Get all univariate contrasts
#'
#' Get all univariate contrasts
#'
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula 
#'
#' @return
#'  Matrix testing each variable one at a time.  Contrasts are on rows
#'
#' @keywords internal
.getAllUniContrasts = function( formula, data){
  fixef_names <- .getFixefNames(formula, data)
  ## For large designs, might be worth using Matrix::Diagonal to make
  ## a sparse diagonal matrix
  L <- diag(x = 1, nrow = length(fixef_names))
  dimnames(L) <- list(
    Levels = fixef_names,
    Contrasts = fixef_names
  )
  L
}

.getContrastInit = function( formula, data ){
  fixef_names <- .getFixefNames(formula, data)
  L <- rep(0, length(fixef_names))
  names(L) <- fixef_names
  L
}