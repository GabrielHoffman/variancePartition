#' Extract contrast matrix for linear mixed model
#'
#' Extract contrast matrix, L, testing a single variable.  Contrasts involving more than one variable can be constructed by modifying L directly
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
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
#' form <- ~ Batch + (1 | Individual) + (1 | Tissue)
#' L <- getContrast(geneExpr, form, info, "Batch3")
#'
#' # get contrast matrix testing if Batch3 - Batch2 = 0
#' form <- ~ Batch + (1 | Individual) + (1 | Tissue)
#' L <- getContrast(geneExpr, form, info, c("Batch3", "Batch2"))
#'
#' # To test against Batch1 use the formula:
#' # 	~ 0 + Batch + (1|Individual) + (1|Tissue)
#' # to estimate Batch1 directly instead of using it as the baseline
#'
#' @export
# @docType methods
#' @rdname getContrast-method
getContrast <- function(exprObj, formula, data, coefficient) {
  if (length(coefficient) > 2) {
    stop("Length of coefficient array limited to 2")
  } else if (length(coefficient) == 0) {
    stop("Need at least one coefficient")
  } else if (any(is.na(coefficient))) {
    stop("Coefficient must not be NA")
  }
  formula <- as.formula(formula)

  L <- .getContrastInit(formula, data)

  # assign coefficient coding
  if (any(!coefficient %in% names(L))) {
    stop("coefficient is not in the formula.  Valid coef are:\n", paste(names(L), collapse = ", "))
  }
  L[coefficient[1]] <- 1

  if (length(coefficient) == 2) {
    L[coefficient[2]] <- -1
  }

  L
}


#' Construct Matrix of Custom Contrasts
#'
#' Construct the contrast matrix corresponding to specified contrasts of a set of parameters. Each specified set of contrast weights must sum to 1.
#'
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula
#' @param ... expressions, or character strings which can be parsed to expressions, specifying contrasts
#' @param contrasts character vector specifying contrasts
#' @param suppressWarnings (default FALSE). suppress warnings for univariate contrasts
#' @param nullOnError (default FALSE). When a contrast entry is invalid, throw warning and return NULL for that contrast entry
#'
#' @return
#' matrix of linear contrasts between regression coefficients
#'
#' @details
#' This function expresses contrasts between a set of parameters as a numeric matrix. The parameters are usually the coefficients from a linear (mixed) model fit, so the matrix specifies which comparisons between the coefficients are to be extracted from the fit. The output from this function is usually used as input to \code{dream()}.
#'
#' This function creates a matrix storing the contrasts weights that are applied to each coefficient.
#'
#' Consider a variable \code{v} with levels \code{c('A', 'B', 'C')}.  A contrast comparing \code{A} and \code{B} is \code{'vA - vB'} and tests whether the difference between these levels is different than zero. Coded for the 3 levels this has weights \code{c(1, -1, 0)}.  In order to compare \code{A} to the other levels, the contrast is \code{'vA - (vB + vC)/2'} so that \code{A} is compared to the average of the other two levels. This is encoded as \code{c(1, -0.5, -0.5)}.  This type of proper matching in testing multiple levels is enforced by ensuring that the contrast weights sum to 1. Based on standard regression theory only weighted sums of the estimated coefficients are supported.
#'
#' This function is inspired by \code{limma::makeContrasts()} but is designed to be compatible with linear mixed models for \code{dream()}
#'
#' Names in ... and contrasts will be used as column names in the returned value.
#'
#' @examples
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
#' form <- ~ 0 + Batch + (1 | Individual) + (1 | Tissue)
#'
#' # Define contrasts
#' # Note that for each contrass, the weights sum to 1
#' L <- makeContrastsDream(form, info, contrasts = c(Batch1_vs_2 = "Batch1 - Batch2", Batch3_vs_4 = "Batch3 - Batch4", Batch1_vs_34 = "Batch1 - (Batch3 + Batch4)/2"))
#'
#' # show contrasts matrix
#' L
#'
#' # Plot to visualize contrasts matrix
#' plotContrasts(L)
#'
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit <- dream(geneExpr[1:10, ], form, info, L = L)
#'
#' # examine contrasts after fitting
#' head(coef(fit))
#'
#' # show results from first contrast
#' topTable(fit, coef = "Batch1_vs_2")
#'
#' # show results from second contrast
#' topTable(fit, coef = "Batch3_vs_4")
#'
#' # show results from third contrast
#' topTable(fit, coef = "Batch1_vs_34")
#'
#' @importFrom rlang new_environment eval_tidy caller_env
#' @seealso \code{plotContrasts()}
#' @export
makeContrastsDream <- function(formula, data, ..., contrasts = NULL, suppressWarnings = FALSE, nullOnError = FALSE) {
  if (all(is.na(contrasts))) contrasts <- NULL

  # address special case of matching terms with colon (:)
  # since colon is implicitly converted to dot internally somewhere
  contrasts <- gsub(":", ".", contrasts)

  # droplevels
  data <- as.data.frame(data)
  data <- droplevels(data)

  e <- .getContrastExpressions(..., contrasts = contrasts)

  if (length(e) == 0) {
    stop("Must specify at least one contrast")
  }

  formula <- as.formula(formula)

  L_uni <- .getAllUniContrasts(formula, data)

  # address special case of matching terms with colon (:)
  colnames(L_uni) <- gsub(":", ".", colnames(L_uni))

  L_uni_env <- new_environment(
    c(asplit(L_uni, 2)),
    caller_env()
  )
  # L <- do.call(cbind, lapply(e, eval_tidy, env = L_uni_env))

  # add ability to drop problematic contrasts
  # if a contrast is invalid, throw warning and return null for just that entry
  eval_tidy_safe <- function(expr, data = NULL, env = caller_env(), nullOnError = TRUE) {
    if (nullOnError) {
      tryCatch(eval_tidy(expr, data, env), error = function(err) {
        warning(paste("makeContrastsDream:", err$message), immediate. = TRUE, call. = FALSE)
        NULL
      })
    } else {
      eval_tidy(expr, data, env)
    }
  }
  # comment out 12/13/23 to handle failure with only one contrast
  # nullOnError can only by TRUE if more than 1 contrast is specified
  # nullOnError <- nullOnError & (length(contrasts) > 1)
  res <- lapply(e, eval_tidy_safe, env = L_uni_env, nullOnError = nullOnError)
  L <- do.call(cbind, res)

  # check validity of contrasts
  if (!is.null(L)) {
    rownames(L) <- rownames(L_uni)
    names(dimnames(L)) <- c("Levels", "Contrasts")

    # Check that contrast coefs sum to 0
    #----
    s <- colSums(L)
    # Check even with floating point error
    # s!=0
    iseq <- abs(s) < sqrt(.Machine$double.eps)
    if (any(!iseq)) {
      txt <- paste0("Each contrast should sum to 0.\n  Not satisified for contrasts: ", paste(names(s)[!iseq], collapse = ", "))
      warning(txt, immediate. = TRUE)
    }
  }

  L
}

#' @importFrom rlang enexprs parse_expr
.getContrastExpressions <- function(..., contrasts) {
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

#' @importFrom reformulas nobars
.getFixefNames <- function(formula, data, ...) {
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
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula
#'
#' @return
#'  Matrix testing each variable one at a time.  Contrasts are on rows
#'
#' @keywords internal
.getAllUniContrasts <- function(formula, data) {
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

.getContrastInit <- function(formula, data) {
  fixef_names <- .getFixefNames(formula, data)
  L <- rep(0, length(fixef_names))
  names(L) <- fixef_names
  L
}
