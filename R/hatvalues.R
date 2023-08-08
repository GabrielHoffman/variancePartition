#' Compute hatvalues
#'
#' Compute hatvalues from dream fit
#'
#' @param model model fit from \code{dream()}
#' @param vobj \code{EList} returned by \code{voom()} or \code{voomWithDreamWeights()}.
#' @param ... other arguments, currently ignored
#'
#' @export
#' @rdname hatvalues-method
#' @aliases hatvalues,MArrayLM-method
setMethod(
  "hatvalues", "MArrayLM",
  function(model, vobj, ...) {
    if (!is.null(model$hatvalues)) {
      return(model$hatvalues)
    }

    if (missing(vobj)) {
      stop("Must specify voom fit: vobj")
    }

    if (!identical(nrow(model), nrow(vobj))) {
      stop("model and vobj must have same number of rows")
    }

    # Note: rownames(model) gives NULL if there is 1 gene
    if ((nrow(model) > 1) && !identical(rownames(model), rownames(vobj))) {
      stop("rownames between model and vobj must match")
    }

    if (is(vobj, "EList")) {
      weights <- vobj$weights
    } else {
      weights <- matrix(1, nrow(vobj), ncol(vobj))
    }

    # With *unweighted* fixed effect model, hat values don't depend on response
    # But weights do affect the hat values
    # Compute hat values for each response
    hv <- lapply(seq(nrow(vobj)), function(j) {
      X <- sqrt(weights[j, ]) * model$design

      # compure full H matrix
      # H = X %*% solve(crossprod(X), t(X))
      # diag(H)
      # Just compute diagonals
      rowSums(X * t(solve(crossprod(X), t(X))))
    })
    hv <- do.call(rbind, hv)
    rownames(hv) <- rownames(vobj)

    hv
  }
)


#' @export
#' @rdname hatvalues-method
#' @aliases hatvalues,MArrayLM2-method
setMethod(
  "hatvalues", "MArrayLM2",
  function(model, ...) {
    model$hatvalues
  }
)
