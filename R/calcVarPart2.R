
#' Compute variance statistics
#'
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#'
#' @param fit regression model
#' @param ... other arguments passed to \code{fastglmm::varpart()}
#' 
#' @details this is now an iterface to \code{fastglmm::varpart()}
#'
#' @seealso \code{fastglmm::varpart()}
#' @examples
#' library(lme4)
#' data(varPartData)
#'
#' # Linear mixed model
#' fit <- lmer(geneExpr[1, ] ~ (1 | Tissue) + Age, info)
#' calcVarPart(fit)
#'
#' # Linear model
#' # Note that the two models produce slightly different results
#' # This is expected: they are different statistical estimates
#' # of the same underlying value
#' fit <- lm(geneExpr[1, ] ~ Tissue + Age, info)
#' calcVarPart(fit)
#
#' @export
#' @importFrom fastglmm varpart
setGeneric("calcVarPart",
  signature = "fit",
  function(fit, ...) {
    fastglmm::varpart(fit,...)
  }
)

