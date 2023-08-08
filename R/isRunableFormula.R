# Gabriel Hoffman
# April 6, 2021

#' Test if formula is full rank on this dataset
#'
#' Test if formula is full rank on this dataset
#'
#' @param exprObj expression object
#' @param formula formula
#' @param data data
#'
#' @importFrom lme4 lFormula lmerControl
#' @export
isRunableFormula <- function(exprObj, formula, data) {
  isRunable <- TRUE
  tryCatch(
    {
      control <- lmerControl(check.rankX = "stop.deficient")
      lFormula(formula = formula, data = data, control = control)
    },
    error = function(e) {
      mesg <- "the fixed-effects model matrix is column rank deficient"
      if (any(grepl(mesg, e$message))) {
        isRunable <<- FALSE
      }
    }
  )
  isRunable
}
