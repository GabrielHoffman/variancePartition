# Check that lm/lmer model is valid
# Throw warning if
# 	1) Intercept is ommited
# 	2) Any coefficient is NA
# 	3) a categorical variable is modeled as a fixed effect
setGeneric("checkModelStatus",
  signature = "fit",
  function(fit, dreamCheck = FALSE, ...) {
    standardGeneric("checkModelStatus")
  }
)

setMethod(
  "checkModelStatus", "lm",
  function(fit, dreamCheck = FALSE, ...) {
    # if no intercept is specified
    if (length(which(names(coef(fit)) == "(Intercept)")) == 0) {
      txt <- "No Intercept term was specified in the formula"
      stop(txt)
    }

    # if any coefficient is NA
    if (any(is.na(coef(fit)))) {
      stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
    }

    # Check condition number of covariance matrix
    condNum <- kappa(cov2cor(as.matrix(vcov(fit)))) # exact=TRUE

    if (condNum > 1e8) {
      stop(paste0("Condition number (", format(condNum, digits = 1), ") is very high.\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
    }
  }
)

setMethod(
  "checkModelStatus", "lmerMod",
  function(fit, dreamCheck = FALSE, ...) {
    run_model_check_mixed(fit, dreamCheck)
  }
)


setMethod(
  "checkModelStatus", "glmerMod",
  function(fit, dreamCheck = FALSE, ...) {
    run_model_check_mixed(fit, dreamCheck)
  }
)

#' @importFrom aod negbin
setMethod(
  "checkModelStatus", "negbin",
  function(fit, dreamCheck = FALSE, ...) {
    # run_model_check( fit, showWarnings, dream, colinearityCutoff, immediate )
  }
)

#' @importFrom lme4 isSingular
run_model_check_mixed <- function(fit, dreamCheck = FALSE) {
  # if no intercept is specified, give warning
  if (!dreamCheck && length(which(colnames(fit@pp$X) == "(Intercept)")) == 0) {
    txt <- "No Intercept term was specified in the formula."
    stop(txt)
  }

  # if any coefficient is NA
  if (dreamCheck && any(is.na(coef(fit)))) {
    stop("The variables specified in this model are redundant,\nthe design matrix is not full rank")
  }

  # Check condition number of covariance matrix
  condNum <- kappa(cov2cor(as.matrix(vcov(fit)))) # exact=TRUE

  if (condNum > 1e8) {
    stop(paste0("Condition number (", format(condNum, digits = 1), ") is very high.\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
  }

  # check that factors are random and continuous variables are fixed
  ###################################################################

  # remove backticks with gsub manually
  # solve issue that backticks are conserved is some but not all parts of lmer()

  # Simplified testing of random versus fixed effects
  # allows (A|B) only where A is continuous

  # variables fit by regression
  testVar <- attr(attr(fit@frame, "terms"), "term.labels")
  testVar <- gsub("`", "", testVar)

  # get type for each variable
  # keep only tested variables
  varType <- attr(attr(fit@frame, "terms"), "dataClasses")[-1]
  varType <- varType[testVar]

  # random effects
  randVar <- names(fit@flist)

  # fixed effects
  # starting with all variables, remove random variables
  fixedVar <- setdiff(testVar, randVar)

  for (i in 1:length(varType)) {
    # if factor is not random
    if (!dreamCheck && varType[i] %in% c("factor", "character") && (!names(varType)[i] %in% randVar)) {
      txt <- paste("Categorical variables modeled as fixed effect:", paste(names(varType)[i], collapse = ", "), "\n  Must model either _all_ or _no_ categorical variables as random effects here")
      stop(txt)
    }

    # If numeric/double is not fixed
    if (!dreamCheck && varType[i] %in% c("numeric", "double") && (!names(varType)[i] %in% fixedVar)) {
      stop(paste("Continuous variable cannot be modeled as a random effect:", names(varType)[i]))
    }
  }

  # show convergance message, if model is not singular
  if (!is.null(fit@optinfo$conv$lme4$messages) && !isSingular(fit)) {
    stop(fit@optinfo$conv$lme4$messages)
  }
}
