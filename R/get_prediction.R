# Gabriel Hoffman
# June 2, 2020
#

#' Compute predicted value of formula for linear mixed model
#'
#' Compute predicted value of formula for linear mixed model for with lmer
#'
#' @param fit model fit with lmer
#' @param formula formula of fixed and random effects to predict
#'
#' @return Predicted values from formula using parameter estimates from fit linear mixed model
#'
#' @details Similar motivation as \code{lme4:::predict.merMod()}, but that function cannot use just a subset of the fixed effects: it either uses none or all.  Note that the intercept is included in the formula by default.  To exclude it from the prediction use \code{~ 0 + ...} syntax
#' 
#' @examples
#' 
#' library(lme4)
#' 
#' # fit model
#' fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
#' 
#' # predict Days, but exclude intercept
#' get_prediction( fm1, ~ 0 + Days)
#' 
#' # predict Days and (Days | Subject) random effect, but exclude intercept
#' get_prediction( fm1, ~ 0 + Days +  (Days | Subject))
#' 
#' @import lme4
#' @export
get_prediction = function( fit, formula){

	# initialize to zeros
	pred_fixed = pred_rand = rep(0, length(fit@resp$y))

	# if a random effect is specified
	if( !is.null(findbars(formula)) ){
		# RANDOM
		#-------
		# Get sum of BLUP's for the specified random effects

		ran_form = lme4:::reOnly(formula)

		pred_rand = predict(fit, re.form=ran_form, random.only=TRUE)
	}

	# FIXED 
	#------
	# Compute X\beta for the specified set of fixed effects

	# formula of only fixed effects
	fixed_form = nobars(formula)

	# create design matrix from formula
	# dsgn = lme4::getME(fit, "X")
	dsgn = model.matrix( fixed_form, fit@frame)

	# only extract coefficients matching design 
	beta = fixef(fit)[colnames(dsgn)]

	if( length(beta) > 0){		
		# prediction of fixed effects
		pred_fixed = as.numeric(dsgn %*% beta)
	}

	# combine
	y_pred = pred_rand + pred_fixed
	names(y_pred) = rownames(dsgn)
	y_pred
}