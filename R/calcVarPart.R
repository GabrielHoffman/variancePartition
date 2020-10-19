# Gabriel Hoffman
# 
# October 13, 2020
# Simplify calculations of variance fractions and 
# add compatability with glm's


#' Compute variance statistics
#' 
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#'
#' @param fit model fit from lm() or lmer()
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param ... additional arguments (not currently used)
#' 
#' @return
#' fraction of variance explained / ICC for each variable in the linear model 
#' 
#' @description
#' For linear model, variance fractions are computed based on the sum of squares explained by each component.  For the linear mixed mode, the variance fractions are computed by variance component estimates for random effects and sum of squares for fixed effects.  
#'
#' For a generalized linear model, the variance fraction also includes the contribution of the link function so that fractions are reported on the linear (i.e. link) scale rather than the observed (i.e. response) scale. For linear regression with an identity link, fractions are the same on both scales.  But for logit or probit links, the fractions are not well defined on the observed scale due to the transformation imposed by the link function.  
#'
#' The variance implied by the link function is the variance of the corresponding distribution:
#' 	logit -> logistic distribution -> variance is pi^2/3
#'	probit -> standard normal distribution -> variance is 1
#'
#' Reviewed by
#' Nakagawa and Schielzeth. 2012. A general and simple method for obtaining R2 from generalized linear mixed-effects models.  https://doi.org/10.1111/j.2041-210x.2012.00261.x
#'
#' Proposed by
#' McKelvey and Zavoina. A statistical model for the analysis of ordinal level dependent variables. The Journal of Mathematical Sociology 4(1) 103-120 https://doi.org/10.1080/0022250X.1975.9989847
#'
#' Also see
#' DeMaris. Explained Variance in Logistic Regression: A Monte Carlo Study of Proposed Measures. Sociological Methods & Research 2002 https://doi.org/10.1177/0049124102031001002
#'
#' We note that Nagelkerke's pseudo R^2 evaluates the variance explained by the full model.  Instead, a variance partitioning approach evaluates the variance explained by each term in the model, so that the sum of each systematic plus random term sums to 1 (Hoffman and Schadt, 2016, Nakagawa and Schielzeth, 2012).
#'
#' @examples
#' library(lme4)
#' data(varPartData)
#'
#' # Linear mixed model
#' fit <- lmer( geneExpr[1,] ~ (1|Tissue) + Age, info)
#' calcVarPart( fit )
#'
#' # Linear model
#' # Note that the two models produce slightly different results
#' # This is expected: they are different statistical estimates 
#' # of the same underlying value
#' fit <- lm( geneExpr[1,] ~ Tissue + Age, info)
#' calcVarPart( fit )
#' 
#' @export
#' @docType methods
#' @rdname calcVarPart-method
setGeneric("calcVarPart", signature="fit",
  function(fit, showWarnings=TRUE, ...)
      standardGeneric("calcVarPart")
)


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lm-method
setMethod("calcVarPart", "lm",
function(fit, showWarnings=TRUE, ...){

	# check validity of model fit
	checkModelStatus( fit, showWarnings, ...)

	an = anova(fit)
	varFrac = an[['Sum Sq']] / sum( an[['Sum Sq']] )
	names(varFrac) = rownames(an)
	varFrac
})

#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lmerMod-method
setMethod("calcVarPart", "lmerMod",
function(fit, showWarnings=TRUE,...){

	# check validity of model fit
	checkModelStatus( fit, showWarnings, ...)

	# extract variance components
	vc = unlist(getVarianceComponents(fit))

	# create fractions
	varFrac = vc / sum(vc)

	# remove ".(Intercept)" string
	names(varFrac) = gsub("\\.\\(Intercept\\)", "", names(varFrac))

	varFrac
})


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,glm-method
setMethod("calcVarPart", "glm",
function(fit, showWarnings=TRUE, ...){

	N = nrow(fit$model)

	# Compute eta for each term
	# predicted value in linear space for each term
	Eta = predict(fit, type="terms")

	# compute residual term for each link
	distVar = switch( fit$family$link , 
			"logit" = (pi^2)/3,
			"probit" = 1,
			"identity" = var(fit$residuals) )

	if( is.null(distVar) ){
		stop("glm link not supported: ", fit$family$link )
	}

	# variance on linear scale
	# get variance of each term
	# append with variance due to link function
	var_term = c(apply(Eta, 2, var), distVar)

	# compute fraction
	varFrac_linear = var_term / sum( var_term )
	names(varFrac_linear) = c(colnames(Eta), "Residuals")

	return( varFrac_linear )
})


getVarianceComponents = function( fit ){
	
	varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)

	# order variables by name
	# essential so that all models are ordered the same
	varComp = varComp[order(names(varComp))]

	# extract predictor for each fixed effect
	idx = which(colnames(fit@pp$X) != "(Intercept)")

	# if there are fixed effects
	if( length(idx) > 0){

		# this part is now in 1.5.2
		# better estimates of fixed effects
		fxeff = sapply( idx, function(i){
			fit@pp$X[,i] * lme4::fixef(fit)[i]
		})
		colnames(fxeff) = colnames(fit@pp$X)[idx]

		# compute variance of each fixed effect
		N = nrow(fxeff)

		# variance of eahc term
		fixedVar = apply(fxeff, 2, var) * (N-1) / N

		# variance of sum of terms
		varFixedTotal = var(rowSums(fxeff)) * (N-1) / N

		# scale variance by the sum  and the total variance
		for( key in names(fixedVar)){
			varComp[[key]] = as.array(fixedVar[[key]] / sum(fixedVar) * varFixedTotal)
		}
	}

	# get residuals
	varComp$Residuals = attr(lme4::VarCorr(fit), 'sc')^2
	names(varComp$Residuals) = ''
	
	return(varComp)
}
