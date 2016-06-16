

#' Compute variance statistics
#' 
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#'
#' @param fit model fit from lm() or lmer()
#' @param adjust remove variation from specified variables from the denominator.  This computes the adjusted ICC with respect to the specified variables
#' @param adjustAll adjust for all variables.  This computes the adjusted ICC with respect to all variables
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param ... additional arguments (not currently used)
#' 
#' @return
#' fraction of variance explained / ICC for each variable in the model
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
  function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
      standardGeneric("calcVarPart")
)


# calcVarPart <- function( fit ){
# 	if( is.null(fit) ){
# 		stop("model fit is NULL.  There was an error fitting the regression model")
# 	}
# 	UseMethod("calcVarPart")
# }


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lm-method
setMethod("calcVarPart", "lm",
function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
{

	# check validity of model fit
	checkModelStatus( fit, showWarnings,...)

	# Get ANOVA
	a = anova(fit) 

	# get variables to remove from denominator
	# also, check variables
	adjust = getAdjustVariables( rownames(a), adjust, adjustAll)

	# get Sum of Squares
	if( is.null(adjust) ){
		varFrac = a[['Sum Sq']] / sum( a[['Sum Sq']] )
	}else{

		v = a[['Sum Sq']]
		names(v) = rownames(a)

		varFrac = c()
		for( i in 1:(length(v)-1) ){
			# total variance minus components to remove, but add back the current variable
			varFrac[i] = v[i] / get_denom( v, adjust, i)
		}
		# Residual variance
		varFrac[length(v)] = v[["Residuals"]] / sum(v)

	}

	# set names 
	names(varFrac) = rownames(a)

	return( varFrac )
}
)

get_denom = function( v, adjust, i){

	currentSet = setdiff(adjust, names(v)[i])

	if( length(currentSet) > 0){
		denom = sum(v) - sum(sapply(currentSet, function(x) v[x]))
	}else{
		denom = sum(v)
	}
	return(denom)
}

isMultipleVaryingCoefficientTerms = function( fit ){

	# get variance components values
	varComp = getVarianceComponents( fit )

	# get varying coefficient terms
	res = which(sapply(varComp, length) > 1)

	# if there are more than 1
	if( length(res) > 1){
		warning(paste("Cannot have more than one varying coefficient term:", paste(names(res), collapse=", "), "\nThe results will not behave as expected and may be very wrong!!"))
	}
}

#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lmerMod-method
setMethod("calcVarPart", "lmerMod",
function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE,...)
{
	# check validity of model fit
	checkModelStatus( fit, showWarnings, ...)

	# get variance components values
	varComp = getVarianceComponents( fit )

	# get variables to remove from denominator
	# also, check variables
	adjust = getAdjustVariables( names(varComp), adjust, adjustAll)

	if( max(sapply(varComp, length)) > 1 && !is.null(adjust) ){
		stop("The adjust and adjustAll arguments are not currently supported for varying coefficient models")
	}

	variableLevels = list()
	for(key in colnames(fit@frame)){
		key2 = paste(paste(key, levels(fit@frame[[key]]), sep=''), collapse=',')
		variableLevels[[key2]] = key
	}

 	# Get variance terms for each variable
 	# for standard terms, return the variance
 	# for varying coefficient terms return weighted sum of variance
 	# 	weighted by the sample sizes corresponding to the subsets of the data
	get_total_variance = function(varComp){
		sapply(varComp, function(x){
		if( length(x) == 1){
			return( x )
		}else{			
			key = variableLevels[[paste(names(x), collapse=',')]]

			# get total variance from varying coefficient model
			weights = (table(fit@frame[[key]]) / nrow(fit@frame))
			x %*% weights / sum(weights)
		}
		})
	}

	varPart = c()

	for( key in names(varComp) ){
		for( i in 1:length(varComp[[key]]) ){

			# get variance contributed by other variables
			varOtherVariables = varComp[-which(names(varComp) == key)]

			# remove this variance from the denominator
			currentSet = setdiff(adjust, key)

			if( length(currentSet) > 0 && key != "Residuals"){
				adjustVariance = sum(sapply(currentSet, function(x) varComp[[x]]))
			}else{
				adjustVariance = 0
			}

			# demoninator: remove variables in this class (i.e. same key)
			# from the variables in the current class, only consider 1 variable
			totVar = tryCatch({
			     get_total_variance(varOtherVariables)
			}, error = function(e) {
			    stop("Problem with varying coefficient model in formula: should have form (A+0|B)")
			}, finally = {
			})
			denom = sum(totVar) + varComp[[key]][i] - adjustVariance

			# compute fraction
			frac = varComp[[key]][i] / denom 

			# name variable based on two levels
			if( names(varComp[[key]])[i] %in% c("(Intercept)", '') ){
				names(frac) = key
			}else{
				names(frac) = paste( names(varComp[[key]])[i], key, sep='/')
			}

			# save result
			varPart = c(varPart, frac)
		}
	}

	return( varPart )
}
)


#' Extract variance terms
#' 
#' Extract variance terms from a model fit with lm() or lmer()
#'
#' @param fit list of lmer() model fits
#'  
#' @return 
#'  variance explained by each variable
# @details
#' @examples
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' cl <- makeCluster(4)
#' registerDoParallel(cl)
#' # or by using the doSNOW package
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' modelList <- fitVarPartModel( geneExpr, form, info )
#' 
#' fit <- modelList[[1]]
#' getVarianceComponents( fit )
#' 
#' # stop cluster
#' stopCluster(cl)
#' 
#' @export
getVarianceComponents = function( fit ){
	# compute ICC, but don't divide by variances in the same class
	varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)

	# order variables by name
	# essential so that all models are ordered the same
	varComp = varComp[order(names(varComp))]

	# variance of fixed effects, if they exist
	# including fixed effects in the sum of variances makes a big difference, 
	#	especially when the fixed effect makes a big contribution
	# this approach minimized the difference between modelling
	# 	as a fixed or random effect
	# although, there is many be a substantial difference in estimates

	idx = which(colnames(fit@pp$X) != "(Intercept)")
	if( length(idx) > 0){

		# get predicted fixed effects
		fxeff = sapply( idx, function(i){
			fit@pp$X[,i] * lme4::fixef(fit)[i]
		})
		colnames(fxeff) = colnames(fit@pp$X)[idx]

		fixedVar = apply(fxeff, 2, var)

		for( i in 1:length(fixedVar) ){

			key = names(fixedVar)[i]
			varComp[[key]] = fixedVar[i]
			names(varComp[[key]]) = ''
		}

	}

	# get residual variance
	varComp$Residuals = attr(lme4::VarCorr(fit), 'sc')^2
	names(varComp$Residuals) = ''

	return( varComp )
}

# Construct list of variable names to remove from the denominator
getAdjustVariables = function( variables, adjust, adjustAll){

	if( adjustAll ){
		adjust = variables
	}

	# Check that variables in adjust array are present
	idx = match(adjust, variables)

	if( any(is.na(idx)) ){
		stop(paste("The following variables cannot be used in the 'adjust' argument\nbecause they are not present in the model fit:", paste(adjust[which(is.na(idx))], collapse='; ')))
	}

		# remove Residuals from the adjust list
	if( "Residuals" %in% adjust){
		adjust = adjust[-match("Residuals", adjust)]
	}

	return( adjust )
}
