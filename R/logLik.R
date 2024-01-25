

#' Log-likelihood from model fit
#'
#' Log-likelihood from model fit
#'
#' @param object result of \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} used to fit model
#' @param ... See \code{?stats::logLik}
#'
#' @rawNamespace S3method("logLik", MArrayLM2)
#' @importFrom stats logLik
#' @export
logLik.MArrayLM = function(object, vobj, ...){

	if( !is.null(object$logLik) ){
		return( object$logLik )
	}

	if( is.null(object$residuals) ){
		object$residuals = residuals.MArrayLM(object, vobj)
	}

	values = sapply(seq(nrow(object)), function(i){
		obj = list(residuals = object$residuals[i,,drop=TRUE], 
			rank = object$rank, 
			weights = vobj$weights[i,,drop=TRUE])
		
		class(obj) = "lm"

		logLik( obj )
	})
	names(values) = rownames(object)

	values
}

#' Log-likelihood from model fit
#'
#' Log-likelihood from model fit
#'
#' @param object result of \code{lmFit()} or \code{dream()}
#' @param ... See \code{?stats::logLik}
#'
#' @rawNamespace S3method("logLik", MArrayLM)
#' @export
logLik.MArrayLM2 = function(object, ...){

	object$logLik
}

#' BIC from model fit
#'
#' BIC from model fit
#'
#' @param object result of \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} used to fit model
#' @param ... See \code{?stats::BIC}
#'
#' @rawNamespace S3method("BIC", MArrayLM2)
#' @importFrom stats BIC
#' @export
BIC.MArrayLM = function(object, vobj, ...){

	if( !is.null(object$BIC) ){
		return( object$BIC )
	}

	if( is.null(object$residuals) ){
		object$residuals = residuals.MArrayLM(object, vobj)
	}

	# n = nrow(object$design)
	# df = object$rank + 1

	values = sapply(seq(nrow(object)), function(i){
		obj = list(residuals = object$residuals[i,,drop=TRUE], 
			rank = object$rank, 
			weights = vobj$weights[i,,drop=TRUE])
		
		class(obj) = "lm"

		BIC(obj)
		# ll = logLik( obj )

		# browser()
		# -2*ll + df * log(n)
	})
	names(values) = rownames(object)

	values
}

#' BIC from model fit
#'
#' BIC from model fit using edf
#'
#' @param object result of \code{dream()}
#' @param vobj \code{EList} used to fit model
#' @param ... See \code{?stats::BIC}
#'
#' @rawNamespace S3method("BIC", MArrayLM)
#' @export
BIC.MArrayLM2 = function(object, vobj, ...){

	if( !is.null(object$BIC) ){
		return( object$BIC )
	}

	df = object$edf
	n = ncol(object$residuals)
	ll = object$logLik

	-2*ll + df * log(n)
}












