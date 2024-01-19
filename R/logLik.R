

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

	object$residuals = residuals.MArrayLM(object, vobj)

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
