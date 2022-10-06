# Gabriel Hoffman 
# October 1, 2022
# 
# Evaluate vcov() on result of dream() to 
# compute covariance between estimated coefficients

# TODO:
# Write vcov for two lm() or lmer() fits


#' Co-variance matrix for \code{dream()} fit
#'
#' Define generic \code{vcov()} for result of \code{lmFit()} and \code{dream()}
# 
#' @param object \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} object returned by \code{voom()}
#' @param coef name of coefficient to be extracted
#
#' @importFrom stats coefficients
#' @export
setMethod('vcov', c("MArrayLM"), function(object, vobj, coef){

	if( missing(vobj) ){
		stop("Must include response object as argument")
	}

	# Check that model fit and vobj have the same number of responses
	if( nrow(object) != nrow(vobj) ){
		stop("Model fit and data must have the same number of responses")
	}

	# Check that the responses have the same name
	if( ! identical(rownames(object), rownames(vobj)) ){		
		stop("Model fit and data must have the same responses")
	}

	if( object$method != "ls" ){
		stop("Only valid for models fit with lmFit()")
	}

	if( is( vobj, "EList") ){
		weights = t(vobj$weights)
	}else{
		weights = matrix(1, nrow(object$design), nrow(object))
	}

	# check that coef is valid
	if( ! missing(coef) ){
		i = match(coef, colnames(object$design))
		if( any(is.na(i)) ){
			txt = paste("Coefficients not valid:", paste(coef[is.na(i)], collapse=', '))
			stop(txt)
		}
	}

	# subsetting MArrayLM objects keep all residuals
	# so subset manually here
	features = rownames(coefficients(object))

	# use exact calculation for linear model
	eval_vcov( resids = t(residuals(object)[features,,drop=FALSE]), 
				X = object$design, 
				W = weights, 
				rdf = object$df.residual[1],
				coef = coef)
})



#' Co-variance matrix for \code{dream()} fit
#'
#' Define generic \code{vcov()} for result of \code{lmFit()} and \code{dream()}
# 
#' @param object \code{MArrayLM} object return by \code{lmFit()} or \code{dream()}
#' @param vobj \code{EList} object returned by \code{voom()}
#' @param coef name of coefficient to be extracted
#'
#' @importFrom stats coefficients
#' @export
setMethod('vcov', c("MArrayLM2"), function(object, vobj, coef){

	if( missing(vobj) ){
		stop("Must include response object as argument")
	}

	# Check that model fit and vobj have the same number of responses
	if( nrow(object) != nrow(vobj) ){
		stop("Model fit and data must have the same number of responses")
	}

	# Check that the responses have the same name
	if( ! identical(rownames(object), rownames(vobj)) ){		
		stop("Model fit and data must have the same responses")
	}

	if( is( vobj, "EList") ){
		weights = t(vobj$weights)
	}else{
		weights = matrix(1, nrow(object$design), nrow(object))
	}

	# check that coef is valid
	if( ! missing(coef) ){
		i = match(coef, colnames(object$cov.coefficients.list[[1]]))
		if( any(is.na(i)) ){
			txt = paste("Coefficients not valid:", paste(coef[is.na(i)], collapse=', '))
			stop(txt)
		}
	}

	# subsetting MArrayLM objects keep all residuals
	# so subset manually here
	features = rownames(coefficients(object))

	# use approximate calculation for linear mixed model
	eval_vcov_approx( 	resids = t(residuals(object)[features,,drop=FALSE]), 
						W = weights,
					  	ccl = object$cov.coefficients.list, 
					  	coef = coef)
})



# Evaluate variance-covariance matrix is multivariate regression
# 
# This method is exact for linear regression, even when each response has its own weight vector
# 
# @param resids matrix of residuals from regression
# @param X design matrix
# @param W matrix of precision weights
# @param rdf residual degrees of freedom
# @param coef name of coefficient to be extracted
#
eval_vcov = function(resids, X, W, rdf, coef){

	# With no weights:
	# kronecker(crossprod(res), solve(crossprod(X)), make.dimnames=TRUE) / rdf

	# pre-compute square root of W
	sqrtW = sqrt(W)

	# all pairs of responses
	Sigma = crossprod(resids * sqrtW) 

	# transpose X just once
	tX = t(X)

	# store dimensions of data
	k = ncol(X)
	m = ncol(resids)

	# matrix to store results
	Sigma_vcov = matrix(0, m*k, m*k)

	# Equivalent to kronecker product when weights are shared
	# outer loop
	for(i in seq(1, m)){

		# define positions in output matrix
		idx1 = seq((i-1)*k+1, i*k)

		# scale X by weights
		X_i = X*sqrtW[,i]

		# evaluate (X^TX)^-1 X^T for i 
		A = solve(crossprod(X_i), t(X_i))

		# inner loop
		for(j in seq(i, m)){

			idx2 = seq((j-1)*k+1, j*k)
			X_j = X*sqrtW[,j]
			B = solve(crossprod(X_j), t(X_j))

			value = (Sigma[i,j] / rdf) * tcrossprod(A,B)

			Sigma_vcov[idx1,idx2] = value 
			Sigma_vcov[idx2,idx1] = value 
		}
	}

	# assign names
	colnames(Sigma_vcov) = c(outer(colnames(X), colnames(resids), function(a,b) paste(b,a, sep=':')))
	rownames(Sigma_vcov) = colnames(Sigma_vcov)

	# select coefficients
	if( missing(coef) ){
		# all coefficients
		i = seq(1, ncol(X))
	}else{
		# subsect to selected coefs
		i = match(coef, colnames(X))
	}

	# names of coefficients to retain
	keep = c(outer(colnames(X)[i], colnames(resids), function(a,b) paste(b,a, sep=':')))

	# subset covariance matrix
	Sigma_vcov[keep,keep]
}


# Evaluate variance-covariance matrix is multivariate regression
# 
# This method is approximate since part of the calculations assume equal weights.  This is useful for the linear mixed model where the exact calculation is very challanging
# 
# @param resids matrix of residuals from regression
# @param W matrix of precision weights
# @param ccl list of vcov matrices for each response
# @param coef name of coefficient to be extracted
#
eval_vcov_approx = function(resids, W, ccl, coef){

	# store dimensions of data
	k = ncol(ccl[[1]])
	m = ncol(resids)

	# residual covariance
	scale_res = scale(resids*sqrt(W)) / sqrt(nrow(resids) - 1)
	Sigma = crossprod(scale_res)

	# matrix to store results
	Sigma_vcov = matrix(0, m*k, m*k)

	# Equivalent to kronecker product when weights are shared
	# outer loop
	for(i in seq(1, m)){

		# define positions in output matrix
		idx1 = seq((i-1)*k+1, i*k)

		chol_cov_i = chol(ccl[[i]])
	
		# inner loop
		for(j in seq(i, m)){

			idx2 = seq((j-1)*k+1, j*k)

			chol_cov_j = chol(ccl[[j]])
	
			value = Sigma[i,j] * crossprod(chol_cov_i, chol_cov_j)

			Sigma_vcov[idx1,idx2] = value 
			Sigma_vcov[idx2,idx1] = value 
		}
	}

	# assign names
	colnames(Sigma_vcov) = c(outer(colnames(ccl[[1]]), colnames(resids), function(a,b) paste(b,a, sep=':')))
	rownames(Sigma_vcov) = colnames(Sigma_vcov)

	# select coefficients
	if( missing(coef) ){
		# all coefficients
		i = seq(1, ncol(ccl[[1]]))
	}else{
		# subsect to selected coefs
		i = match(coef, colnames(ccl[[1]]))
	}

	# names of coefficients to retain
	keep = c(outer(colnames(ccl[[1]])[i], colnames(resids), function(a,b) paste(b,a, sep=':')))

	# subset covariance matrix
	Sigma_vcov[keep,keep]
}


