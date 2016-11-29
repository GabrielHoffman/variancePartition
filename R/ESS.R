# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai

#' Effective sample size
#' 
#' Compute effective sample size based on correlation structure in linear mixed model
#'
#' @param fit model fit from lmer()
#' @param method "full" uses the full correlation structure of the model. The "approximate" method makes the simplifying assumption that the study has a mean of m samples in each of k groups, and computes m based on the study design.  When the study design is evenly balanced (i.e. the assumption is met), this gives the same results as the "full" method.  
#' 
#' @return
#' effective sample size for each random effect in the model
#' 
#' @details
#'
#' Effective sample size calculations are based on:

#' Liu, G., and Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. Biometrics, 53(3), 937-47.
#'
#' "full" method: if V_x = var(Y;x) is the variance-covariance matrix of Y, the response, based on the covariate x, then the effective sample size corresponding to this covariate is \\Sigma_\{i,j\} (V_x^\{-1\})_\{i,j\}.  In R notation, this is: sum(solve(V_x)).  In practice, this can be evaluted as sum(w), where R %*% w == One and One is a column vector of 1's
#
# "fast" method: takes advantage of the fact that the eigen decompostion of a sparse, low rank, symmetric matrix. May be faster for large datasets.
#'
#' "approximate" method: Letting m be the mean number of samples per group, k be the number of groups, and rho be the intraclass correlation, the effective sample size is m*k / (1+rho*(m-1))
#'
#' Note that these values are equal when there are exactly m samples in each group.  If m is only an average then this an approximation.
#'
#' @examples
#' library(lme4)
#' data(varPartData)
#'
#' # Linear mixed model
#' fit <- lmer( geneExpr[1,] ~ (1|Individual) + (1|Tissue) + Age, info)
#'
#' # Effective sample size
#' ESS( fit )
#' 
#' @export
#' @docType methods
#' @rdname ESS-method
setGeneric("ESS", signature="fit",
  function(fit, method="full")
      standardGeneric("ESS")
)

#' @export
#' @rdname ESS-method
#' @aliases ESS,lmerMod-method
setMethod("ESS", "lmerMod",
	function( fit, method="full" ){

		if( !(method %in% c("full", "approximate")) ){
			stop(paste("method is not valid:", method))
		}

		# get correlation terms
		vp = calcVarPart( fit )

		n_eff = c()

		if( method %in% c('full') ){

			# get structure of study design
			sigG = get_SigmaG( fit )

			# number of samples 
			N = nrow(sigG$Sigma)

			# create vector of ones
			One = matrix(1, N, 1)

			ids = names(coef(fit))
			for( key in ids){
				i = which( key == ids)

				# guarantee that fraction is positive
				# by adding small value
				# the fast sum_of_ginv() fails if this value is exactly zero
				fraction = vp[[key]] + 1e-10

				if( method == "full" ){
					C = sigG$G[[i]] * fraction
					diag(C) = 1 # set diagonals to 1
					# n_eff[i] = sum(ginv(as.matrix(C)))

					# Following derivation at https://golem.ph.utexas.edu/category/2014/12/effective_sample_size.html
					# sum(solve(R)) is equivalent to sum(w) where R %*% w = 1
					# This uses sparse linear algebra and is extremely fast: ~10 x faster than 
					# sum_of_ginv with sparse pseudo inverse
					# November 29th, 2016
					n_eff[i] = sum(solve(C, One))

				}
				# Disabled in this version 
				# "fast" is exact and can be much faster for > 500 samples.
				# else{
				# 	# November 22, 2016
				# 	A = sigG$G[[i]] * fraction
				# 	value = 1 - A[1,1]
				# 	k = nlevels(fit@frame[[key]])
				# 	n_eff[i] = sum_of_ginv( A, value, k)
				# }
			}
			names(n_eff) = ids

		}else{

			ids = names(coef(fit))
			for( key in ids){
				i = which( key == ids)
				rho = vp[[key]]
				k = nlevels(fit@frame[[key]])
				m = nrow(fit@frame) / k 
				n_eff[i] = m*k / (1+rho*(m-1))
			}
			names(n_eff) = ids
		}

		# I think summing these values doesn't make sense
		# March 5, 2015
		#n_eff = c(n_eff, Total = sum(n_eff))

		return( n_eff )
	}
)



# # Compute sum( solve(A) + diag(value)) for low rank, sparse, symmetric A
# # sum( solve(A + diag(value, nrow(A)))) 
# sum_of_ginv = function(A, value, k){

# 	# # full rank
# 	# decompg = svd(A)
# 	# sum((decompg$u) %*% (((1/(decompg$d +value)))* t(decompg$v)))

# 	# # low rank, dense matrix
# 	# decompg = svd(A, nu=k, nv=k)
# 	# sum((decompg$u[,1:k]) %*% (((1/(decompg$d[1:k] +value)))* t(decompg$v[,1:k])))

# 	# low rank, sparse matrix
# 	decomp = eigs_sym(A, k)
# 	sum((decomp$vectors) %*% ((1/(decomp$values +value))* t(decomp$vectors)))
# }
