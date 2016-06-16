# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai

#' Effective sample size
#' 
#' Compute effective sample size based on correlation structure in linear mixed model
#'
#' @param fit model fit from lmer()
#' @param method "full" uses the full correlation structure of the model. The "approximate" method makes the simplifying assuption that the study has a mean of m samples in each of k groups, and computes m based on the study design.  When the study design is evenly balanced (i.e. the assumption is met), this gives the same results as the "full" method.  
#' 
#' @return
#' effective sample size for each random effect in the model
#' 
#' @details
#'
#' Effective sample size calculations are based on:

#' Liu, G., and Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. Biometrics, 53(3), 937-47.
#'
#' "full" method: if V_x = var(Y;x) is the variance-covariance matrix of Y, the respose, based on the covariate x, then the effective sample size corresponding to this covariate is \\Sigma_\{i,j\} (V_x^\{-1\})_\{i,j\}.  In R notation, this is: sum(solve(V_x)).
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

		if( method == 'full'){

			# get structure of study design
			sigG = get_SigmaG( fit )

			ids = names(coef(fit))
			for( key in ids){
				i = which( key == ids)
				C = as.matrix(sigG$G[[i]]) * vp[[key]]
				diag(C) = 1 # set diagonals to 1
				n_eff[i] = sum(ginv(as.matrix(C)))
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

