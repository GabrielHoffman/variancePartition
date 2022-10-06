# Gabriel Hoffman 
# October 4, 2022
# 
# Perform multivariate tests on results from dream() using vcov()






#' Multivariate tests on results from \code{dream()}
#'
#' Evaluate multivariate tests on results from \code{dream()} using \code{vcov()} to compute the covariance between estimated regression coefficients across multiple responses.  A joint test to see if the coefficients are jointly different from zero is performed using meta-analysis methods that account for the covariance.
#' 
#' @param fit \code{MArrayLM} or \code{MArrayLM2} returned by \code{dream()}
#' @param vobj matrix or \code{EList} object returned by \code{voom()}
#' @param features indeces or names of features to perform multivariate test on 
#' @param coef name of coefficient or contrast to be tested
#' @param method statistical method used to perform multivariate test.  See details. \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'LS'} is a fixed effect test that models the covariance between coefficients.  \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#'  
#' @details See package \link{\code{remaCor}} for details about \code{remaCor::RE2C()} and \code{remaCor::LS()} methods.  When only 1 feature is selected, the original t-statistic and p-value are returned.
#' 
#' @examples
#' # library(variancePartition)
#' library(edgeR)
#' library(BiocParallel)
#' 
#' data(varPartDEdata)
#' 
#' # normalize RNA-seq counts
#' dge = DGEList(counts = countMatrix)
#' dge = calcNormFactors(dge)
#' 
#' # specify formula with random effect for Individual
#' form <- ~ Disease + (1|Individual) 
#' 
#' # compute observation weights
#' vobj = voomWithDreamWeights( dge[1:20,], form, metadata)
#' 
#' # fit dream model 
#' fit = dream( vobj, form, metadata)
#' fit = eBayes(fit)
#' 
#' # Multivariate test of features 1 and 2
#' mvTest(fit, vobj, 1:2, coef="Disease1")
#'
#' @importFrom remaCor RE2C LS 
#' @importFrom stats coefficients
#' @export
mvTest = function(fit, vobj, features, coef, method = c("RE2C", "LS", "tstat", "sidak", "fisher")){

	method = match.arg(method)

	i = match(coef, colnames(coefficients(fit)))
	if( any(is.na(i)) ){
		txt = paste("Coefficients not valid:", paste(coef[is.na(i)], collapse=', '))
		stop(txt)
	}

	if( missing(features) ){
		features = rownames(fit)
	}

	# check that all features are valid
	# This didn't handle index features, so comment out
	# idx = match(features, rownames(fit))
	# if( any(is.na(idx)) ){	

	# 	txt = paste("features not found:", paste(features[is.na(idx)], collapse=', '))
	# 	stop(txt)
	# }

	# extract coefficients from features
	tab = topTable(fit[features,], coef=coef, sort.by="none", number=Inf)
	beta = tab$logFC

	if( length(features) == 1){
		df = data.frame(stat = tab$t, 
						pvalue = tab$P.Value, 
						method = method)
		return(df)
	}

	# extract covariance
	Sigma = vcov(fit[features,], vobj[features,], coef)

	if( method == "LS"){
		res = LS(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

		df = data.frame(stat = res$beta / res$se,
						pvalue = res$p, 
						method = method)

	}else if( method == "RE2C"){
		res = RE2C(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

		df = data.frame(stat = res$stat1, 
						pvalue = res$RE2Cp, 
						method = method)

	}else if( method == "tstat"){
		tstat = crossprod(tab$t, solve(cov2cor(Sigma))) %*% tab$t 
		pv = pchisq( tstat, length(tab$t), lower.tail=FALSE)

		df = data.frame(stat = tstat, 
						pvalue = pv, 
						method = method)
	}else if( method == "sidak"){
		pv = 1 - (1 - min(tab$P.Value))^nrow(tab)

		df = data.frame(stat = NA,
						pvalue = pv, 
						method = method)
	}else if( method == "fisher"){
		stat = -2 * sum(log(tab$P.Value))
		k = nrow(tab)
		pv = pchisq(stat, 2*k, lower.tail=FALSE)

		df = data.frame(stat = stat,
						pvalue = pv, 
						method = method)
	}

	df
}

