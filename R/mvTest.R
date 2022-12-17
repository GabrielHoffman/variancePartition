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
#' @param method statistical method used to perform multivariate test.  See details.  \code{'FE'} is a fixed effect test that models the covariance between coefficients.  \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#'  
#' @details See package \code{remaCor} for details about the \code{remaCor::RE2C()} test, and see \code{remaCor::LS()} for details about the fixed effect test.  When only 1 feature is selected, the original p-value is returned and the test statistic is set to \code{NA}.
#' 
#' For the \code{"RE2C"} test, the final test statistic is the sum of a test statistic for the mean effect (\code{stat.FE}) and heterogeneity across effects (\code{stat.het}).
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
#' @importFrom stats coefficients pchisq cov2cor
#' @export
mvTest = function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher")){

	method = match.arg(method)

	i = match(coef, colnames(coefficients(fit)))
	if( any(is.na(i)) ){
		txt = paste("Coefficients not valid:", paste(coef[is.na(i)], collapse=', '))
		stop(txt)
	}

	if( missing(features) ){
		features = rownames(fit)
	}else{
		# only test features that are characters
		feat.test = features[is.character(features)]
		i = match(feat.test, rownames(fit))
		if(any(is.na(i))){
			txt = paste("Features not found:", paste(feat.test[is.na(i)], collapse=', '))
			stop(txt)
		}
	}

	# extract coefficients from features
	tab = topTable(fit[features,], coef=coef, sort.by="none", number=Inf)
	beta = tab$logFC

	n_features = length(features)

	# extract covariance
	Sigma = vcov(fit[features,], vobj[features,], coef)

	if( method == "FE"){
		if( n_features == 1){
			# for one test, return estimated t-stat as stat
			df = data.frame(stat = tab$t, 
							pvalue = tab$P.Value, 
							n_features = 1,
							method = method)
		}else{
			res = LS(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

			df = data.frame(stat = res$beta / res$se,
							pvalue = res$p,
							n_features = n_features, 
							method = method)
		}

	}else if( method == "RE2C"){

		if( n_features == 1){
			# for one test, heterogeneity is zero
			df = data.frame(stat.FE = tab$t^2, 
							stat.het = 0,
							pvalue = tab$P.Value, 
							n_features = 1,
							method = method)
		}else{
			res = RE2C(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

			df = data.frame(stat.FE = res$stat1, 
							stat.het = res$stat2,
							pvalue = res$RE2Cp,
							n_features = n_features, 
							method = method)
		}

	}else if( method == "tstat"){

		if( n_features == 1){
			
			df = data.frame(stat = tab$t, 
							pvalue = tab$P.Value, 
							n_features = n_features,
							method = method)
		}else{
			Sigma_corr = cov2cor(Sigma)
			# tstat = crossprod(tab$t, solve(Sigma_corr, tab$t))

			# in case Sigma_corr is not invertable
			tstat = tryCatch( crossprod(tab$t, solve(Sigma_corr, tab$t)), error = function(e){
				warning("Covariance matrix is not invertable. Returning NA values.")
				NA
				})  

			pv = pchisq( tstat, length(tab$t), lower.tail=FALSE)

			df = data.frame(stat = tstat, 
							pvalue = pv, 
							n_features = n_features,
							method = method)
		}
	}else if( method == "sidak"){
		pv = 1 - (1 - min(tab$P.Value))^nrow(tab)

		df = data.frame(stat = NA,
						pvalue = pv, 
						n_features = n_features,
						method = method)
	}else if( method == "fisher"){
		stat = -2 * sum(log(tab$P.Value))
		k = nrow(tab)
		pv = pchisq(stat, 2*k, lower.tail=FALSE)

		df = data.frame(stat = stat,
						pvalue = pv, 
						n_features = n_features,
						method = method)
	}

	df
}

