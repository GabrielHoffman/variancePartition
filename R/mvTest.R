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
#' @param features a) indeces or names of features to perform multivariate test on, b) list of indeces or names.  If missing, perform joint test on all features.
#' @param coef name of coefficient or contrast to be tested
#' @param method statistical method used to perform multivariate test.  See details.  \code{'FE'} is a fixed effect test that models the covariance between coefficients.  \code{'RE2C'} is a random effect test of heterogeneity of the estimated coefficients that models the covariance between coefficients, and also incorporates a fixed effects test too. \code{'tstat'} combines the t-statistics and models the covariance between coefficients. \code{'sidak'} returns the smallest p-value and accounting for the number of tests. \code{'fisher'} combines the p-value using Fisher's method assuming independent tests.
#' @param shrink.cov shrink the covariance matrix between coefficients using the Schafer-Strimmer method 
#' @param progressbar if TRUE, show progress bar
#' @param ... other arugments
#'  
#' @details See package \code{remaCor} for details about the \code{remaCor::RE2C()} test, and see \code{remaCor::LS()} for details about the fixed effect test.  When only 1 feature is selected, the original p-value is returned and the test statistic is set to \code{NA}.
#' 
#' For the \code{"RE2C"} test, the final test statistic is the sum of a test statistic for the mean effect (\code{stat.FE}) and heterogeneity across effects (\code{stat.het}).  \code{mvTest()} returns 0 if \code{stat.het} is negative in extremely rare cases.
#' 
#' @return 
#' Returns a \code{data.frame} with the statistics from each test, the \code{pvalue} from the test, \code{n_features},  \code{method}, and \code{lambda} from the Schafer-Strimmer method to shrink the estimated covariance.  When \code{shrink.cov=FALSE}, \code{lambda = 0}.
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
#' # Test multiple sets of features
#' lst = list(a = 1:2, b=3:4)
#' mvTest(fit, vobj, lst, coef="Disease1")
#' @export
#' @docType methods
#' @rdname mvTest-method
setGeneric("mvTest", signature=c("fit", "vobj", 'features'),
  function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov=TRUE, progressbar=TRUE,...)
      standardGeneric("mvTest")
)


#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,integer-method
#' @importFrom remaCor RE2C LS 
#' @importFrom stats coefficients pchisq cov2cor
#' @importFrom corpcor estimate.lambda
#' @export
setMethod("mvTest", c('MArrayLM', "EList", "vector"),
function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov=TRUE, progressbar=TRUE,...){

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
	#-------------------
	# get Sigma directly
	# Sigma = vcov(fit[features,], vobj[features,], coef)

	# Instead, get sqrt of covariance so that Sigma is crossprod(P)
	# Then estimate shrinkage intensity
	# and shrink covariance
	# this is important when p approaches n
	# Note that the test below does not model uncertainy in lambda
	P = vcovSqrt(fit[features,], vobj[features,], coef, approx=TRUE)
	Sigma = crossprod(P)

	if( shrink.cov == FALSE) shrink.cov = "FALSE"
	lambda = switch( shrink.cov, 
						Schafer = estimate.lambda(P, verbose=FALSE),
						'FALSE' = 0)

	Sigma = (1-lambda) * Sigma + lambda * diag(diag(Sigma), ncol(Sigma))

	# Meta-analyis method
	#####################
	if( method == "FE"){
		if( n_features == 1){
			# for one test, return estimated t-stat as stat
			df = data.frame(stat = tab$t, 
							pvalue = tab$P.Value, 
							n_features = 1,
							lambda = lambda,
							method = method)
		}else{

			# if( is(fit, "MArrayLM") ){
			# 	# fixed effect model
			# 	nu = fit$df.residual
			# }else{
			# 	# mixed model
			# 	nu = fit$rdf
			# }	

			# if sample size is large enough
			# if( all(nu > 50) ){
				# Use asymptotic normal null distribution 
				res = LS(beta, sqrt(diag(Sigma)), cov2cor(Sigma))
			# }else{
			# 	# 
			# 	res = LS.empirical(beta, sqrt(diag(Sigma)), cov2cor(Sigma), nu)
			# }

			df = data.frame(stat = res$beta / res$se,
							pvalue = res$p,
							n_features = n_features,
							lambda = lambda, 
							method = method)
		}

	}else if( method == "RE2C"){

		if( n_features == 1){
			# for one test, heterogeneity is zero
			df = data.frame(stat.FE = tab$t^2, 
							stat.het = 0,
							pvalue = tab$P.Value, 
							n_features = 1,
							lambda = lambda,
							method = method)
		}else{
			res = RE2C(beta, sqrt(diag(Sigma)), cov2cor(Sigma))

			df = data.frame(stat.FE = res$stat1, 
							stat.het = max(0, res$stat2),
							pvalue = res$RE2Cp,
							n_features = n_features,
							lambda = lambda, 
							method = method)
		}

	}else if( method == "tstat"){

		if( n_features == 1){
			
			df = data.frame(stat = tab$t, 
							pvalue = tab$P.Value, 
							n_features = n_features,
							lambda = lambda,
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
						lambda = lambda,
							method = method)
		}
	}else if( method == "sidak"){
		pv = 1 - (1 - min(tab$P.Value))^nrow(tab)

		df = data.frame(stat = NA,
						pvalue = pv, 
						n_features = n_features,
						lambda = lambda,
						method = method)
	}else if( method == "fisher"){
		stat = -2 * sum(log(tab$P.Value))
		k = nrow(tab)
		pv = pchisq(stat, 2*k, lower.tail=FALSE)

		df = data.frame(stat = stat,
						pvalue = pv, 
						n_features = n_features,
						lambda = lambda,
						method = method)
	}

	df
})

# if no features are specified, test all features
#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,missing-method
#' @export
setMethod("mvTest", c('MArrayLM', "EList", "missing"),
function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov=TRUE, progressbar=TRUE,...){ 

	mvTest(fit, vobj, features=seq(nrow(vobj)), coef, method, shrink.cov=shrink.cov, progressbar=progressbar,... )
})



#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,EList,list-method
#' @import progress
#' @importFrom stats runif
#' @export
setMethod("mvTest", c('MArrayLM', "EList", "list"),
function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov=TRUE, progressbar=TRUE,...){ 

	if( is.null(names(features)) ){
		stop("features list must have non-null names(features)")
	}

	# set up progress bar
	pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta", total = length(features), width= 60, clear=FALSE)

	# features is list
	res = lapply( names(features), function(id){
		if( progressbar & runif(1) < .05){
			ratio = match(id, names(features)) / length(features)
			pb$update(ratio = ratio )
		}
		res = mvTest(fit, vobj, features[[id]], coef, method, shrink.cov=shrink.cov)
		data.frame(ID = id, res)
		})
	res = do.call(rbind, res)

	if( progressbar & ! pb$finished ){
		pb$update(1.0)
	}
	pb$terminate()

	res
})

#' @rdname mvTest-method
#' @aliases mvTest,MArrayLM,matrix-method
#' @export
setMethod("mvTest", c('MArrayLM', "matrix"),
function(fit, vobj, features, coef, method = c("FE", "RE2C", "tstat", "sidak", "fisher"), shrink.cov=TRUE, progressbar=TRUE,...){
		
	method = match.arg(method)	

	# create weights
	W = matrix(1, nrow(vobj), ncol(vobj))

	vobj = new("EList", list(E = vobj, weights = W)) 

	mvTest(fit, vobj, features, coef, method, shrink.cov=shrink.cov, progressbar=progressbar,...  )
})


