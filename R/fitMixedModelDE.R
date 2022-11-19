
#' Class MArrayLM2
#'
#' Class \code{MArrayLM2} 
#'
#' @name MArrayLM2-class
#' @rdname MArrayLM2-class
#' @exportClass MArrayLM2
setClass("MArrayLM2",
#  Linear model fit
representation("MArrayLM")#, varComp="data.frame", sigGStruct='list')
)

setIs("MArrayLM2","LargeDataObject")
setAs(from='MArrayLM', to='MArrayLM2', function(from){
	structure(from, class="MArrayLM2")
	})


# define S3 version of these functions

#' Residuals for result of dream
#' 
#' Residuals for result of dream
#' 
#' @param object See \code{?stats::residuals}
#' @param ... See \code{?stats::residuals}
#'
#' @rawNamespace S3method("residuals", MArrayLM2)
#' @export   
residuals.MArrayLM2 = function( object, ...){
	if( is.null(object$residuals) ){
		stop( "Residuals were not computed, must run:\n dream(...,computeResiduals=TRUE)")
	}
	if( nargs() > 1 ){#& is.null(suppressWarnings) ){
		warning("\n Second argument is ignored here,\n but can be passed for compatability with limma.\n Results are the same either way")
	}
	object$residuals
}


# #' @importFrom limma residuals.MArrayLM
# # @export  
# residuals.MArrayLM = function( object, ...){

# 	if( is.null(object$residuals) ){
# 		# use residuals computed by limma
# 		res = limma::residuals.MArrayLM( object, ...)
# 	}else{
# 		# use precomputed residuals
# 		res = object$residuals
# 	}
# 	res
# }


# S4 methpds

#' residuals for MArrayLM
#'
#' residuals for MArrayLM
#'
#' @param object MArrayLM object from dream
#' @param ... other arguments, currently ignored
#'
#' @return results of residuals
#' @export
setMethod("residuals", "MArrayLM",
	function( object, ...){
		if( nargs() == 1 ){
			result = object$residuals
		}else{
			result = residuals.MArrayLM(object,...)
		}
		result
	})


#' residuals for MArrayLM2
#'
#' residuals for MArrayLM2
#'
#' @param object MArrayLM2 object from dream
#' @param ... other arguments, currently ignored
#'
#' @return results of residuals
#' @export
setMethod("residuals", "MArrayLM2",
	function( object, ...){
		residuals.MArrayLM2(object,...)
		})




# Evaluate contrasts for linear mixed model
# 
# Evaluate contrasts for linear mixed model
#
# @param fit model fit
# @param L contrast matrix
# @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough approximation for most datasets.
# 
# @return
# df, sigma, beta, SE of model
#
# @details If the Kenward-Roger covariance matrix is not positive definite, the Satterthwaite method is used
#
# @export
# @docType methods
# @rdname eval_lmm-method
#' @import foreach
#' @importFrom lme4 fixef 
#' @importFrom stats sigma
.eval_lmm = function( fit, L, ddf ){

	j = 1
	# evaluate each contrast
	# cons = lmerTest::contest(fit, L, ddf=ddf)
	cons = foreach( j = 1:ncol(L), .combine=rbind) %do% {
		lmerTest::contest(fit, L[,j], ddf=ddf)
	}

	df = as.numeric(cons[,'DenDF'])

	if(ddf == "Kenward-Roger"){
		# KR
		V = pbkrtest::vcovAdj.lmerMod(fit, 0)

		# if matrix is not PSD
		if( min(diag(as.matrix(V))) < 0){			
			warning("The adjusted Kenward-Roger covariance matrix is not positive definite.\nUsing Satterthwaite approximation instead")

			# Satterthwaite
			V = vcov(fit)
		}
		# df = pbkrtest::get_Lb_ddf(fit, L)
	}else{
		# Satterthwaite
		V = vcov(fit)
		# df = as.numeric(contest(fit, L, ddf="Sat")['DenDF'])
	}

	# sigma = attr(lme4::VarCorr(fit), "sc")

	# get contrasts
	# beta = as.matrix(sum(L * fixef(fit)), ncol=1)
	# colnames(beta) = "logFC"

	beta = foreach( j = 1:ncol(L), .combine=rbind) %do% {
		as.matrix(sum(L[,j] * fixef(fit)), ncol=1)
	}
	colnames(beta) = "logFC"
	rownames(beta) = colnames(L)

	# SE = as.matrix(sqrt(sum(L * (V %*% L))), ncol=1)		
	# colnames(SE) = "logFC"
	SE = foreach( j = 1:ncol(L), .combine=rbind) %do% {
		as.matrix(sqrt(sum(L[,j] * (V %*% L[,j]))), ncol=1)
	}
	colnames(SE) = "logFC"
	rownames(SE) = colnames(L)

	# pValue = 2*pt(as.numeric(abs(beta / SE)), df, lower.tail=FALSE)
	pValue = as.numeric(cons[,'Pr(>F)'])

	list(	cons 	= cons,
			df		= df,
			sigma	= sigma(fit),
			beta	= beta,
			SE		= SE,
			pValue	= pValue,
			vcov 	= V )
}


#' @importFrom methods is
.checkNA = function(exprObj){

	if( is(exprObj, "sparseMatrix") || is( exprObj, "matrix") ){
		countNA = sum(!is.finite(exprObj))
	}else{		
		# is.finite is not defined for data.frames, so convert to matrix first
		countNA = sum(!is.finite(as.matrix(exprObj)))

		# check if values are NA
		# countNA = sum(!is.finite(exprObj)) # sum(is.nan(exprObj))
	}

	if( countNA > 0 ){
		stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
	}

	# check if all genes have variance
	rv = apply( exprObj, 1, var)
	if( any( rv == 0) ){
		idx = which(rv == 0)
		stop(paste("Response variable", idx[1], 'has a variance of 0'))
	}		
}

#' Differential expression with linear mixed model
#' 
#' Fit linear mixed model for differential expression and preform hypothesis test on fixed effects as specified in the contrast matrix \code{L}
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.
#' @param data data.frame with columns corresponding to formula 
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough approximation for most datasets. "adaptive" (Default) uses KR for <= 10 samples.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param weightsMatrix matrix the same dimension as \code{exprObj} with observation-level weights from \code{voom()}.  Used only if useWeights is TRUE 
#' @param control control settings for \code{lmer()}
#' @param suppressWarnings if TRUE, do not stop because of warnings or errors in model fit
#' @param quiet suppress message, default FALSE
#' @param BPPARAM parameters for parallel evaluation
#' @param computeResiduals if TRUE, compute residuals and extract with \code{residuals(fit)}.  Setting to FALSE saves memory
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is TRUE.  See Details.
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#' 
#' @return 
#' MArrayLM2 object (just like MArrayLM from limma), and the directly estimated p-value (without eBayes)
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression (Hoffman and Roussos, 2021).  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is 
#'
#' \code{fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)}
#'
#' \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run code in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lmer.
#'
#' Hypothesis tests and degrees of freedom are producted by \code{lmerTest} and \code{pbkrtest} pacakges
#'
#' While \code{REML=TRUE} is required by \code{lmerTest} when ddf='Kenward-Roger', ddf='Satterthwaite' can be used with \code{REML} as \code{TRUE} or \code{FALSE}.  Since the Kenward-Roger method gave the best power with an accurate control of false positive rate in our simulations, and since the Satterthwaite method with REML=TRUE gives p-values that are slightly closer to the Kenward-Roger p-values, \code{REML=TRUE} is the default.  See Vignette "3) Theory and practice of random effects and REML"
#'
#' @references{
#'   \insertRef{hoffman2021dream}{variancePartition}
#' }
#' @examples
#' # library(variancePartition)
#' library(BiocParallel)
#'
#' # load simulated data:
#' # geneExpr: matrix of *normalized* gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' 
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' # NOTE: dream() runs on *normalized* data
#' fit = dream( geneExpr[1:10,], form, info)
#' fit = eBayes(fit)
#'
#' # view top genes
#' topTable( fit, coef="Batch2", number=3 )
#'
#' # get contrast matrix testing if the coefficient for Batch3 is 
#' # different from coefficient for Batch2
#' # The variable of interest must be a fixed effect
#' L = makeContrastsDream(form, info, contrasts=c("Batch3 - Batch2"))
#' 
#' # plot contrasts
#' plotContrasts( L )
#' 
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit2 = dream( geneExpr[1:10,], form, info, L)
#' fit = eBayes(fit)
#' 
#' # view top genes
#' topTable( fit2, coef="Batch3 - Batch2", number=3)
#' 
#' # Parallel processing using multiple cores with reduced memory usage
#' param = SnowParam(4, "SOCK", progressbar=TRUE)
#' fit3 = dream( geneExpr[1:10,], form, info, L, BPPARAM = param)
#' fit = eBayes(fit)
#'
#' # Fit fixed effect model for each gene
#' # Use lmFit in the backend
#' form <- ~ Batch 
#' fit4 = dream( geneExpr[1:10,], form, info, L)
#' fit4 = eBayes( fit4 )
#' 
#' # view top genes
#' topTable( fit4, coef="Batch3 - Batch2", number=3 )
#'
#' # Compute residuals using dream
#' residuals(fit4)[1:4, 1:4]
#' 
#' @export
# @docType methods
#' @rdname dream-method
#' @importFrom pbkrtest get_SigmaG
#' @importFrom BiocParallel bpiterate bpok SerialParam
#' @importFrom lme4 VarCorr 
#' @importFrom stats hatvalues as.formula
#' @importFrom foreach foreach
#' @importFrom methods as
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @import doParallel
#'
dream <- function( exprObj, formula, data, L, ddf = c("adaptive", "Satterthwaite", "Kenward-Roger"), useWeights=TRUE, weightsMatrix=NULL, control = vpcontrol,suppressWarnings=FALSE, quiet=FALSE, BPPARAM=SerialParam(), computeResiduals=TRUE, REML=TRUE, ...){

	exprObjInit = exprObj

	exprObjMat = as.matrix( exprObj )
	formula = as.formula( formula )

	# convert to data.frame
	data = as.data.frame(data)

	# only retain columns used in the formula
	# This reduces overhead for parallel processing with large datasets
	data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]
	data = droplevels(data)

	ddf = match.arg(ddf)
	colinearityCutoff = 0.999

	if( ! is.data.frame(data) ){
		stop("data argument must be a data.frame")
	}

	# check dimensions of reponse and covariates
	if( ncol(exprObj) != nrow(data) ){		
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
	}

	# check if variables in formula has NA's
	hasNA = hasMissingData(formula, data)

	if( any(hasNA) ){
		warning(paste("\nVariables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\nSamples with missing data will be dropped.\n"), immediate.=TRUE, call.=FALSE)

		# drop samples with missing data in formula variables
	    idx = sapply(all.vars(formula), function(v) {
	        which(is.na(data[[v]]))
	    })
	    idx = unique(unlist(idx))
	    
	    data = droplevels(data[-idx,,drop=FALSE])
	    exprObj = exprObj[,-idx,drop=FALSE]
	    exprObjMat = as.matrix( exprObj )
	}

	# assign weightsMatrix from exprObj
	if( is( exprObj, "EList") ){
		if( useWeights ){
			weightsMatrix = exprObj$weights
		}
		.checkNA( exprObj$E )
	}else{
		.checkNA( exprObj )
	}

	# "adaptive" (Default) uses KR for <= 10 samples.
	if( ddf == "adaptive" ){
		ddf = ifelse( ncol(exprObj) <= 10, "Kenward-Roger", 'Satterthwaite')
	}

	if( !(ddf %in% c("Kenward-Roger", 'Satterthwaite')) ){
		stop("Specify ddf correctly")
	}

	if( ddf == "Kenward-Roger" & ! REML ){
		stop("Kenward-Roger must be used with REML")
	}

	# if weightsMatrix is not specified, set useWeights to FALSE
	if( useWeights && is.null(weightsMatrix) ){
		# warning("useWeights was ignored: no weightsMatrix was specified")
		useWeights = FALSE
	}

	# if useWeights, and (weights and expression are the same size)
	if( useWeights && !identical( dim(exprObj), dim(weightsMatrix)) ){
		 stop( "exprObj and weightsMatrix must be the same dimensions" )
	}

	if( ! useWeights ){
		weightsMatrix = NULL	
	}
	
	# If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
	if( ! identical(colnames(exprObj), rownames(data)) ){
		 warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
	}

	if( .isDisconnected() ){
		stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
	}

	# Contrasts
	###########

	univariateContrasts = FALSE
	if( missing(L) || is.null(L) ){
		# all univariate contrasts
		L = .getAllUniContrasts( formula, data)
		univariateContrasts = TRUE
	}else{
		# format contrasts 
		if( is(L, "numeric") ){
			L = as.matrix(L, ncol=1)
		}else if( is(L, 'data.frame') ){
			L = as.matrix(L)
		}
		if( is.null(colnames(L)) ){
			colnames(L) = paste0('L', seq_len(ncol(L)))
		}

		# check columns that only have a single 1
		tst = apply(L, 2, function(x){
			length(x[x!=0]) == 1
			})
		if( any(tst) ){
			warning("Contrasts with only a single non-zero term are already evaluated by default.")
		}
		# remove univariate contrasts
		# L = L[,!tst,drop=FALSE]

		# add all univariate contrasts
		Luni = .getAllUniContrasts( formula, data)
		L = cbind(L, Luni)
	}

	if( ncol(L) == 0){
		stop( "Must include fixed effect in the model for hypothesis testing")
	}

	# check rownames of contrasts
	if( length(unique(colnames(L))) != ncol(L) ){
		stop(paste("Contrast names must be unique: ", paste(colnames(L), collapse=', ')))
	}

	# Evaluate model on each gene
	#############################
	
	if( ! .isMixedModelFormula( formula ) ){
		if( !quiet ){
			message("Fixed effect model, using limma directly...")
			message("User can apply eBayes() afterwards...")
		}

		design = model.matrix( formula, data)

		if( useWeights ){			
			# least squares fit
			ret = lmFit( exprObj, design, weights=weightsMatrix )
		}else{
			if( "weights" %in% names(exprObj)){
				exprObj$weights = NULL
			}
			ret = lmFit( exprObj, design )
		}

		if( computeResiduals ){
			# Evaluate residuals here, since can't be evaluate after contrasts.fit() is run
			# add this to the standard MArrayLM object for use by custom residuals function
			ret$residuals = residuals( ret, exprObj )
		}

		# apply contrasts
		if( ! univariateContrasts ){
			ret = contrasts.fit( ret, L)
		}

		# compute hat values
		ret$hatvalues = hatvalues(ret, exprObj)

	}else{

		# add response (i.e. exprObj[j,]) to formula
		# Use term 'responsePlaceholder' to store the value of reponse j (exprObj[j,])
		# This is not an ideal way to structure the evaluation of response j.
		# 	The formula is evaluated a different scope, (i.e. within lmer()), and there is no need to pass the
		# 	entire exprObj object into that scope.  With lexical scope, I found it was possible that  
		# 	the value of exprObj[j,] could be different when evaluated in the lower scope
		# This placeholder term addresses that issue
		form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

		# fit first model to initialize other model fits
		# this make the other models converge faster
		responsePlaceholder = nextElem(exprIter(exprObjMat, weightsMatrix, useWeights, scale=FALSE))

		timeStart = proc.time()
		fitInit <- lmerTest::lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )

		# extract covariance matrices  
		sigGStruct = get_SigmaG( fitInit )$G
		
		# check L
		if( ! identical(rownames(L), names(fixef(fitInit))) ){
			stop("Names of entries in L must match fixed effects")
		}
		
		# extract statistics from model
		mod = .eval_lmm( fitInit, L, ddf)
		timediff = proc.time() - timeStart

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=!suppressWarnings, dream=TRUE, colinearityCutoff=colinearityCutoff )

		a = names(fixef(fitInit))
		b = rownames(L)

		if( ! identical( a,b) ){
			stop("Terms in contrast matrix L do not match model:\n  Model: ", paste(a, collapse=',' ), "\n  L: ", paste(b, collapse=',' ), "\nNote thhat order must be the same")
		}

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=responsePlaceholder$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta", total = nrow(exprObj), width= 60, clear=FALSE)

		timeStart = proc.time()

		# Define function for parallel evaluation
		.eval_models = function(responsePlaceholder, data2, form, REML, theta, control, na.action=stats::na.exclude,...){	
			# modify data2 for this gene
			data2$expr = responsePlaceholder$E
 
			# fit linear mixed model
			suppressWarnings({
				fit <- lmerTest::lmer( eval(parse(text=form)), data=data2, REML=REML,..., weights=responsePlaceholder$weights, control=control,na.action=na.action)
				})

			# extract statistics from model
			mod = .eval_lmm( fit, L, ddf)

			res = NULL
			if( computeResiduals ){
				 res = residuals(fit)
			}

			ret = list(	coefficients 	= mod$beta, 
						design 			= fit@pp$X,
						# approximate df of residuals 
						rdf 			= rdf.merMod(fit), 
						# df of test statistics
						df.residual 	= mod$df, 
						Amean 			= mean(fit@frame[,1]), 
						method 			= 'lmer',
						sigma 			= mod$sigma,
						stdev.unscaled 	= mod$SE/mod$sigma,
						pValue 			= mod$pValue,
						residuals 		= res ) 

			# get variance terms for random effects
			varComp <- lapply(lme4::VarCorr(fit), function(fit) attr(fit, "stddev")^2)
			varComp[['resid']] = attr(lme4::VarCorr(fit), 'sc')^2

			if( univariateContrasts ){
				V = mod$vcov
			}else{
				# do expensive evaluation only when L is defined by the user
				V = crossprod(chol(mod$vcov) %*% L)
			}

			# compute hatvalues
			h = hatvalues(fit)

			list( 	ret 	= new("MArrayLM", ret),
					varComp = varComp,
					# effective degrees of freedom as sum of diagonals of hat matrix
					edf		= sum(h), 
					hatvalues = h,
					vcov 	= V)
		}

		.eval_master = function( obj, data2, form, REML, theta, control, na.action=stats::na.exclude,... ){

			# use only 1 OpenMP thread for linear algebra
			omp_set_num_threads(1)

			lapply(seq_len(nrow(obj$E)), function(j){
				.eval_models( list(E=obj$E[j,], weights=obj$weights[j,]), data2, form, REML, theta, control, na.action,...)
			})
		}

		# Evaluate function
		####################

		it = iterBatch(exprObjMat, weightsMatrix, useWeights, scale=FALSE, n_chunks = 100, BPPARAM = BPPARAM)

		if( !quiet ) message(paste0("Dividing work into ",attr(it, "n_chunks")," chunks..."))

		resList <- bpiterate( it, .eval_master,
		                      data2=data2, form=form, REML=REML, theta=fitInit@theta, control=control,...,
		                      BPPARAM=BPPARAM)

		# if there is an error in evaluating fxn (usually in parallel backend)
		if( !bpok(list(resList)) ){
      		stop("Error evaluating fxn:\n\n", resList)
		}
		# It can also return a list of errors, or a list where only some elements are errors
		if( !all(bpok(resList)) ){
		      first_error <- resList[[which(!bpok(resList))[1]]]
		      stop("Error evaluating fxn:\n\n", first_error)
		}

		# If no errors, then it's safe to concatenate all the results together.
		resList <- do.call(c, resList)

		names(resList) = seq_len(length(resList))

		if( !quiet ) message("\nTotal:", paste(format((proc.time() - timeStart)[3], digits=1, scientific = FALSE), "s"))

		x = 1
		# extract results
		coefficients = foreach( x = resList, .combine=cbind ) %do% {x$ret$coefficients}
		df.residual = foreach( x = resList, .combine=cbind ) %do% {	x$ret$df.residual}
		rdf = foreach( x = resList, .combine=cbind ) %do% {	x$ret$rdf}
		pValue = foreach( x = resList, .combine=cbind ) %do% { x$ret$pValue }
		stdev.unscaled  = foreach( x = resList, .combine=cbind ) %do% {x$ret$stdev.unscaled} 
		if( computeResiduals ){
			residuals = foreach( x = resList, .combine=cbind ) %do% { x$ret$residuals }
		}
		hat = foreach( x = resList, .combine=cbind ) %do% {x$hatvalues}
		if( ! is.matrix(hat) ){ hat = matrix(hat,ncol=1)}
		colnames(hat) = rownames(exprObj)

		# transpose
		coefficients = t( coefficients )
		df.residual = t( df.residual )
		pValue = t( pValue )
		stdev.unscaled = t( stdev.unscaled )
		hat = t( hat )
		
		if( computeResiduals ){
			residuals = t(residuals)
			rownames(residuals) = rownames(exprObj)
		}

		colnames(coefficients) = colnames(L)
		rownames(coefficients) = rownames(exprObj)

		design = resList[[1]]$ret$design

		colnames(df.residual) = colnames(L)
		rownames(df.residual) = rownames(exprObj)

		colnames(pValue) = colnames(L)
		rownames(pValue) = rownames(exprObj)

		Amean = sapply( resList, function(x) x$ret$Amean) 
		names(Amean ) = rownames(exprObj)
		method = "lmer"
		sigma = sapply( resList, function(x) x$ret$sigma)
		names(sigma) = rownames(exprObj)

		varComp = lapply(resList, function(x){
			x = unlist(x$varComp)
			names(x) = gsub("\\.\\(Intercept\\)", '', names(x))
			as.data.frame(t(x))
		})
		varComp = do.call("rbind", varComp)
		rownames(varComp) = rownames(coefficients)

		edf = sapply(resList, function(x) x$edf)

		colnames(stdev.unscaled) = colnames(L)
		rownames(stdev.unscaled) = rownames(exprObj)

		ret = list( coefficients 	= coefficients,
		 			design 			= design, 
		 			rdf 			= c(rdf), 
		 			df.residual		= df.residual, 
		 			hatvalues 		= hat,
		 			Amean 			= Amean, 
		 			method 			= method, 
		 			sigma 			= sigma, 
		 			contrasts 		= L,
		 			stdev.unscaled 	= stdev.unscaled)

		if( 'genes' %in% names(exprObjInit) ){
			ret$genes = exprObjInit$genes
		}

		ret = new("MArrayLM", ret)
		# ret$pValue = pValue

		# set covariance between covariates
		# C = solve(crossprod(ret$design))
		# C = chol2inv(chol(crossprod(ret$design)))
		V = chol2inv(qr(ret$design)$qr)
		rownames(V) = colnames(ret$design)
		colnames(V) = colnames(ret$design)

		# remove intercept term if it exists
		# fnd = colnames(C) == "(Intercept)" 
		# if( any(fnd) ){
		# 	C = C[!fnd,!fnd,drop=FALSE]
		# }

		if( ! univariateContrasts ){
			# do expensive evaluation only when L is defined by the user
			V = crossprod(chol(V) %*% L)
		}

		# covariance used by limma.  
		# Assumes covariance is the same for all genes
		ret$cov.coefficients = V

		# allows covariance to differ for each gene based on variance components
		ret$cov.coefficients.list = lapply(resList, function(x) as.matrix(x$vcov))

		if( computeResiduals ){
			ret$residuals = residuals
		}

		ret = as(ret, "MArrayLM2")
		# add additional information for pinnacle
		# ret = new("MArrayLM2", ret, varComp, sigGStruct)
		attr(ret, "varComp") = varComp
		attr(ret, "sigGStruct") = sigGStruct
		attr(ret, "edf") = edf

		# compute standard values for t/F/p without running eBayes
		# eBayes can be run afterwards, if wanted
		ret = .standard_transform( ret )
	}

	# return formula, data
	ret$formula = formula
	ret$data = data

	ret
}


#' Compute standard post-processing values
#' 
#' These values are typically computed by eBayes
#' 
#' @param fit result of dream (MArrayLM2)
#' @param sigma vector of standard errors used to compute t-statistic. Can be maximum likelihood estimates, or posterior means
#'
#' @return MArrayLM2 object with values computed
#' 
#' @importFrom stats pchisq pf
#' @keywords internal
.standard_transform = function(fit, sigma = fit$sigma){

	# If fit$df.prior is not defined, set df.prior to zero
	if( ! is.null(fit$df.prior) ){
		fit$df.total = fit$df.residual + fit$df.prior
	}else{		
		fit$df.total = fit$df.residual
	}

	# t-test
	out = fit
	out$t <- fit$coefficients / fit$stdev.unscaled / sigma
	out$p.value <- 2*pt(abs(out$t), df=fit$df.total, lower.tail=FALSE )
	# out$p.value.loge <- log(2) + pt(abs(out$t), df=fit$df.total, lower.tail=FALSE, log.p=TRUE )

	# F-test
	if(!is.null(out$design) && is.fullrank(out$design)) {

		# only evaluate F-stat on real coefficients, not contrasts
		realcoef = colnames(out)[colnames(out) %in% colnames(out$design)]
		realcoef = realcoef[realcoef!="(Intercept)"]

		if( is.null(realcoef) || (length(realcoef) == 0) ){

			# this happends when only the intercept term is included
			warning("No testable fixed effects were included in the model.\n  Running topTable() will fail.")
		}else{
			df = rowMeans(out[,realcoef]$df.total)

			F.stat <- classifyTestsF(out[,realcoef], df=df, fstat.only=TRUE)
			out$F <- as.vector(F.stat)
			df1 <- attr(F.stat,"df1")
			df2 <- attr(F.stat,"df2")
			if(df2[1] > 1e6){ # Work around bug in R 2.1
				out$F.p.value <- pchisq(df1*out$F,df1,lower.tail=FALSE)
			}else{
				out$F.p.value <- pf(out$F,df1,df2,lower.tail=FALSE)
			}
		}
	}

	# if fit$df.prior does not exist, then remove the df.total term
	if( is.null(fit$df.prior) ){
		out$df.total = NULL
	}

	out
}




#' Subseting for MArrayLM2
#'
#' Enable subsetting on MArrayLM2 object. Same as for MArrayLM, but apply column subsetting to df.residual and cov.coefficients.list
#'
#' @param object MArrayLM2
#' @param i row
#' @param j col
#'
#' @name [.MArrayLM2
#' @return subset
#' @export
# @S3method [ MArrayLM2
#' @rawNamespace S3method("[", MArrayLM2)
#' @importFrom stats p.adjust
#' @rdname subset.MArrayLM2-method
#' @aliases subset.MArrayLM2,MArrayLM2-method
#' @keywords internal
assign("[.MArrayLM2",
	function(object, i, j){
	if(nargs() != 3){
		stop("Two subscripts required",call.=FALSE)
	}

	# apply standard MArrayLM subsetting
	obj = as(object, 'MArrayLM')

	if(!missing(j)){
		obj = obj[,j]
	}
	if(!missing(i)){
		obj = obj[i,]
	}

	# custom code to deal with df.total, df.residual and rdf
	if( is.null(ncol(object$df.total)) ){
		if(!missing(i)){
			obj$df.total = object$df.total[i]
		}else{
			obj$df.total = object$df.total
		}
	}else{			
		tmp = object$df.total
		if(!missing(i)){
			tmp = object$df.total[i,,drop=FALSE]
		}		
		if(!missing(j)){
			tmp = tmp[,j,drop=FALSE]
		}	
		obj$df.total = tmp
	}

	if( is.null(ncol(object$df.residual)) ){
		if(!missing(i)){
			obj$df.residual = object$df.residual[i]
		}else{			
			obj$df.residual = object$df.residual
		}
	}else{	
		tmp = object$df.residual
		if(!missing(i)){
			tmp = object$df.residual[i,,drop=FALSE]
		}		
		if(!missing(j)){
			tmp = tmp[,j,drop=FALSE]
		}	
		obj$df.residual = tmp
	}

	if( ! is.null(object$rdf) ){
		if(!missing(i)){
			obj$rdf = object$rdf[i]
		}else{
			obj$rdf = object$rdf
		}
	} 

	# obj$pValue = object$pValue[i,j]
	obj$s2.prior = object$s2.prior
	obj$df.prior = object$df.prior

	obj = as(obj, "MArrayLM2")

	#  copy gene-specific covariance, if it exists
	if( ! is.null(object$cov.coefficients.list) ){
		if(!missing(i)){
			if( is.numeric(i) ){
				# extract by index
				obj$cov.coefficients.list = object$cov.coefficients.list[i]
			}else{				
				# extract by matching feature name
				idx = match(i, rownames(object))
				obj$cov.coefficients.list = object$cov.coefficients.list[idx]
			}
		}else{
			obj$cov.coefficients.list = object$cov.coefficients.list
		}
		# name cov.coefficients.list using names of the whole object
		names(obj$cov.coefficients.list) = rownames(obj)
	}	

	if( is.null(obj$df.total)){
		obj$df.total = rowMeans(obj$df.residual)
	}

	# the F-statistic and p-value are evaluated when subsetting is applied
	# so need to apply df2 here
	# If columns have been subsetted, need to re-generate F
	if(!is.null(obj[["F"]]) && !missing(j)) {

		F.stat <- classifyTestsF(obj,df=obj$df.total,fstat.only=TRUE)
		obj$F <- as.vector(F.stat)
		df1 <- attr(F.stat,"df1")
		df2 <- attr(F.stat,"df2")
		if (df2[1] > 1e6){ 
			obj$F.p.value <- pchisq(df1*obj$F,df1,lower.tail=FALSE)
		}else{
			obj$F.p.value <- pf(obj$F,df1,df2,lower.tail=FALSE)
		}
	}
	obj	
})


setGeneric("eBayes", function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1) ){

	eBayes(fit, proportion, stdev.coef.lim, trend, robust, winsor.tail.p  )
})


#' eBayes for MArrayLM2
#'
#' eBayes for result of linear mixed model for with \code{dream()} using residual degrees of freedom approximated with \code{rdf.merMod()}
#'
#' @param fit fit
#' @param proportion proportion
#' @param stdev.coef.lim stdev.coef.lim
#' @param trend trend
#' @param robust robust
#' @param winsor.tail.p winsor.tail.p 
#'
#' @return results of eBayes using approximated residual degrees of freedom
#'
#' @export
#' @rdname eBayes-method
#' @aliases eBayes,MArrayLM2-method
#' @importFrom limma eBayes
#' @seealso dream rdf.merMod
setMethod("eBayes", "MArrayLM2",
function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)){

	# limma::eBayes() uses df.residual as the residual degrees of freedom,
	# 	while dream() uses rdf.
	# For linear models these values are always equal,
	# but for linear mixed models, rdf must be computed separately 
	# save df for test statistics
	df.test = fit$df.residual
	fit$df.residual = fit$rdf

	# Use limma::eBayes() but with new df.residual values
	fit_eb = limma::eBayes(	fit 			= fit, 
							proportion 		= proportion, 
							stdev.coef.lim 	= stdev.coef.lim, 
							trend 			= trend, 
							robust 			= robust, 
							winsor.tail.p 	= winsor.tail.p )

	# re-set to the df.residual of the test statistics
	fit_eb$df.residual = df.test
	fit_eb$rdf = fit$rdf

	# Calculate p-values using estimated degrees of freedom
	# and posterior variance estimates
	.standard_transform( fit_eb, sigma = sqrt(fit_eb$s2.post) ) 
})



# # change_t_df
# #
# # change_t_df
# #
# # @param t t-statistic
# # @param df current degrees of freedom
# # @param d_target target degrees of freedom
# #
# # @return Given a t-statistic with degrees of freedom df, return a new t-statistic corresponding to a target degrees of freedom.  Both t-statistics give the same p-value when compared to their respective degrees of freedom
# .change_t_df = function(t, df, d_target){
#   p = pt(abs(t), df, lower.tail=FALSE)
#   qt(p, d_target, lower.tail=FALSE)
# }

# # standardized_t_stat
# #
# # standardized_t_stat
# #
# # @param fit model fit from dream()
# #
# # @return Given a t-statistic with degrees of freedom df, return a new t-statistic corresponding to a target degrees of freedom.  Both t-statistics give the same p-value when compared to their respective degrees of freedom
# #
# # @details Since dream() used degrees of freedom estimated from the data, the t-statistics across genes have different df values used to compute p-values.  Therefore the order of the raw t-statistics doesn't correspond to the order of the p-values since the df is different for each gene.  This function resolved this issue by transofrming the t-statistics to all have the same df.  Under a fixed effects model, df = N - p, where N is the sample size and p is the number of covariates.  Here, df is the target degrees of freedom used in the transformation
# .standardized_t_stat = function( fit ){

#   res = data.frame(t=as.numeric(fit$t), df.sat=as.numeric(fit$df.residual))

#   # d_target = max(res$df.sat)
#   d_target = nrow(fit$design) - ncol(fit$design)

#   res$t_dot = sign(res$t) * .change_t_df( res$t, res$df.sat, d_target )

#   fit2 = fit
#   fit2$t = matrix(res$t_dot, nrow=nrow(fit$t))
#   rownames( fit2$t) = rownames( fit$t)
#   colnames( fit2$t) = colnames( fit$t)

#   # fit2$df.residual = rep(d_target, nrow(res))

#   # if df.prior is defined, set new df.total, since df.residual has changed
#   if( !is.null(fit2$df.prior) ){
# 	fit2$df.total = fit2$df.residual + fit2$df.prior
#   }

#   fit2
# }


#' Compare p-values from two analyses
#'
#' Plot -log10 p-values from two analyses and color based on donor component from variancePartition analysis
#'
#' @param p1 p-value from first analysis
#' @param p2 p-value from second analysis
#' @param vpDonor donor component for each gene from variancePartition analysis
#' @param dupcorvalue scalar donor component from duplicateCorrelation
#' @param fraction fraction of highest/lowest values to use for best fit lines
#' @param xlabel for x-axis
#' @param ylabel label for y-axis
#'
#' @return ggplot2 plot
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' library(BiocParallel)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Perform very simple analysis for demonstration
#' 
#' # Analysis 1
#' form <- ~ Batch 
#' fit = dream( geneExpr, form, info)
#' fit = eBayes( fit )
#' res = topTable( fit, number=Inf, coef="Batch3" )
#' 
#' # Analysis 2
#' form <- ~ Batch + (1|Tissue)
#' fit2 = dream( geneExpr, form, info)
# fitEB2 = eBayes( fit2 )
#' res2 = topTable( fit2, number=Inf, coef="Batch3" )
#' 
#' # Compare p-values
#' plotCompareP( res$P.Value, res2$P.Value, runif(nrow(res)), .3 )
#' 
# # stop cluster
# stopCluster(cl)
#'
#' @export
# @docType methods
#' @rdname plotCompareP-method
plotCompareP = function( p1, p2, vpDonor, dupcorvalue, fraction=.2, xlabel=bquote(duplicateCorrelation~(-log[10]~p)), ylabel=bquote(dream~(-log[10]~p))){

	if( length(unique(c(length(p1), length(p2), length(vpDonor)))) != 1){
		stop("p1, p2 and vpDonor must have the same number of entries")
	}
	if( length(dupcorvalue) != 1){
		stop("dupcorvalue must be a scalar")
	}

	df2 = data.frame(p1=-log10(p1), p2=-log10(p2), vpDonor=vpDonor, delta = vpDonor - dupcorvalue)

	N = nrow(df2)

	c1 = sort(df2$delta)[N*fraction]
	c2 = sort(df2$delta)[length(df2$delta) - N*fraction]

	l1 = lm(p2 ~ p1, df2[df2$delta >= c2,])
	l2 = lm(p2 ~ p1, df2[df2$delta <= c1,])

	df_line = data.frame(rbind(coef(l1), coef(l2)))
	colnames(df_line) = c('a', 'b')
	df_line$type = c('darkred', 'navy')

	lim = c(0, max(max(df2$p1), max(df2$p2)))

	# xlab("duplicateCorrelation (-log10 p)") + ylab("dream (-log10 p)")
	ggplot(df2, aes(p1, p2, color = vpDonor)) + geom_abline() + geom_point(size=2) + theme_bw(17) +  theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5)) + xlim(lim) + ylim(lim) + xlab(xlabel) + ylab(ylabel) + 
		geom_abline( intercept=df_line$a, slope=df_line$b, color=df_line$type, linetype=2) + scale_color_gradientn(name = "Donor", colours = c("blue","green","red"), 
	                       values = rescale(c(0, dupcorvalue, 1)),
	                       guide = "colorbar", limits=c(0, 1))
}

#' Multiple Testing Genewise Across Contrasts
#' 
#'  For each gene, classify a series of related t-statistics as up, down or not significant.
#' 
#' @param object numeric matrix of t-statistics or an 'MArrayLM2' object from which the t-statistics may be extracted.
#' @param ... additional arguments
#' 
#' @details Works like limma::classifyTestsF, except object can have a list of covariance matrices object$cov.coefficients.list, instead of just one in object$cov.coefficients
#' @seealso \code{limma::classifyTestsF}
# @export
setGeneric("classifyTestsF", signature="object",
  function(object, ...)
      standardGeneric("classifyTestsF")
)



#' Multiple Testing Genewise Across Contrasts
#' 
#'  For each gene, classify a series of related t-statistics as up, down or not significant.
#' 
#' @param object numeric matrix of t-statistics or an 'MArrayLM2' object from which the t-statistics may be extracted.
#' @param cor.matrix covariance matrix of each row of t-statistics.  Defaults to the identity matrix.
#' @param df numeric vector giving the degrees of freedom for the t-statistics.  May have length 1 or length equal to the number of rows of tstat.
#' @param p.value numeric value between 0 and 1 giving the desired size of the test
#' @param fstat.only logical, if 'TRUE' then return the overall F-statistic as for 'FStat' instead of classifying the test results
#' 
#' @details Works like limma::classifyTestsF, except object can have a list of covariance matrices object$cov.coefficients.list, instead of just one in object$cov.coefficients
#' @seealso \code{limma::classifyTestsF}
#' @importFrom stats qf
# @export
setMethod("classifyTestsF", "MArrayLM2",
  function(object,cor.matrix=NULL, df=Inf, p.value=0.01, fstat.only=FALSE) {
#	Use F-tests to classify vectors of t-test statistics into outcomes
#	Gordon Smyth
#	20 Mar 2003.  Last revised 6 June 2009.

#	Method intended for MArrayLM objects but accept unclassed lists as well
	if(is.list(object)) {
		if(is.null(object$t)) stop("tstat cannot be extracted from object")
		computeCorrMat = ifelse(is.null(cor.matrix), TRUE, FALSE)
		if(missing(df) && !is.null(object$df.prior) && !is.null(object$df.residual)){
			df <- object$df.prior + object$df.residual
		}
		tstat <- as.matrix(object$t)
	} else {
		tstat <- as.matrix(object)
	}
	ngenes <- nrow(tstat)
	ntests <- ncol(tstat)
	if(ntests == 1) {
		if(fstat.only) {
			fstat <- tstat^2
			attr(fstat,"df1") <- 1
			attr(fstat,"df2") <- df
			return(fstat)
		} else {
			p <- 2 * pt(abs(tstat), df, lower.tail=FALSE)
			return(new("TestResults", sign(tstat) * (p < p.value) ))
		}
	}

	# F test of multiple coefficients
	#-------------------------------

	if( ngenes != length(object$cov.coefficients.list) ){
		msg = paste0("Number of genes does not equal number of elements in cov.coefficients.list\n", ngenes, " != ", length(object$cov.coefficients.list))
		stop(msg)
	}

	fstat <- rep(NA, ngenes) 
	names(fstat) = rownames(tstat)

	result <- matrix(0,ngenes,ntests,dimnames=dimnames(tstat))

	for(i in seq_len(ngenes) ){

		if( computeCorrMat ){
			if( is.null(object$cov.coefficients.list) ){
				C <- object$cov.coefficients
			}else{			
				C <- object$cov.coefficients.list[[i]]
			}
			# subset based on coefficient names
			C = C[colnames(object),colnames(object)]
			cor.matrix <- cov2cor( C )
		}else{
			cor.matrix = cov2cor(cor.matrix)
		}

		# cor.matrix is estimated correlation matrix of the coefficients
		# and also the estimated covariance matrix of the t-statistics
		if(is.null(cor.matrix)){
			r <- ntests
			Q <- diag(r)/sqrt(r)
		}else{
			E <- eigen(cor.matrix,symmetric=TRUE)
			r <- sum(E$values/E$values[1] > 1e-8)
			Q <- limma:::.matvec( E$vectors[,1:r], 1/sqrt(E$values[1:r]))/sqrt(r)
		}

		# Return overall moderated F-statistic only
		if(fstat.only) {
			fstat[i] <- drop( (tstat[i,,drop=FALSE] %*% Q)^2 %*% array(1,c(r,1)) )
		}
		if( i == 1){					
			attr(fstat,"df1") <- r
			attr(fstat,"df2") <- df[i]
		}

		# Return TestResults matrix
		qF <- qf(p.value, r, df[i], lower.tail=FALSE)
		if(length(qF)==1) qF <- rep(qF,ngenes) 
		x <- tstat[i,]
		if(any(is.na(x)))
			result[i,] <- NA
		else
			if( crossprod(crossprod(Q,x)) > qF[i] ) {
				ord <- order(abs(x),decreasing=TRUE)
				result[i,ord[1]] <- sign(x[ord[1]])
				for (j in 2:ntests) {
					bigger <- ord[1:(j-1)]
					x[bigger] <- sign(x[bigger]) * abs(x[ord[j]])
					if( crossprod(crossprod(Q,x)) > qF[i] )
						result[i,ord[j]] <- sign(x[ord[j]])
					else
						break
				}
			}
	}

	if(fstat.only){
		return( fstat )
	}

	new("TestResults",result)
})









