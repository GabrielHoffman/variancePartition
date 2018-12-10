
#' Class MArrayLM2.
#'
#' Class \code{MArrayLM2} 
#'
#' @name MArrayLM2-class
#' @rdname MArrayLM2-class
#' @exportClass MArrayLM2
setClass("MArrayLM2",
#  Linear model fit
representation("MArrayLM")
)
setIs("MArrayLM2","LargeDataObject")
setAs(from='MArrayLM', to='MArrayLM2', function(from){
	structure(from, class="MArrayLM2")
	})

#' Extract contrast matrix for linear mixed model
#' 
#' Extract contrast matrix, L, testing a single variable.  Contrasts involving more than one variable can be constructed by modifying L directly
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)  Formulas with only fixed effects also work
#' @param data data.frame with columns corresponding to formula 
#' @param coefficient the coefficient to use in the hypothesis test
#' 
#' @return
#' Contrast matrix testing one variable
#'
#' @examples
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # get contrast matrix testing if the coefficient for Batch2 is zero 
#' # The variable of interest must be a fixed effect
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L = getContrast( geneExpr, form, info, "Batch3")
#'
#' # get contrast matrix testing if Batch3 - Batch2 = 0
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L = getContrast( geneExpr, form, info, c("Batch3", "Batch2"))
#'
#' # To test against Batch1 use the formula:
#' # 	~ 0 + Batch + (1|Individual) + (1|Tissue) 
#' # to estimate Batch1 directly instead of using it as the baseline
#'
#' @export
#' @docType methods
#' @rdname getContrast-method
getContrast = function( exprObj, formula, data, coefficient){ 

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

	if( length(coefficient) > 2){
		stop("Length of coefficient array limited to 2")
	}

	REML=TRUE
	useWeights=TRUE
	weightsMatrix=NULL
	showWarnings=FALSE
	dream=TRUE
	fxn=identity
	colinearityCutoff=.999
	control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )

	# check dimensions of reponse and covariates
	if( ncol(exprObj) != nrow(data) ){		
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e. rows)" )
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

	# If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
	if( ! identical(colnames(exprObj), rownames(data)) ){
		 warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
	}

	# add response (i.e. exprObj[,j] to formula
	form = paste( "gene14643$E", paste(as.character( formula), collapse=''))

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene14643 = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,control=control ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	method = ''
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){
		
		design = model.matrix( formula, data)

		L = rep(0, ncol(design))
		names(L) = colnames(design)

	}else{
		fit = lmer( eval(parse(text=form)), data=data,control=control, REML=TRUE )
		
		L = rep(0, length(fixef(fit)))
		names(L) = names(fixef(fit))
	}

	# assign coefficient coding
	if( any(!coefficient %in% names(L)) ){
		stop("coefficient is not in the formula.  Valid coef are:\n", paste(names(L), collapse=', '))
	}
	L[coefficient[1]] = 1

	if( length(coefficient) == 2){			
		L[coefficient[2]] = -1
	}
	
	L
}

# Evaluate contrasts for linear mixed model
# 
# Evaluate contrasts for linear mixed model
#
# @param fit model fit
# @param L contrast matrix
# @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough exproximation for most datasets.
# 
# @return
# df, sigma, beta, SE of model
#
# @export
# @docType methods
# @rdname eval_lmm-method
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
		# df = pbkrtest::get_Lb_ddf(fit, L)
	}else{
		# Satterthwaite
		V = vcov(fit)
		# df = as.numeric(contest(fit, L, ddf="Sat")['DenDF'])
	}

	sigma = attr(lme4::VarCorr(fit), "sc")

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

	list(cons = cons,
		df		= df,
		sigma	= sigma,
		beta	= beta,
		SE		= SE,
		pValue	= pValue)
}


#' Differential expression with linear mixed model
#' 
#' Fit linear mixed model for differential expression and preform hypothesis test on fixed effects as specified in the contrast matrix L
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)  Formulas with only fixed effects also work, and lmFit() followed by contrasts.fit() are run.
#' @param data data.frame with columns corresponding to formula 
#' @param L contrast matrix specifying a linear combination of fixed effects to test
#' @param ddf Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough exproximation for most datasets.
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is TRUE.  Strongly discourage against changing this option
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param control control settings for lmer()
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return 
#' MArrayLM2 object (just like MArrayLM from limma), and the directly estimated p-value (without eBayes)
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is ~ a + b + (1|c), then to model is 
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' useWeights=TRUE causes weightsMatrix[j,] to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using foreach/dopar to run loops in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lmer.
#'
#' Hypothesis tests and degrees of freedom are producted by lmerTest and pbkrtest pacakges
#' @examples
#'
#' # load library
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
#' # get contrast matrix testing if the coefficient for Batch2 is zero 
#' # The variable of interest must be a fixed effect
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L = getContrast( geneExpr, form, info, "Batch3")
#' 
#' # Fit linear mixed model for each gene
#' # run on just 10 genes for time
#' fit = dream( geneExpr[1:10,], form, info, L)
#' 
#' # Run empirical Bayes post processing from limma
#' fitEB = eBayes( fit )
#' 
#' # view top genes
#' topTable( fitEB )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
#' @rdname dream-method
dream <- function( exprObj, formula, data, L, ddf = c("Satterthwaite", "Kenward-Roger"), REML=TRUE, useWeights=TRUE, weightsMatrix=NULL,control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...){ 

	exprObjInit = exprObj

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )
	ddf = match.arg(ddf)
	colinearityCutoff=.999

	# check dimensions of reponse and covariates
	if( ncol(exprObj) != nrow(data) ){		
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
	}

	# check if all genes have variance
	rv = apply( exprObj, 1, var)
	if( any( rv == 0) ){
		idx = which(rv == 0)
		stop(paste("Response variable", idx[1], 'has a variance of 0'))
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

	# If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
	if( ! identical(colnames(exprObj), rownames(data)) ){
		 warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
	}

	# format contrasts 
	if( class(L) == "numeric" ){
		L = as.matrix(L, ncol=1)
	}

	# add response (i.e. exprObj[,j] to formula
	form = paste( "gene14643$E", paste(as.character( formula), collapse=''))

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene14643 = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,...,control=control ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	method = ''
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){
		cat("Fixed effect model, using limma directly...\n")

		design = model.matrix( formula, data)
		fit = lmFit( exprObj, design )
		ret = contrasts.fit( fit, L)

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene14643 = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmerTest::lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )
		
		# check L
		if( ! identical(rownames(L), names(fixef(fitInit))) ){
			stop("Names of entries in L must match fixed effects")
		}
		
		# extract statistics from model
		mod = .eval_lmm( fitInit, L, ddf)
		timediff = proc.time() - timeStart

		# check size of stored objects
		objSize = object.size( fitInit ) * nrow(exprObj)

		# total time = (time for 1 gene) * (# of genes) / 60 / (# of threads)
		showTime = timediff[3] * nrow(exprObj) / 60 / getDoParWorkers()

		# cat("Projected memory usage: >", format(objSize, units = "auto"), "\n")

		if( showTime > .01 ){
			cat("Projected run time: ~", paste(format(showTime, digits=1), "min"), "\n")
		}

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=FALSE, dream=TRUE, colinearityCutoff )

		a = names(fixef(fitInit))
		b = rownames(L)

		if( ! identical( a,b) ){
			stop("Terms in contrast matrix L do not match model:\n  Model: ", paste(a, collapse=',' ), "\n  L: ", paste(b, collapse=',' ), "\nNote thhat order must be the same")
		}

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=gene14643$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",,
			total = nrow(exprObj), width= 60, clear=FALSE)

		pids = .get_pids()

		timeStart = proc.time()

		# loop through genes
		resList <- foreach(gene14643=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4", "lmerTest", "pbkrtest"), .export='.eval_lmm' ) %dopar% {

			# modify data2 for this gene
			data2$expr = gene14643$E
 
			# fit linear mixed model
			fit = lmerTest::lmer( eval(parse(text=form)), data=data2, REML=REML,..., weights=gene14643$weights, start=fitInit@theta, control=control,na.action=stats::na.exclude)

			# extract statistics from model
			mod = .eval_lmm( fit, L, ddf)

			# progressbar
			if( Sys.getpid() == pids[1]){
				pb$update( gene14643$n_iter / gene14643$max_iter )
			}

			ret = list(coefficients = mod$beta, 
				design = fit@pp$X, 
				df.residual = mod$df, 
				Amean = mean(fit@frame[,1]), 
				method = 'lmer',
				sigma = mod$sigma,
				stdev.unscaled = mod$SE/mod$sigma,
				pValue = mod$pValue)

			new("MArrayLM", ret)
		}
		cat("\nFinished...")
		cat("\nTotal:", paste(format((proc.time() - timeStart)[3], digits=0), "s\n"))		

		x = 1
		# extract results
		coefficients = foreach( x = resList, .combine=cbind ) %do% {x$coefficients}
		df.residual = foreach( x = resList, .combine=cbind ) %do% {	x$df.residual}
		pValue = foreach( x = resList, .combine=cbind ) %do% { x$pValue }
		stdev.unscaled  = foreach( x = resList, .combine=cbind ) %do% {x$stdev.unscaled} 
		
		# transpose
		coefficients = t( coefficients )
		df.residual = t( df.residual )
		pValue = t( pValue )
		stdev.unscaled = t( stdev.unscaled )		

		colnames(coefficients) = colnames(L)
		rownames(coefficients) = rownames(exprObj)

		design = resList[[1]]$design

		colnames(df.residual) = colnames(L)
		rownames(df.residual) = rownames(exprObj)

		colnames(pValue) = colnames(L)
		rownames(pValue) = rownames(exprObj)

		Amean =sapply( resList, function(x) x$Amean) 
		names(Amean ) = rownames(exprObj)
		method = "lmer"
		sigma = sapply( resList, function(x) x$sigma)
		names(sigma) = rownames(exprObj)

		colnames(stdev.unscaled) = colnames(L)
		rownames(stdev.unscaled) = rownames(exprObj)

		ret = list( coefficients = coefficients,
		 			design = design, 
		 			df.residual = df.residual, 
		 			Amean = Amean, 
		 			method = method, 
		 			sigma = sigma, 
		 			contrasts = L,
		 			stdev.unscaled = stdev.unscaled)

		if( 'genes' %in% names(exprObjInit) ){
			ret$genes = exprObjInit$genes
		}

		ret = new("MArrayLM", ret)
		ret$pValue = pValue
		ret = as(ret, "MArrayLM2")
	}
	ret	
}


# 
# 

#' Subseting for MArrayLM2
#'
#' Enable subsetting on MArrayLM2 object. Same as for MArrayLM, but apply column subsetting to df.residual and pValue
#'
#' @param object MArrayLM2
#' @param i row
#' @param j col
#'
#' @name [.MArrayLM2
#' @return subset
#' @export
# @rdname [.MArrayLM2-method
# @aliases [.MArrayLM2,MArrayLM2-method
assign("[.MArrayLM2",
	function(object, i, j){
		if(nargs() != 3){
			stop("Two subscripts required",call.=FALSE)
		}

		obj = as(object, 'MArrayLM')

		if(!missing(j)){
			obj = obj[,j]
		}
		if(!missing(i)){
			obj = obj[i,]
		}
		obj$df.residual = object$df.residual[i,j]
		obj$pValue = object$pValue[i,j]
		obj
		})

#' eBayes for MArrayLM2
#'
#' eBayes for MArrayLM2
#'
#' @param fit fit
#' @param proportion proportion
#' @param stdev.coef.lim stdev.coef.lim
#' @param trend trend
#' @param robust robust
#' @param winsor.tail.p winsor.tail.p 
#'
#' @return resold of eBayes
#' @export
#' @rdname eBayes-method
#' @aliases eBayes,MArrayLM2-method
setMethod("eBayes", "MArrayLM2",
function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)){

	i = 1
	retList = foreach( i = 1:ncol(fit) ) %do% {

		ret = limma::eBayes( fit[,i], proportion=proportion, stdev.coef.lim =stdev.coef.lim, trend=trend, robust=robust, winsor.tail.p =winsor.tail.p )

		# transform moderated t-statistics to have same degrees of freedom
		.standardized_t_stat( ret )	
	}

	fit2 = retList[[1]]

	fit2$coefficients = do.call("cbind", lapply(retList, function(fit) fit$coefficients))
	fit2$contrasts = do.call("cbind", lapply(retList, function(fit) fit$contrasts))
	fit2$stdev.unscaled = do.call("cbind", lapply(retList, function(fit) fit$stdev.unscaled))
	fit2$pValue = do.call("cbind", lapply(retList, function(fit) fit$pValue))
	fit2$df.prior = do.call("cbind", lapply(retList, function(fit) fit$df.prior))
	fit2$s2.prior = do.call("cbind", lapply(retList, function(fit) fit$s2.prior))
	fit2$var.prior = do.call("cbind", lapply(retList, function(fit) fit$var.prior))
	fit2$s2.post = do.call("cbind", lapply(retList, function(fit) fit$s2.post))
	fit2$t = do.call("cbind", lapply(retList, function(fit) fit$t))
	fit2$df.total = do.call("cbind", lapply(retList, function(fit) fit$df.total))
	fit2$p.value = do.call("cbind", lapply(retList, function(fit) fit$p.value))
	fit2$lods = do.call("cbind", lapply(retList, function(fit) fit$lods))

	colnames(fit2$pValue) = colnames(fit2$coefficients)
	colnames(fit2$t) = 	colnames(fit2$coefficients)
	colnames(fit2$df.total) = colnames(fit2$coefficients)

	fit2 
})


# #' topTable for MArrayLMM_lmer
# #'
# #' topTable for MArrayLMM_lmer
# #'
# #' @param fit fit
# #' @param coef coef
# #' @param number number
# #' @param genelist genelist
# #' @param adjust.method adjust.method
# #' @param sort.by sort.by
# #' @param resort.by resort.by
# #' @param p.value p.value
# #' @param lfc lfc
# #' @param confint confint
# #'
# #' @return resold of topTable
# #' @export
# #' @rdname topTable-method
# #' @aliases topTable,MArrayLMM_lmer-method
# setMethod("topTable", "MArrayLMM_lmer",
# function(fit, coef=NULL, number=10, genelist=fit$genes, adjust.method="BH",
#               sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE){
# 	topTable( fit@object, number=number,  adjust.method=adjust.method, sort.by= sort.by, resort.by=resort.by, p.value=p.value, lfc=lfc, confint=confint  )	
# })


# #' toptable for MArrayLMM_lmer
# #'
# #' toptable for MArrayLMM_lmer
# #'
# #' @param fit fit
# #' @param coef coef
# #' @param number number
# #' @param genelist genelist
# #' @param A A
# #' @param eb eb
# #' @param adjust.method adjust.method
# #' @param sort.by sort.by
# #' @param resort.by resort.by
# #' @param p.value p.value
# #' @param lfc lfc
# #' @param confint confint
# #' @param ... ...
# #'
# #' @return resold of toptable
# #' @export
# #' @rdname toptable-method
# #' @aliases toptable,MArrayLMM_lmer-method
# setMethod("toptable", "MArrayLMM_lmer",
# function(fit, coef=1, number=10, genelist=NULL, A=NULL, eb=NULL, adjust.method="BH",
#               sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE,...){
# 	topTable( fit@object, number=number,  adjust.method=adjust.method, sort.by= sort.by, resort.by=resort.by, p.value=p.value, lfc=lfc, confint=confint  )	
# })

# #' Example dataset for dream
# #'
# #' Simulated RNA-seq counts in 'countMatrix' and information in 'metadata'
# #'
# #' @docType data
# #'
# #' @usage data(varPartDEdata)
# #'
# #' @format matrix of RNA-seq counts
# #'
# #' @examples
# #' data(varPartDEdata)
# "varPartDEdata"



# change_t_df
#
# change_t_df
#
# @param t t-statistic
# @param df current degrees of freedom
# @param d_target target degrees of freedom
#
# @return Given a t-statistic with degrees of freedom df, return a new t-statistic corresponding to a target degrees of freedom.  Both t-statistics give the same p-value when compared to their respective degrees of freedom
.change_t_df = function(t, df, d_target){
  p = pt(abs(t), df, lower.tail=FALSE)
  qt(p, d_target, lower.tail=FALSE)
}

# standardized_t_stat
#
# standardized_t_stat
#
# @param fit model fit from dream()
#
# @return Given a t-statistic with degrees of freedom df, return a new t-statistic corresponding to a target degrees of freedom.  Both t-statistics give the same p-value when compared to their respective degrees of freedom
#
# @details Since dream() used degrees of freedom estimated from the data, the t-statistics across genes have different df values used to compute p-values.  Therefore the order of the raw t-statistics doesn't correspond to the order of the p-values since the df is different for each gene.  This function resolved this issue by transofrming the t-statistics to all have the same df.  Under a fixed effects model, df = N - p, where N is the sample size and p is the number of covariates.  Here, df is the target degrees of freedom used in the transformation
.standardized_t_stat = function( fit ){
  res = data.frame(t=as.numeric(fit$t), df.sat=fit$df.residual)

  # d_target = max(res$df.sat)
  d_target = nrow(fit$design) - ncol(fit$design)

  res$t_dot = sign(res$t) * .change_t_df( res$t, res$df.sat, d_target )

  fit$t = res$t_dot 
  fit$df.residual = rep(d_target, nrow(res))
  fit
}



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
#' # Perform very simple analysis for demonstration
#' 
#' # Analysis 1
#' form <- ~ Batch 
#' L = getContrast( geneExpr, form, info, "Batch3")
#' fit = dream( geneExpr, form, info, L)
#' fitEB = eBayes( fit )
#' res = topTable( fitEB, number=Inf )
#' 
#' # Analysis 2
#' form <- ~ Batch + (1|Tissue)
#' L = getContrast( geneExpr, form, info, "Batch3")
#' fit = dream( geneExpr, form, info, L)
#' fitEB = eBayes( fit )
#' res2 = topTable( fitEB, number=Inf )
#' 
#' # Compare p-values
#' plotCompareP( res$P.Value, res2$P.Value, runif(nrow(res)), .3 )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
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
	ggplot(df2, aes(p1, p2, color = vpDonor)) + geom_abline() + geom_point(size=2) + theme_bw(12) +  theme(aspect.ratio=1,plot.title = element_text(hjust = 0.5)) + xlim(lim) + ylim(lim) + xlab(xlabel) + ylab(ylabel) + 
		geom_abline( intercept=df_line$a, slope=df_line$b, color=df_line$type, linetype=2) + scale_color_gradientn(name = "Donor", colours = c("blue","green","red"), 
	                       values = rescale(c(0, dupcorvalue, 1)),
	                       guide = "colorbar", limits=c(0, 1))
}


