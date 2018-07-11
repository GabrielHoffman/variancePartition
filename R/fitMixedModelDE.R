# #' Class MArrayLMM_lmer.
# #'
# #' Class \code{MArrayLMM_lmer} 
# #'
# #' @name MArrayLMM_lmer-class
# #' @rdname MArrayLMM_lmer-class
# #' @exportClass MArrayLMM_lmer
# # setClass("MArrayLMM_lmer", representation(object="MArrayLM", contrast="numeric", pValue="numeric"))
# # setClass("MArrayLMM_lmer", representation(pValue="numeric"), contains="MArrayLM")

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
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
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
#' object containing MArrayLM object from limma, along with the contrast matrix, and the directly estimated p-value (without eBayes)
# MArrayLMM_lmer object containing MArrayLM object from limma, along with the contrast matrix, and the directly estimated p-value (without eBayes)
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
		if( ! identical(names(L), names(fixef(fitInit))) ){
			stop("Names of entries in L must match fixed effects")
		}

		cons = lmerTest::contest(fitInit, L, ddf=ddf)
		df = as.numeric(cons['DenDF'])

		if(ddf == "Kenward-Roger"){
			# KR
			V = pbkrtest::vcovAdj.lmerMod(fitInit, 0)
		}else{
			# Satterthwaite
			V = vcov(fitInit)
		}

		sigma = attr(lme4::VarCorr(fitInit), "sc")	
		beta = as.matrix(sum(L * fixef(fitInit)), ncol=1)
		SE = as.matrix(sqrt(sum(L * (V %*% L))), ncol=1)
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
		b = names(L)

		if( ! identical( a,b) ){
			stop("Terms in contrast matrix L do not match model:\n  Model: ", paste(a, collapse=',' ), "\n  L: ", paste(b, collapse=',' ), "\nNote thhat order must be the same")
		}

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=gene14643$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		resList <- foreach(gene14643=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4", "lmerTest", "pbkrtest") ) %dopar% {

			# modify data2 for this gene
			data2$expr = gene14643$E
 
			# fit linear mixed model
			fit = lmerTest::lmer( eval(parse(text=form)), data=data2, ..., REML=REML, weights=gene14643$weights, start=fitInit@theta, control=control,na.action=stats::na.exclude)

			cons = lmerTest::contest(fit, L, ddf=ddf)
			df = as.numeric(cons['DenDF'])

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
			beta = as.matrix(sum(L * fixef(fit)), ncol=1)
			colnames(beta) = "logFC"

			SE = as.matrix(sqrt(sum(L * (V %*% L))), ncol=1)		
			colnames(SE) = "logFC"

			# pValue = 2*pt(as.numeric(abs(beta / SE)), df, lower.tail=FALSE)
			pValue = as.numeric(cons['Pr(>F)'])
			
			ret = list(coefficients = beta, 
				design = fit@pp$X, 
				df.residual = df, 
				Amean = mean(fit@frame[,1]), 
				method = 'lmer',
				sigma = sigma,
				stdev.unscaled = SE/sigma,
				pValue = pValue)

			new("MArrayLM", ret)
		}

		coefficients = t(t(sapply( resList, function(x) x$coefficients)))
		colnames(coefficients) = colnames(resList[[1]]$coefficients)
		rownames(coefficients) = rownames(exprObj)

		design = resList[[1]]$design

		df.residual = sapply( resList, function(x) x$df.residual)
		names(df.residual ) = rownames(exprObj)

		pValue = sapply( resList, function(x) x$pValue)
		names(pValue) = rownames(exprObj)

		Amean = sapply( resList, function(x) x$Amean)
		names(Amean ) = rownames(exprObj)
		method = "lmer"
		sigma = sapply( resList, function(x) x$sigma)
		names(sigma) = rownames(exprObj)

		stdev.unscaled  = t(t(sapply( resList, function(x) x$stdev.unscaled )))
		colnames(stdev.unscaled ) = colnames(resList[[1]]$stdev.unscaled )
		rownames(stdev.unscaled ) = rownames(exprObj)

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
	}
	
	# new("MArrayLMM_lmer", object=ret, contrast=L, pValue=pValue)	
	# new("MArrayLMM_lmer", ret, pValue=pValue)
	ret
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
plotCompareP = function( p1, p2, vpDonor, dupcorvalue, fraction=.2){

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
	ggplot(df2, aes(p1, p2, color = vpDonor)) + geom_abline() + geom_point(size=2) + theme_bw(12) +  
		theme(aspect.ratio=1) + xlab("p1 (-log10 p)") + ylab("p2 (-log10 p)") + xlim(lim) + ylim(lim) +
		geom_abline( intercept=df_line$a, slope=df_line$b, color=df_line$type, linetype=2) + scale_color_gradientn(name = "Donor", colours = c("blue","green","red"), 
	                       values = rescale(c(0, dupcorvalue, 1)),
	                       guide = "colorbar", limits=c(0, 1))
}

# #' eBayes for MArrayLMM_lmer
# #'
# #' eBayes for MArrayLMM_lmer
# #'
# #' @param fit fit
# #' @param proportion proportion
# #' @param stdev.coef.lim stdev.coef.lim
# #' @param trend trend
# #' @param robust robust
# #' @param winsor.tail.p winsor.tail.p 
# #'
# #' @return resold of eBayes
# #' @export
# #' @rdname eBayes-method
# #' @aliases eBayes,MArrayLMM_lmer-method
# setMethod("eBayes", "MArrayLMM_lmer",
# function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
#     trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)){
# 	ret = eBayes( fit@object, proportion=proportion, stdev.coef.lim =stdev.coef.lim, trend=trend, robust=robust, winsor.tail.p =winsor.tail.p   )	
# 	new("MArrayLMM_lmer", object=ret, contrast=fit@contrast)
# })

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

#' Example dataset for dream
#'
#' Simulated RNA-seq counts in 'countMatrix' and information in 'metadata'
#'
# "countMatrix"











