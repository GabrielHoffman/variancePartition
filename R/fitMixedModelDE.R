
setClass("MArrayLMM_lmer", representation(object="MArrayLM", contrast="numeric", pValue="numeric"))

#' @export
getContrast = function( exprObj, formula, data, coefficient){ 

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

	REML=FALSE
	useWeights=TRUE
	weightsMatrix=NULL
	showWarnings=TRUE
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
		stop("For fixed effect model, use limma directly")
	}else{
		fit = lmer( eval(parse(text=form)), data=data,control=control )
		
		coef = 'TissueB'

		L = rep(0, length(fixef(fit)))
		names(L) = names(fixef(fit))
		if( !(coefficient %in% names(L)) ){
			stop("coefficient is not in the formula.  Valid coef are:\n", paste(names(L), collapse=', '))
		}
		L[coefficient] = 1
	}
	
	L
}



#' @export
fitMixedModelDE <- function( exprObj, formula, data, L, REML=FALSE, ddf = c("Satterthwaite", "Kenward-Roger"), useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, colinearityCutoff=.999,control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...){ 

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

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
		stop("For fixed effect model, use limma directly")
	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene14643 = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmerTest::lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )
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
		objSize = object.size( fxn(fitInit) ) * nrow(exprObj)

		# total time = (time for 1 gene) * (# of genes) / 60 / (# of threads)
		showTime = timediff[3] * nrow(exprObj) / 60 / getDoParWorkers()

		# cat("Projected memory usage: >", format(objSize, units = "auto"), "\n")

		if( showTime > .01 ){
			cat("Projected run time: ~", paste(format(showTime, digits=1), "min"), "\n")
		}

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff )

		a = names(fixef(fitInit))
		b = names(L)

		if( ! identical( a,b) ){
			stop("Terms in contrast matrix L do not match model:\n  Model: ", paste(a, collapse=',' ), "\n  L: ", paste(b, collapse=',' ), "\nNote thhat order must be the same")
		}

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=gene14643$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		resList <- foreach(gene14643=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4") ) %dopar% {

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
		method = "ls"
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
		 			stdev.unscaled = stdev.unscaled)

		ret = new("MArrayLM", ret)
	}
	
	new("MArrayLMM_lmer", object=ret, contrast=L, pValue=pValue)	
}


setMethod("eBayes", "MArrayLMM_lmer",
function(fit, proportion = 0.01, stdev.coef.lim = c(0.1, 4), 
    trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)){
	ret = eBayes( fit@object, proportion=proportion, stdev.coef.lim =stdev.coef.lim, trend=trend, robust=robust, winsor.tail.p =winsor.tail.p   )	
	new("MArrayLMM_lmer", object=ret, contrast=fit@contrast)
})


setMethod("topTable", "MArrayLMM_lmer",
function(fit, coef=NULL, number=10, genelist=fit$genes, adjust.method="BH",
              sort.by="B", resort.by=NULL, p.value=1, lfc=0, confint=FALSE){
	topTable( fit@object, number=number,  adjust.method=adjust.method, sort.by= sort.by, resort.by=resort.by, p.value=p.value, lfc=lfc, confint=confint  )	
})