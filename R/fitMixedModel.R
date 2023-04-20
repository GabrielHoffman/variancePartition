
## default control for lmer fits
vpcontrol <- lme4::lmerControl(calc.derivs = FALSE,
                               check.rankX = "stop.deficient",
                               check.conv.singular =
                                 lme4::.makeCC("ignore", tol = 1e-4))

#' @importFrom lme4 lmer refit
#' @importFrom stats update
run_lmm_on_gene = function(obj, formula, data, control, na.action, REML, fxn, fit.init=NULL, dreamCheck = FALSE, varTol= 1e-5){

	# if sparseMatrix, convert to matrix
	if( is(obj$E, "dgCMatrix") ){
		obj$E = as.matrix(obj$E)
	}

	# append data
	w.local = NULL
	data$y.local = c(obj$E)
	data$w.local = c(obj$weights)

	# scale regression weights
	data$w.local = data$w.local / mean(data$w.local)

	form.local = update(formula, y.local ~ .)

	if( var(data$y.local, na.rm=TRUE) < varTol ){
		stop(paste0("response variance < ", varTol))
	} 

	if( .isMixedModelFormula( formula ) ){
		if( ! is.null(fit.init) ){
			# if fit.init is passed, use refit
			fit = refit(fit.init, 
				newresp = data$y.local, 
				newweights = data$w.local,
				control = control)
		}else{
			# fit linear mixed model from scratch
			fit = lmer(form.local, data, 
				weights = w.local, 
				control = control, 
				na.action = na.action,
				REML = REML)
		}
	}else{
		fit = lm(form.local, data, weights = w.local, 
				na.action = na.action)
	}	

	# check fit for errors
	checkModelStatus( fit, dreamCheck)

	# run post-processing function on model fit
	res = fxn( fit )

	res
}


# run analysis on each batch
#' @importFrom BiocParallel bpstopOnError<- bptry
run_lmm_on_batch = function(obj, form, data, control, na.action, REML, fxn, fit.init=NULL, dreamCheck = FALSE, varTol= 1e-5, BPPARAM = SerialParam()){

	# create list with one gene per entry
	exprList = lapply(seq(nrow(obj)), function(j){
			new("EList", list(
				E = obj$E[j,,drop=FALSE], 
				weights = obj$weights[j,,drop=FALSE]))
			})

	# run iterator	
	# bplapply is much faster for small chunks than 
	# bpiterate using an iterator
	bpstopOnError(BPPARAM) <- FALSE
	res = bptry(bplapply(exprList, run_lmm_on_gene, 
		form = form, 
		data = data, 
		control = control,
		na.action = na.action,
		REML = REML,
		fxn = fxn,
		fit.init = fit.init, 
		dreamCheck = dreamCheck,
		varTol = varTol,
		BPPARAM = SerialParam()))

	# get failed jobs
	failedJobs = (!bpok(res)) | sapply(res, is.null)

	names(res) = rownames(obj)

	# get error messages as text
	errorText = attr(res, "errors")
	errorText = sapply(errorText, function(x) x$message)

	# assign gene names to the error text
	if( !is.null(errorText) ){
		idx = as.numeric(names(errorText))
		names(errorText) = names(res)[idx]
	}

	list(succeeded = res[!failedJobs],
		 errors = errorText)
}

#' @importFrom BiocParallel bpstopOnError<-
#' @importFrom RhpcBLASctl omp_set_num_threads
run_lmm = function( obj, form, data, control = variancePartition:::vpcontrol, fxn, REML = FALSE, useInitialFit = TRUE, dreamCheck = FALSE, varTol= 1e-5, BPPARAM=SerialParam(),... ){

	# only use 1 thread internally
	omp_set_num_threads(1)

	# catch warnings like errors
	options(warn = 2)

	if(any(obj$weights < 0)) stop("All weights must be positive")

	# run single fit to recycle the internal data structure
	# also to run checks on data
	# it.init = iterRows(obj$E, obj$weights, sizeOfChunk=1)
	it.init = iterRows( matrix(seq(ncol(obj)), nrow=1),  sizeOfChunk=1)

	fit.init = run_lmm_on_gene(it.init(), form, data, 
				control = control,
			    na.action = stats::na.exclude,
			    REML = REML,
			    fxn = identity,
			    dreamCheck = dreamCheck)

	# checkModelStatus( fit.init, dreamCheck )

	if( ! useInitialFit ) fit.init = NULL

	# initialize iterator
	it = iterRows(obj$E, obj$weights, 
				sizeOfChunk = ceiling(nrow(obj) / BPPARAM$workers) )

	# iterate over batches	
	bpstopOnError(BPPARAM) <- FALSE
	resList <- bpiterate( it, run_lmm_on_batch,
			                      form = form, 
			                      data = data,
			                      control = control,
			                      na.action = stats::na.exclude,
			                      REML = REML,
			                      fit.init = fit.init,
			                      dreamCheck = dreamCheck,
			                      varTol = varTol,
			                      fxn = fxn,
			                      BPPARAM = BPPARAM,
			                      reduce.in.order = TRUE)

	# reset to enable warnings
	options(warn = 0)

	list( succeeded = lapply(resList, function(x) x$succeeded),
		  errors = unlist(sapply(resList, function(x) x$errors)))
}



filterInputData = function(exprObj, formula, data, useWeights, isCounts = FALSE){
	# convert to data.frame
	data = as.data.frame(data)
	keep = colnames(data) %in% all.vars(formula)
	data = droplevels(data[, keep,drop=FALSE])

	# make sure form is a formula
	formula = stats::as.formula( formula )

	# check that variables in the formula are all in the data
	idx = unique(all.vars(formula)) %in% colnames(data)
	if( any(!idx) ){
		txt = paste(unique(all.vars(formula))[!idx], collapse=', ')
		stop("Variable in formula not found in data: ", txt)
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
	}

	if( ! isCounts ){
		# convert exprObj to EList if not already
		if( ! is(exprObj, "EList") ){
			W = matrix(1, nrow(exprObj), ncol(exprObj))
	    	exprObj = new("EList", list(E = exprObj, weights=W))
	    }

	    if( ! useWeights ){
	    	exprObj$weights[] = 1
	    }
	}

    # If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
	if( ! identical(colnames(exprObj), rownames(data)) ){
		 warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
	}

	list(exprObj = exprObj,
		formula = formula,
		data = data)
}








