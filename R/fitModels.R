
#' Fit linear (mixed) model
#' 
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula 
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  See Details.  
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from voom().  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from \code{voom()}.  Used only if useWeights is TRUE 
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param fxn apply function to model fit for each gene.  Defaults to identify function so it returns the model fit itself
#' @param control control settings for \code{lmer()}
#' @param BPPARAM parameters for parallel evaluation
#' @param quiet suppress message, default FALSE
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#' 
#' @return
#' \code{list()} of where each entry is a model fit produced by \code{lmer()} or \code{lm()}
#' 
#' @importFrom MASS ginv
# @importFrom RSpectra eigs_sym
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom graphics abline axis hist image layout lines mtext par plot plot.new rect text title
#' @importFrom stats anova as.dendrogram as.dist cancor coef cov2cor density dist fitted.values hclust lm median model.matrix order.dendrogram quantile reorder residuals sd terms var vcov pt qt
#' @importFrom scales rescale
#' @importFrom iterators nextElem
#' @import Rdpack
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is 
#'
#' \code{fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)}
#'
#' If there are no random effects, so formula is \code{~ a + b + c}, a 'standard' linear model is used:
#'
#' \code{fit <- lm( exprObj[j,] ~ a + b + c, data=data)}
#'
#' In both cases, \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lm/lmer.
#'
#' Since this function returns a list of each model fit, using this function is slower and uses more memory than \code{fitExtractVarPartModel()}.
#'
#' \code{REML=FALSE} uses maximum likelihood to estimate variance fractions.  This approach produced unbiased estimates, while \code{REML=TRUE} can show substantial bias.  See Vignette "3) Theory and practice of random effects and REML"
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
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expression
#' # If categorical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' # also sort columns
#' plotVarPart( sortCols( varPart ) )
#'
#' # Advanced: 
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results <- fitVarPartModel( geneExpr, form, info )
#' 
#' # Step 2: extract variance fractions
#' varPart <- extractVarPart( results )
#'
#' # Note: fitVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form <- ~ (1|sex) + (1|type) + score
#' info2 <- Biobase::pData(sample.ExpressionSet)
#' results2 <- fitVarPartModel( sample.ExpressionSet, form, info2 )
#' 
# # Parallel processing using multiple cores with reduced memory usage
# param <- SnowParam(4, "SOCK", progressbar=TRUE)
# results2 <- fitVarPartModel( sample.ExpressionSet, form, info2, BPPARAM=param)
#'
#' @export
#' @docType methods
#' @rdname fitVarPartModel-method
setGeneric("fitVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(),...)
      standardGeneric("fitVarPartModel")
)

# internal driver function
#' @importFrom BiocParallel SerialParam bpiterate bplapply bpok
#' @importFrom methods is new
#' @importFrom lme4 lmer 
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom progress progress_bar 
#' @importFrom utils object.size
#' @import foreach
#' @import lme4
.fitVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, colinearityCutoff=.999,control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=quiet, BPPARAM=SerialParam(), ...){

	# convert to data.frame
	data = as.data.frame(data)

	# if( ! is(exprObj, "sparseMatrix")){
	# 	exprObj = as.matrix( exprObj )	
	# }
	formula = stats::as.formula( formula )

	# only retain columns used in the formula
	data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]

	# check dimensions of reponse and covariates
	if( ncol(exprObj) != nrow(data) ){		
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
	}

	# check if variables in formula has NA's
	hasNA = hasMissingData(formula, data)

	if( any(hasNA) ){
		warning(paste("Variables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\nSamples with missing data will be dropped.\n"), immediate.=TRUE)
	}

	# check if all genes have variance
	if( ! is(exprObj, "sparseMatrix")){
		# check if values are NA
		countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
		if( countNA > 0 ){
			stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
		}

		rv = apply( exprObj, 1, var)
	}else{
		rv = c()
		for( i in seq_len(nrow(exprObj)) ){
			rv[i] = var( exprObj[i,])
		}
	}
	if( any( rv == 0) ){
		idx = which(rv == 0)
		stop(paste("Response variable", idx[1], 'has a variance of 0'))
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
	if( .isDisconnected() ){
		stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
	}

	# If samples names in exprObj (i.e. columns) don't match those in data (i.e. rows)
	if( ! identical(colnames(exprObj), rownames(data)) ){
		 warning( "Sample names of responses (i.e. columns of exprObj) do not match\nsample names of metadata (i.e. rows of data).  Recommend consistent\nnames so downstream results are labeled consistently." )
	}

	# add response (i.e. exprObj[j,]) to formula
	# Use term 'responsePlaceholder' to store the value of reponse j (exprObj[j,])
	# This is not an ideal way to structure the evaluation of response j.
	# 	The formula is evaluated a different scope, (i.e. within lmer()), and there is no need to pass the
	# 	entire exprObj object into that scope.  With lexical scope, I found it was possible that  
	# 	the value of exprObj[j,] could be different when evaluated in the lower scope
	# This placeholder term addresses that issue
	form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,...,control=control ), error = function(e) e)

	# detect error when variable in formula does not exist
	if( inherits(possibleError, "error") ){
		err = grep("object '.*' not found", possibleError$message)
		if( length(err) > 0 ){
			stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
		}
	}
	
	pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",,
			total = nrow(exprObj), width= 60, clear=FALSE)

	timeStart = proc.time()

	mesg <- "No random effects terms specified in formula"
	method = ''
	if( isTRUE(inherits(possibleError, "error") && identical(possibleError$message, mesg)) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

		resList <- foreach(responsePlaceholder=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4") ) %do% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=responsePlaceholder$weights,na.action=stats::na.exclude,...)

			# apply function
			fxn( fit )
		}

		method = "lm"

	}else{

		if( isTRUE(inherits(possibleError, "error") &&  grep('the fixed-effects model matrix is column rank deficient', possibleError$message) == 1) ){
			stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
		} 

		# fit first model to initialize other model fits
		# this make the other models converge faster
		responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )

		timediff = proc.time() - timeStart

		# check size of stored objects
		objSize = object.size( fxn(fitInit) ) * nrow(exprObj)

		if( !quiet ) message("Memory usage to store result: >", format(objSize, units = "auto"))

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=responsePlaceholder$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		# Define function for parallel evaluation
		.eval_models = function(responsePlaceholder, data2, form, REML, theta, fxn, control, na.action=stats::na.exclude,...){

			# modify data2 for this gene
			data2$expr = responsePlaceholder$E
 
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data2, ..., REML=REML, weights=responsePlaceholder$weights, control=control,na.action=na.action)

			# apply function
			fxn( fit )
		}

		.eval_master = function( obj, data2, form, REML, theta, fxn, control, na.action=stats::na.exclude,... ){

			# use only 1 OpenMP thread for linear algebra
			omp_set_num_threads(1)

			lapply(seq_len(nrow(obj$E)), function(j){
				.eval_models( list(E=obj$E[j,], weights=obj$weights[j,]), data2, form, REML, theta, fxn, control, na.action,...)
			})
		}

		# Evaluate function
		###################
		
		it = iterBatch(exprObj, weightsMatrix, useWeights, n_chunks = 100, BPPARAM = BPPARAM)
		
		if( !quiet ) message(paste0("Dividing work into ",attr(it, "n_chunks")," chunks..."))

    resList <- bpiterate( it, .eval_master,
			data2=data2, form=form, REML=REML, theta=fitInit@theta, fxn=fxn, control=control,...,
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

		method = "lmer"
	}

	# pb$update( responsePlaceholder$max_iter / responsePlaceholder$max_iter )
	if( !quiet ) message("\nTotal:", paste(format((proc.time() - timeStart)[3], digits = 1, scientific = FALSE), "s"))
	# set name of each entry
	names(resList) <- rownames( exprObj )

 	new( "VarParFitList", resList, method=method )
}

## matrix
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,matrix-method
setMethod("fitVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitVarPartModel(exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM,...)
  }
)

# data.frame
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,data.frame-method
setMethod("fitVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitVarPartModel( as.matrix(exprObj), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
  }
)

## EList 
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,EList-method
setMethod("fitVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitVarPartModel( as.matrix(exprObj$E), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM,...)
  }
)

## ExpressionSet
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,ExpressionSet-method
#' @importFrom Biobase ExpressionSet exprs
setMethod("fitVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitVarPartModel( as.matrix(exprs(exprObj)), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
  }
)

# sparseMatrix
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,sparseMatrix-method
#' @importFrom Matrix sparseMatrix
setMethod("fitVarPartModel", "sparseMatrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitVarPartModel( exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, quiet=quiet, BPPARAM=BPPARAM, ...)
  }
)

#' Fit linear (mixed) model, report variance fractions
#' 
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Report fraction of variance attributable to each variable 
#'
#' @param exprObj matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by \code{voom()} from the \code{limma} package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: \code{~ a + b + (1|c)}
#' @param data \code{data.frame} with columns corresponding to formula 
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.   See Details.
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from \code{voom()}.  Used only if \code{useWeights} is TRUE 
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param control control settings for \code{lmer()}
#' @param quiet suppress message, default FALSE
#' @param BPPARAM parameters for parallel evaluation
#' @param ... Additional arguments for \code{lmer()} or \code{lm()}
#' 
#' @return
#' list() of where each entry is a model fit produced by \code{lmer()} or \code{lm()}
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is 
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' If there are no random effects, so formula is ~ a + b + c, a 'standard' linear model is used:
#'
#' \code{fit <- lm( exprObj[j,] ~ a + b + c, data=data)}
#'
#' In both cases, \code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to \code{lm}/\code{lmer}.
#'
#' \code{REML=FALSE} uses maximum likelihood to estimate variance fractions.  This approach produced unbiased estimates, while \code{REML=TRUE} can show substantial bias.  See Vignette "3) Theory and practice of random effects and REML"
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
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expression
#' # If categorical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#'
#' # Note: fitExtractVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form <- ~ (1|sex) + (1|type) + score
#' info2 <- Biobase::pData(sample.ExpressionSet)
#' varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2 )
#' 
# # Parallel processing using multiple cores with reduced memory usage
# param = SnowParam(4, "SOCK", progressbar=TRUE)
# varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2, BPPARAM = param)
#'
#'
#' @export
#' @docType methods
#' @rdname fitExtractVarPartModel-method
#' @importFrom BiocParallel SerialParam bpiterate bplapply bpok
setGeneric("fitExtractVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
      standardGeneric("fitExtractVarPartModel")
)

# internal driver function
#' @importFrom methods is new
#' @importFrom lme4 lmer 
#' @importFrom RhpcBLASctl omp_set_num_threads
#' @importFrom progress progress_bar
#' @import foreach
.fitExtractVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE, colinearityCutoff=.999, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(),...){ 

	# convert to data.frame
	data = as.data.frame(data)

	# exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

	# only retain columns used in the formula
	data = data[, colnames(data) %in% unique(all.vars(formula)), drop=FALSE]

	# check dimensions of reponse and covariates
	if( ncol(exprObj) != nrow(data) ){		
		stop( "the number of samples in exprObj (i.e. cols) must be the same as in data (i.e rows)" )
	}

	# check if variables in formula has NA's
	hasNA = hasMissingData(formula, data)

	if( any(hasNA) ){
		warning(paste("Variables contain NA's:", paste(names(hasNA[hasNA]), collapse=', '), "\nSamples with missing data will be dropped.\n"), immediate.=TRUE)
	}

	if( ! is(exprObj, "sparseMatrix")){
		# check if values are NA
		countNA = sum(is.nan(exprObj)) + sum(!is.finite(exprObj))
		if( countNA > 0 ){
			stop("There are ", countNA, " NA/NaN/Inf values in exprObj\nMissing data is not allowed")
		}
		
		rv = apply( exprObj, 1, var)
	}else{
		# if exprObj is a sparseMatrix, this method will compute row-wise
		# variances with using additional memory
		rv = c()
		for( i in seq_len(nrow(exprObj)) ){
			rv[i] = var( exprObj[i,])
		}
	}
	if( any( rv == 0) ){
		idx = which(rv == 0)
		stop(paste("Response variable", idx[1], 'has a variance of 0'))
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
	if( .isDisconnected() ){
		stop("Cluster connection lost. Either stopCluster() was run too soon\n, or connection expired")
	}

	# add response (i.e. exprObj[,j] to formula
	form = paste( "responsePlaceholder$E", paste(as.character( formula), collapse=''))

	# run lmer() to see if the model has random effects
	# if yes run lmer() in the loop
	# else run lm()
	responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data, control=control,... ), error = function(e) e)

	# detect error when variable in formula does not exist
	if( inherits(possibleError, "error") ){
		err = grep("object '.*' not found", possibleError$message)
		if( length(err) > 0 ){
			stop("Variable in formula is not found: ", gsub("object '(.*)' not found", "\\1", possibleError$message) )
		}
	}

	if( !quiet) pb <- progress_bar$new(format = ":current/:total [:bar] :percent ETA::eta",
			total = nrow(exprObj), width= 60, clear=FALSE)

	if( ! .isMixedModelFormula( formula ) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

		testValue = calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )		

		timeStart = proc.time()

		varPart <- foreach(responsePlaceholder=exprIter(exprObj, weightsMatrix, useWeights), .packages=c("splines","lme4") ) %do% {

			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=responsePlaceholder$weights,na.action=stats::na.exclude,...)

			calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )		
		}

		modelType = "anova"

	}else{

		if( inherits(possibleError, "error") && grep('the fixed-effects model matrix is column rank deficient', possibleError$message) == 1 ){
			stop(paste(possibleError$message, "\n\nSuggestion: rescale fixed effect variables.\nThis will not change the variance fractions or p-values."))
		}  

		# fit first model to initialize other model fits
		# this make the other models converge faster
		responsePlaceholder = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control)
		timediff = proc.time() - timeStart

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff, immediate=TRUE )

		timeStart = proc.time()

		# Define function for parallel evaluation
		.eval_models = function(responsePlaceholder, data, form, REML, theta, control, na.action=stats::na.exclude,...){
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=responsePlaceholder$weights, control=control,na.action=na.action)

			calcVarPart( fit, showWarnings=showWarnings, colinearityCutoff=colinearityCutoff )
		}

		.eval_master = function( obj, data, form, REML, theta, control, na.action=stats::na.exclude,... ){

			# use only 1 OpenMP thread for linear algebra
			omp_set_num_threads(1)

			lapply(seq_len(nrow(obj$E)), function(j){
				.eval_models( list(E=obj$E[j,], weights=obj$weights[j,]), data, form, REML, theta, control, na.action,...)
			})
		}

		# Evaluate function
		####################

		it = iterBatch(exprObj, weightsMatrix, useWeights, n_chunks = 100, BPPARAM = BPPARAM)

		if( !quiet) message(paste0("Dividing work into ",attr(it, "n_chunks")," chunks..."))

		varPart <- bpiterate( it, .eval_master,
			data=data, form=form, REML=REML, theta=fitInit@theta, control=control,...,
			BPPARAM=BPPARAM)

		# if there is an error in evaluating fxn (usually in parallel backend)
		if( !bpok(list(varPart)) ){
      stop("Error evaluating fxn:\n\n", varPart)
		}
		# It can also return a list of errors, or a list where only some elements are errors
		if( !all(bpok(varPart)) ){
      first_error <- varPart[[which(!bpok(varPart))[1]]]
      stop("Error evaluating fxn:\n\n", first_error)
		}

		# If no errors, then it's safe to concatenate all the results together.
		varPart <- do.call(c, varPart)

		modelType = "linear mixed model"
	}

	if(!quiet) message("\nTotal:", paste(format((proc.time() - timeStart)[3], digits = 1, scientific = FALSE), "s"))

	varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
	colnames(varPartMat) <- names(varPart[[1]])
	rownames(varPartMat) <- rownames(exprObj)

	res <- new("varPartResults", varPartMat, type=modelType, method="Variance explained (%)")
	
	return( res )
}

# matrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,matrix-method
setMethod("fitExtractVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitExtractVarPartModel(exprObj, formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                     	 BPPARAM=BPPARAM, ...)
  }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod("fitExtractVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                     	 BPPARAM=BPPARAM, ...)
  }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod("fitExtractVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj$E), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights,  showWarnings=showWarnings, control=control, quiet=quiet,  
                         BPPARAM=BPPARAM, ...)
  }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod("fitExtractVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprs(exprObj)), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet, 
                     	 BPPARAM=BPPARAM,...)
  }
)

# sparseMatrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,sparseMatrix-method
setMethod("fitExtractVarPartModel", "sparseMatrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,  showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), quiet=FALSE, BPPARAM=SerialParam(), ...)
  {
    .fitExtractVarPartModel( exprObj, formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix,  showWarnings=showWarnings, control=control, quiet=quiet,
                     	 BPPARAM=BPPARAM, ...)
  }
)

