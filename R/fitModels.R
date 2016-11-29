
#' Fit linear (mixed) model
#' 
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)
#' @param data data.frame with columns corresponding to formula 
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  Strongly discourage against changing this option
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param fxn apply function to model fit for each gene.  Defaults to identify function so it returns the model fit itself
#' @param control control settings for lmer()
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return
#' list() of where each entry is a model fit produced by lmer() or lm()
#' 
#' @import gplots colorRamps lme4 pbkrtest ggplot2 limma foreach reshape2 iterators doParallel Biobase methods utils
# dendextend
#' @importFrom MASS ginv
# @importFrom RSpectra eigs_sym
#' @importFrom grDevices colorRampPalette hcl
#' @importFrom graphics abline axis hist image layout lines mtext par plot plot.new rect text title
#' @importFrom stats anova as.dendrogram as.dist cancor coef cov2cor density dist fitted.values hclust lm median model.matrix order.dendrogram quantile reorder residuals sd terms var vcov





# Commented out, for now
###### import colorRamps pbkrtest lmerTest ggplot2 limma foreach reshape2 iterators dendextend doParallel Biobase methods
###### 
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is ~ a + b + (1|c), then to model is 
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' If there are no random effects, so formula is ~ a + b + c, a 'standard' linear model is used:
#'
#' fit <- lm( exprObj[j,] ~ a + b + c, data=data)
#'
#' In both cases, useWeights=TRUE causes weightsMatrix[j,] to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using foreach/dopar to run loops in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lm/lmer.
#'
#' Since this function returns a list of each model fit, using this function is slower and uses more memory than fitExtractVarPartModel().
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
#' info2 <- pData(sample.ExpressionSet)
#' results2 <- fitVarPartModel( sample.ExpressionSet, form, info2 )
#' 
#' # stop cluster
#' stopCluster(cl)
#' @export
#' @docType methods
#' @rdname fitVarPartModel-method
setGeneric("fitVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ),...)
      standardGeneric("fitVarPartModel")
)

# internal driver function
.fitVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, colinearityCutoff=.999,control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...){ 

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

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
	form = paste( "gene$E", paste(as.character( formula), collapse=''))

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,...,control=control ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	method = ''
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff )

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=gene$weights,na.action=stats::na.exclude,...)

			# apply function
			fxn( fit )
		}

		method = "lm"

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control )

		timediff = proc.time() - timeStart

		# check size of stored objects
		objSize = object.size( fxn(fitInit) ) * nrow(exprObj)

		# total time = (time for 1 gene) * (# of genes) / 60 / (# of threads)
		showTime = timediff[3] * nrow(exprObj) / 60 / getDoParWorkers()

		cat("Projected memory usage: >", format(objSize, units = "auto"), "\n")

		if( showTime > .01 ){
			cat("Projected run time: ~", paste(format(showTime, digits=1), "min"), "\n")
		}

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff )

		# specify gene explicitly in data 
		# required for downstream processing with lmerTest
		data2 = data.frame(data, expr=gene$E, check.names=FALSE)
		form = paste( "expr", paste(as.character( formula), collapse=''))

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights), .packages="lme4" ) %dopar% {

			# modify data2 for this gene
			data2$expr = gene$E
 
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data2, ..., REML=REML, weights=gene$weights, start=fitInit@theta, control=control,na.action=stats::na.exclude)

			# apply function
			fxn( fit )
		}

		method = "lmer"
	}

	# set name of each entry
	names(res) <- rownames( exprObj )

 	new( "VarParFitList", res, method=method )
}

## matrix
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,matrix-method
setMethod("fitVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitVarPartModel(exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control,...)
  }
)

# data.frame
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,data.frame-method
setMethod("fitVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitVarPartModel( as.matrix(exprObj), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, ...)
  }
)

## EList 
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,EList-method
setMethod("fitVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitVarPartModel( as.matrix(exprObj$E), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, showWarnings=showWarnings, fxn=fxn, control=control,...)
  }
)

## ExpressionSet
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,ExpressionSet-method
setMethod("fitVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitVarPartModel( as.matrix(exprs(exprObj)), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, control=control, ...)
  }
)

#' Fit linear (mixed) model, report variance fractions
#' 
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Report fraction of variance attributable to each variable 
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)
#' @param data data.frame with columns corresponding to formula 
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  Strongly discourage against changing this option
#' @param useWeights if TRUE, analysis uses heteroskedastic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param adjust remove variation from specified variables from the denominator.  This computes the adjusted ICC with respect to the specified variables
#' @param adjustAll adjust for all variables.  This computes the adjusted ICC with respect to all variables.  This overrides the previous argument, so all variables are include in adjust.
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param control control settings for lmer()
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return
#' list() of where each entry is a model fit produced by lmer() or lm()
#'
#' @details 
#' A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression.  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is ~ a + b + (1|c), then to model is 
#'
#' fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)
#'
#' If there are no random effects, so formula is ~ a + b + c, a 'standard' linear model is used:
#'
#' fit <- lm( exprObj[j,] ~ a + b + c, data=data)
#'
#' In both cases, useWeights=TRUE causes weightsMatrix[j,] to be included as weights in the regression model.
#'
#' Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using foreach/dopar to run loops in parallel.  Parallel processing must be enabled before calling this function.  See below.
#' 
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lm/lmer.
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
#' info2 <- pData(sample.ExpressionSet)
#' varPart2 <- fitExtractVarPartModel( sample.ExpressionSet, form, info2 )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
#' @rdname fitExtractVarPartModel-method
setGeneric("fitExtractVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
      standardGeneric("fitExtractVarPartModel")
)

# internal driver function
.fitExtractVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, colinearityCutoff=.999, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ),...){ 

	exprObj = as.matrix( exprObj )
	formula = stats::as.formula( formula )

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

	# add response (i.e. exprObj[,j] to formula
	form = paste( "gene$E", paste(as.character( formula), collapse=''))

	# control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient") 

	# control = lme4::lmerControl(calc.derivs=FALSE, optCtrl=list(maxfun=1), check.rankX="stop.deficient" )

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data, control=control,... ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff )

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=gene$weights,na.action=stats::na.exclude,...)

			calcVarPart( fit, adjust, adjustAll, showWarnings, colinearityCutoff )
		}
		modelType = "anova"

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))

		timeStart = proc.time()
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=control)
		timediff = proc.time() - timeStart

		# total time = (time for 1 gene) * (# of genes) / 60 / (# of threads)
		showTime = timediff[3] * nrow(exprObj) / 60 / getDoParWorkers()

		if( showTime > .01 ){
			cat("Projected run time: ~", paste(format(showTime, digits=1), "min"), "\n")
		}

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff )

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights), .packages="lme4" ) %dopar% {
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=gene$weights, start=fitInit@theta, control=control,na.action=stats::na.exclude)

			calcVarPart( fit, adjust, adjustAll, showWarnings, colinearityCutoff )
		}

		modelType = "linear mixed model"
	}

	varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
	colnames(varPartMat) <- names(varPart[[1]])
	rownames(varPartMat) <- rownames(exprObj)

	# get list of variation removed from the denominator
	adjust = getAdjustVariables( colnames(varPartMat), adjust, adjustAll)
	if( is.null(adjust) ) adjust = NA

	if( any(!is.na(adjust)) ){
		method = "adjusted intra-class correlation"
	}else{
		method = "Variance explained (%)"
	}	

	res <- new("varPartResults", varPartMat, type=modelType, adjustedFor=array(adjust), method=method)
	
	return( res )
}

# matrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,matrix-method
setMethod("fitExtractVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitExtractVarPartModel(exprObj, formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, control=control, ...)
  }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod("fitExtractVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, control=control, ...)
  }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod("fitExtractVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj$E), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, control=control, ...)
  }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod("fitExtractVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ), ...)
  {
    .fitExtractVarPartModel( as.matrix(exprs(exprObj)), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, control=control,...)
  }
)


