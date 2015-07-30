
# define classes
setClass("varParFrac")
setClass("VarParFitList")


# Iterator over genes
exprIter = function( exprObj, weights, useWeights = TRUE){
    xit <- icountn( nrow(exprObj) )
    nextEl <- function() {
    	j <- nextElem(xit)

    	if( useWeights && !is.null(weights) ){    		
			# scale weights to have mean of 1, otherwise it affects the residual variance too much
    		w = weights[j,] /  mean(weights[j,])
    	}else{
    		w = NULL		
		}

       	list(E = exprObj[j,], weights = w)
    }
    it <- list(nextElem = nextEl)
    class(it) <- c("abstractiter", "iter")
    it
}


#' Simulation dataset for examples
#'
#' A simulated dataset of gene expression and metadata
#'
#' \itemize{    
#'	\item geneCounts gene expression in the form of RNA-seq counts
#'	\item geneExpr gene expression on a continuous scale
#'  \item info metadata about the study design
#' }
#' @docType data
#' @keywords datasets
#' @usage data(varPartData)
#' @format A dataset of 100 samples and 200 genes
#' @name varPartData
NULL

#' Fit linear (mixed) model
#' 
#' Fit linear (mixed) model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)
#' @param data data.frame with columns corresponding to formula 
#' @param REML use restricted maximum likelihood to fit linear mixed model. default is FALSE.  Strongly discourage against changing this option
#' @param useWeights if TRUE, analysis uses heteroskadatistic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
# @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return
#' list() of where each entry is a model fit produced by lmer() or lm()
#' 
#' @import lme4 ggplot2 limma foreach reshape iterators dendextend doParallel Biobase methods
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
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # Can also be invoked with 
#' # cl <- makeCluster(2)
#' # registerDoParallel(cl)
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
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart = fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' # also sort columns
#' plotVarPart( sortCols( varPart ) )
#'
#' # Advanced: 
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results = fitVarPartModel( geneExpr, form, info )
#' 
#' # Step 2: extract variance fractions
#' varPart = extractVarPart( results )
#'
#' # Note: fitVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form = ~ (1|sex) + (1|type) + score
#' info2 = pData(sample.ExpressionSet)
#' results2 = fitVarPartModel( sample.ExpressionSet, form, info2 )
#' 
#' @export
#' @docType methods
#' @rdname fitVarPartModel-method
setGeneric("fitVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
      standardGeneric("fitVarPartModel")
)

# internal driver function
.fitVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL , ...){ 

	exprObj = as.matrix( exprObj )

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

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,... ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# if no intercept is specified, give warning
		if( length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if any coefficient is NA
		if( any(is.na(coef(fit))) ){
			warning("The variables specified in this model are redundant,\nso the design matrix is not full rank:\nThe results will not behave as expected and may be very wrong!!")
		}

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			lm( eval(parse(text=form)), data=data, weights=gene$weights,...)
		}

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=lmerControl(calc.derivs=FALSE))

		# if no intercept is specified, give warning
		if( length(which(colnames(fitInit@pp$X) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if a factor|character is modeled as a fixed effect
		fixedFactors = sapply( attr(terms(fitInit), "term.labels"), function(key) is.factor(fitInit@frame[[key]]) | is.character(fitInit@frame[[key]]))
		if( length(fixedFactors) > 0 && sum(fixedFactors) > 0 ){
			warning(paste("Categorical variables modeled as fixed effect:", paste(names(which(fixedFactors)), collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))
		}

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=gene$weights, control=lmerControl(calc.derivs=FALSE), start=fitInit@theta)
		}
	}

	# set name of each entry
	names(res) <- rownames( exprObj )
	class(res) <- c("VarParFitList", "list")

	return( res )
}

## matrix
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,matrix-method
setMethod("fitVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitVarPartModel(exprObj, formula, data, REML=REML, useWeights=useWeights, ...)
  }
)

# data.frame
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,data.frame-method
setMethod("fitVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitVarPartModel( as.matrix(exprObj), formula, data, REML=REML, useWeights=useWeights, ...)
  }
)

## EList 
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,EList-method
setMethod("fitVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitVarPartModel( as.matrix(exprObj$E), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights,...)
  }
)

## ExpressionSet
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,ExpressionSet-method
setMethod("fitVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitVarPartModel( as.matrix(exprs(exprObj)), formula, data, REML=REML, useWeights=useWeights, ...)
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
#' @param useWeights if TRUE, analysis uses heteroskadatistic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
# @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
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
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # Can also be invoked with 
#' # cl <- makeCluster(2)
#' # registerDoParallel(cl)
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
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart = fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#'
#' # plot a subset of variables
#' plotVarPart( varPart[,c("Individual", "Tissue")] )
#' 
#' # Note: fitExtractVarPartModel also accepts ExpressionSet
#' data(sample.ExpressionSet, package="Biobase")
#'
#' # ExpressionSet example
#' form = ~ (1|sex) + (1|type) + score
#' info2 = pData(sample.ExpressionSet)
#' varPart2 = fitExtractVarPartModel( sample.ExpressionSet, form, info2 )
#'
#' @export
#' @docType methods
#' @rdname fitExtractVarPartModel-method
setGeneric("fitExtractVarPartModel", signature="exprObj",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
      standardGeneric("fitExtractVarPartModel")
)

# internal driver function
.fitExtractVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, ...){ 

	exprObj = as.matrix( exprObj )

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

	# run lmer() to see if the model has random effects
	# if less run lmer() in the loop
	# else run lm()
	gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
	possibleError <- tryCatch( lmer( eval(parse(text=form)), data=data,... ), error = function(e) e)

	mesg <- "No random effects terms specified in formula"
	if( inherits(possibleError, "error") && identical(possibleError$message, mesg) ){

		# fit the model for testing
		fit <- lm( eval(parse(text=form)), data=data,...)

		# if no intercept is specified, give warning
		if( length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if any coefficient is NA
		if( any(is.na(coef(fit))) ){
			warning("The variables specified in this model are redundant,\nso the design matrix is not full rank:\nThe results will not behave as expected and may be very wrong!!")
		}

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=gene$weights,...)

			calcVarPart( fit )
		}

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=lmerControl(calc.derivs=FALSE))

		# if no intercept is specified, give warning
		if( length(which(colnames(fitInit@pp$X) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if a factor|character is modeled as a fixed effect
		fixedFactors = sapply( attr(terms(fitInit), "term.labels"), function(key) is.factor(fitInit@frame[[key]]) | is.character(fitInit@frame[[key]]))
		if( length(fixedFactors) > 0 && sum(fixedFactors) > 0 ){
			warning(paste("Categorical variables modeled as fixed effect:", paste(names(which(fixedFactors)), collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))
		}

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=gene$weights, control=lmerControl(calc.derivs=FALSE), start=fitInit@theta)

			calcVarPart( fit )
		}
	}

	varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
	colnames(varPartMat) <- names(varPart[[1]])
	rownames(varPartMat) <- rownames(exprObj)
	
	class(varPartMat) <- c("varParFrac", "data.frame")
	return( varPartMat )
}

# matrix
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,matrix-method
setMethod("fitExtractVarPartModel", "matrix",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitExtractVarPartModel(exprObj, formula, data,
                     REML=REML, useWeights=useWeights, ...)
  }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod("fitExtractVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj), formula, data,
                     REML=REML, useWeights=useWeights, ...)
  }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod("fitExtractVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj$E), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, ...)
  }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod("fitExtractVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprs(exprObj)), formula, data,
                     REML=REML, useWeights=useWeights, ...)
  }
)


#' Compute variance statistics
#' 
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model
#'
#' @param fit model fit from lm() or lmer()
#' @param ... additional arguments (not currently used)
#' 
#' @return
#' fraction of variance explained / ICC for each variable in the model
#' 
#' @examples
#' data(varPartData)
#'
#' # Linear mixed model
#' fit = lmer( geneExpr[1,] ~ (1|Tissue) + Age, info)
#' calcVarPart( fit )
#'
#' # Linear model
#' # Note that the two models produce slightly different results
#' # This is expected: they are different statistical estimates 
#' # of the same underlying value
#' fit = lm( geneExpr[1,] ~ Tissue + Age, info)
#' calcVarPart( fit )
#' 
#' @export
#' @docType methods
#' @rdname calcVarPart-method
setGeneric("calcVarPart", signature="fit",
  function(fit, ...)
      standardGeneric("calcVarPart")
)


# calcVarPart <- function( fit ){
# 	if( is.null(fit) ){
# 		stop("model fit is NULL.  There was an error fitting the regression model")
# 	}
# 	UseMethod("calcVarPart")
# }


#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lm-method
setMethod("calcVarPart", "lm",
function(fit, ...)
{
	a = anova(fit) 

	# get Sum of Squares
	varFrac = a[['Sum Sq']] / sum( a[['Sum Sq']] )
	# wrong : varFrac = a[['Mean Sq']] / sum( a[['Mean Sq']] )

	# get variable names
	names(varFrac) = rownames(a)

	return( varFrac )
}
)

#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lmerMod-method
setMethod("calcVarPart", "lmerMod",
function(fit, ...)
{
	# compute ICC, but don't divide by variances in the same class
	varComp <- lapply(VarCorr(fit), function(fit) attr(fit, "stddev")^2)

	# order variables by name
	# essential so that all models are ordered the same
	varComp = varComp[order(names(varComp))]

	# variance of fixed effects, if they exist
	# including fixed effects in the sum of variances makes a big difference, 
	#	especially when the fixed effect makes a big contribution
	# this approach minimized the difference between modelling
	# 	as a fixed or random effect
	# although, there is many be a substantial difference in estimates

	idx = which(colnames(fit@pp$X) != "(Intercept)")
	if( length(idx) > 0){

		# get predicted fixed effects
		fxeff = sapply( idx, function(i){
			fit@pp$X[,i] * fixef(fit)[i]
		})
		colnames(fxeff) = colnames(fit@pp$X)[idx]

		fixedVar = apply(fxeff, 2, var)

		for( i in 1:length(fixedVar) ){

			key = names(fixedVar)[i]
			varComp[[key]] = fixedVar[i]
			names(varComp[[key]]) = ''
		}

	}

	# get residual variance
	varComp$Residuals = attr(VarCorr(fit), 'sc')^2
	names(varComp$Residuals) = ''

	varPart = c()

	for( key in names(varComp) ){
		for( i in 1:length(varComp[[key]]) ){

			# get variance contributed by other variables
			varOtherVariables = varComp[-which(names(varComp) == key)]

			# demoninator: remove variables in this class (i.e. same key)
			# from the variables in the current class, only consider 1 variable
			denom =  sum(unlist(varOtherVariables)) + varComp[[key]][i]

			# compute fraction
			frac = varComp[[key]][i] / denom 

			# name variable based on two levels
			if( names(varComp[[key]])[i] %in% c("(Intercept)", '') ){
				names(frac) = key
			}else{
				names(frac) = paste( names(varComp[[key]])[i], key, sep='/')
			}

			# save result
			varPart = c(varPart, frac)
		}
	}

	return( varPart )
}
)


#' Extract variance statistics
#' 
#' Extract variance statistics from list of models fit with lm() or lmer()
#'
#' @param modelList list of lmer() model fits
#'  
#' @return 
#' data.frame of fraction of variance explained by each variable, after correcting for all others.
# @details
#' @examples
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart = fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#'
#' # Advanced: 
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results = fitVarPartModel( geneExpr, form, info )
#' 
#' # Step 2: extract variance fractions
#' varPart = extractVarPart( results )
#'
#' @export
extractVarPart <- function( modelList ){

	# for each model fit, get R^2 values
	entry <- 1
	varPart <- lapply( modelList, function( entry ) 
		calcVarPart( entry )
	)

	varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
	colnames(varPartMat) <- names(varPart[[1]])
	rownames(varPartMat) <- names(modelList)

	class(varPartMat) = c("varParFrac", "data.frame")
	
	return( varPartMat )
}

#' Sort variance parition statistics
#' 
#' Sort columns returned by extractVarPart() or fitExtractVarPartModel()
#'
#' @param x object returned by extractVarPart() or fitExtractVarPartModel()
#' @param FUN function giving summary statistic to sort by.  Defaults to median
#' @param decreasing  logical.  Should the sort be increasing or decreasing?  
#' @param ... other arguments to sort 
#'
#' @return
#' data.frame with columns sorted by mean value, with Residuals in last column
#' 
#' @export
#' @examples
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart = fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' # sort columns by median value
#' plotVarPart( sortCols( varPart ) )
sortCols <- function( x, FUN=median, decreasing = TRUE, ... ){

	# sort by column mean
	i = order(apply(x, 2, FUN), decreasing=decreasing)

	# apply sorting
	x = x[,i,drop=FALSE]

	# find column with Residuals
	i = which(colnames(x) == "Residuals" )

	if( length(i) > 0 ){
		# put residuals at right-most column
		res = cbind(x[,-i,drop=FALSE], Residuals=x[,"Residuals" ])
	}else{
		res = x[,,drop=FALSE]
	}

	return( res)
}




#' Default colors for ggplot
#' 
#' Return an array of n colors the same as the default used by ggplot2
#'
#' @param n number of colors
#'
#' @return
#' array of colors of length n
#' 
#' @examples
#' ggColorHue(4)
#' @export
ggColorHue <- function(n) {
  hues <- seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

#' Violin plot of variance fractions
#' 
#' Violin plot of variance fraction for each gene and each variable
#'
#' @param obj varParFrac object returned by fitExtractVarPart or extractVarPart
#' @param col vector of colors
#' @param label.angle angle of labels on x-axis
#' @param ylim limits of y-axis
#' @param main title of plot
#' 
#' @return
#' Makes violin plots of variance components model
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' varPart = fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#'
#' @export
plotVarPart <- function( obj, col, label.angle=20, ylim=c(0,100), main=""){

	# convert to data.frame
	obj = as.data.frame(obj)

	# if col is not specified
	if( missing( col ) ){
		# get c 
		col = ggColorHue(ncol(obj))
	}

	if( length(col) < ncol(obj) ){
		stop("Not enough colors specified by col")
	}

	# get gene name of each row
	obj$gene <- rownames(obj)

	# convert to data.frame for ggplot
	data <- melt(obj, id="gene")
	data$value <- data$value * 100

	# add to pass R CMD check
	variable <- 1
	value <- 1

	# violin plot
	fig = ggplot(data=data, aes(x=variable, y=value)) + geom_violin( scale="width", aes(fill = factor(variable))) + ylab("Variance explained (%)") + xlab('') + geom_boxplot(width=0.07, fill="grey", outlier.color='black')  + theme_bw() + scale_fill_manual(values=col) + theme(axis.text.x =
	               element_text(size  = 13,
	                            angle = label.angle,
	                            hjust = 1,
	                            vjust = 1)) + theme(legend.position="none") + ylim(ylim)

	if( main != ""){
		fig = fig + ggtitle( main ) +  theme(plot.title = element_text(lineheight=.8, face="bold"))
	}

	return( fig )
}


#' Residuals from model fit
#' 
#' Extract residuals for each gene from model fit with fitVarPartModel()
#'
#' @param object object produced by fitVarPartModel()
#' @param ... other arguments.
#' 
#' @return
#' Residuals extracted from model fits stored in object
#' 
#' @examples
#' # load library
#' # library(variancePartition)
#'
#' # optional step to run analysis in parallel on multicore machines
#' # Here, we used 4 threads
#' library(doParallel)
#' registerDoParallel(4)
#'
#' # load simulated data:
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#'
#' # Fit model
#' modelFit = fitVarPartModel( geneExpr, form, info )
#' 
#' # Extract residuals of model fit
#' res = residuals( modelFit )
#'
#' @export
#' @docType methods
setMethod("residuals", "VarParFitList",
  function(object, ...) {

		# extract residuals
		res <- lapply( object, function(fit)
			residuals( fit )
		)

		# convert to matrix, add names
		resMatrix <- data.frame(matrix(unlist(res), nrow=length(res), byrow=TRUE))
		colnames(resMatrix) <-  names(fitted.values(object[[1]]))
		rownames(resMatrix) <- names(object)
		
		return( resMatrix )
	}
)


# @method residuals VarParFitList


#' Colinearity score
#' 
#' Colinearity score for a regression model indicating if variables are too highly correlated to give meaningful results
#'
#' @param fit regression model fit from lm() or lmer()
#' 
#' @return
#' returns the colinearity score between 0 and 1, where a score > 0.99 means the degree of colinearity is too high
# reports the correlation matrix between coefficients estimates 
# for fixed effects
# the colinearity score is the maximum absolute correlation value of this matrix
#' @export
#' @examples
#' 
# load data 
#' # load library
#' # library(variancePartition)
#'
#' # load simulated data:
#' data(varPartData)
#
#' form = ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' res = fitVarPartModel( geneExpr[1:10,], form, info )
#'  
#' # evaluate the colinearity score on the first model fit
#' # this reports the correlation matrix between coefficients estimates
#' # for fixed effects
#' # the colinearity score is the maximum absolute correlation value
#' # If the colinearity score > .99 then the variance parition 
#' # estimates may be problematic
#' # In that case, a least one variable should be omitted
#' colinearityScore(res[[1]])
#' 
colinearityScore = function(fit){
	 # get correlation matrix
	 C = cov2cor(vcov(fit)); 

	 if( nrow(C) > 1 ){
	 	 # return largest off-diagonal absolute correlation
	 	score = max(abs(C[lower.tri(C)]))
	 }else{
	 	score = 0
	 }

	 attr( score, "vcor") = C
	 return(score)
}

#' plotStratifyBy
#'
#' Plot gene expression stratified by another variable
#'
#' @param geneExpr data.frame of gene expression values and another variable for each sample.  If there are multiple columns, the user can specify which one to use
#' @param xval name of column in geneExpr to be used along x-axis to stratify gene expression
#' @param yval name of column in geneExpr indicating gene expression
#' @param xlab label x-asis. Defaults to value of xval
#' @param ylab label y-asis. Defaults to value of yval
#' @param main main label
#' @param sortBy name of column in geneExpr to sort samples by.  Defaults to xval
#' @param colorBy name of column in geneExpr to color box plots.  Defaults to xval
#' @param sort if TRUE, sort boxplots by median value, else use default ordering
#' @param text plot text on the top left of the plot
#' @param text.y indicate position of the text on the y-axis as a fraction of the y-axis range
#' @param text.size size of text
#' @param pts.cex size of points
#' @param ylim specify range of y-axis
#'
#' @return
#' ggplot2 object
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#' 
#' # load simulated data:
#' data(varPartData)
#' 
#' # Create data.frame with expression and Tissue information for each sample
#' GE = data.frame( Expression = geneExpr[1,], Tissue = info$Tissue)
#' 
#' # Plot expression stratified by Tissue
#' plotStratifyBy( GE, "Tissue", "Expression")
#'
#' # Omit legend and color boxes grey
#' plotStratifyBy( GE, "Tissue", "Expression", colorBy = NULL)
#'
#' # Specify colors
#' col = c( B="green", A="red", C="yellow")
#' plotStratifyBy( GE, "Tissue", "Expression", colorBy=col, sort=FALSE)
#'
#' @export
plotStratifyBy = function( geneExpr, xval, yval, xlab=xval, ylab=yval, main=NULL, sortBy=xval, colorBy=xval, sort=TRUE, text=NULL, text.y=1, text.size=5, pts.cex=1, ylim=NULL ){

	geneExpr = data.frame( geneExpr )
    geneExpr = droplevels( geneExpr )

    sortBy = xval    

    # check that xval and yval exist in geneExpr
    if( !(xval %in% colnames(geneExpr)) ){
    	stop(paste("xval is not found in colnames(geneExpr): xval =", xval))
    }
    if( !(yval %in% colnames(geneExpr)) ){
    	stop(paste("yval is not found in colnames(geneExpr): yval =", yval))
    }

    # check that sortBy exist in geneExpr
    if( !(sortBy %in% colnames(geneExpr)) ){
    	stop(paste("sortBy is not found in colnames(geneExpr): sortBy =", sortBy))
    }

    geneExpr[[yval]] = as.numeric( geneExpr[[yval]] )

    xpos = 0.5 #text.x * nlevels(geneExpr[[xval]])
    ypos = text.y * (max(geneExpr[[yval]]) - min(geneExpr[[yval]])) + min(geneExpr[[yval]])

    if( sort ){    	
	    # sort categories by median expression
	    geneExpr[['reorder']] = reorder(geneExpr[[sortBy]],geneExpr[[yval]], FUN=median)
	    ord = "reorder"
    }else{
    	ord = xval
    }
   
    pOut = ggplot( geneExpr, aes_string(x=ord, y=yval)) + theme_bw() + theme(legend.position="none", axis.ticks.x=element_blank(), axis.text.x = element_blank(), legend.key = element_rect(color = 'white'), plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ylab(ylab) + xlab(xlab) 

    if(  is.null(colorBy) || is.na(colorBy) ){
        pOut = pOut + geom_boxplot(color="grey", fill="grey", outlier.color='black',outlier.shape = NA)
    }else{

    	# if colors are specified and all levels of xval are represented
    	if( sum(levels(geneExpr[[xval]]) %in% names(colorBy)) == nlevels(geneExpr[[xval]]) ){

    		i = match(levels(geneExpr[[ord]]), levels(geneExpr[[xval]]) )
 			pOut = pOut + geom_boxplot(aes_string(fill=xval), color=colorBy[i], outlier.color='black',outlier.shape = NA) + scale_fill_manual( values=array(colorBy))
 		}else{

	        # color boxes by colorBy variable in geneExpr
	        pOut = pOut + geom_boxplot( aes_string(color=colorBy, fill=colorBy), outlier.color='black', outlier.shape = NA)
	    }

	    # add legend
	    pOut = pOut + theme(legend.justification=c(1,0), legend.position=c(1,0))
    }

    # add median bar
    pOut = pOut + stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

    if( !is.null(ylim)){
    	pOut = pOut + ylim(ylim)
    }

    if( !is.null(main) ){
    	 pOut = pOut + ggtitle(main)
   	}

    if( !is.null(text) ){
        pOut = pOut + annotate("text", label = text, x = xpos, y=ypos, text.size = text.size, hjust=0)
    }

    pOut = pOut + geom_jitter(size=pts.cex,height=0, width=0, col="black")

    return( pOut )
}

