# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai


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

# results: the variancePartition statistics
# type: indicate variance fractions or adjusted ICC
# adjustedFor: variables whose variance were removed from the denominator 
# type: lmm or anova
# setClass("varPartResults", representation(results = "data.frame", type = "character", adjustedFor="array", method="character"))

setClass("varPartResults", representation(type = "character", adjustedFor="array", method="character"), contains="data.frame")

# # @export
# setMethod("print", "varPartResults",
#   function(x, ...) {
#   	print( x )
#   }
# )
# setMethod("print", "varPartResults",
#   function(x, ...) {
#   	print( x )
#   }
# )



#' Simulation dataset for examples
#'
#' A simulated dataset of gene expression and metadata
#'
#' \itemize{    
#'	\item geneCounts: gene expression in the form of RNA-seq counts
#'	\item geneExpr: gene expression on a continuous scale
#'  \item info: metadata about the study design
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
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param fxn apply function to model fit for each gene.  Defaults to identify function so it returns the model fit itself
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return
#' list() of where each entry is a model fit produced by lmer() or lm()
#' 
#' @import lme4 ggplot2 limma foreach reshape iterators dendextend doParallel Biobase methods
#' @importFrom pbkrtest get_SigmaG
#' @importFrom MASS ginv
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
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omited by the underlying call to lm/lmer.
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
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
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
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL,showWarnings=TRUE,fxn=identity,...)
      standardGeneric("fitVarPartModel")
)

# internal driver function
.fitVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, colinearityCutoff=.999, ...){ 

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

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff )

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=gene$weights,...)

			# apply function
			fxn( fit )
		}

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML)
		# , control=lmerControl(calc.derivs=FALSE)

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff )

		res <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=gene$weights, start=fitInit@theta)
			# , control=lmerControl(calc.derivs=FALSE)

			# apply function
			fxn( fit )
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
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, ...)
  {
    .fitVarPartModel(exprObj, formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn,...)
  }
)

# data.frame
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,data.frame-method
setMethod("fitVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, ...)
  {
    .fitVarPartModel( as.matrix(exprObj), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, ...)
  }
)

## EList 
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,EList-method
setMethod("fitVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, showWarnings=TRUE,fxn=identity, ...)
  {
    .fitVarPartModel( as.matrix(exprObj$E), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, showWarnings=showWarnings, fxn=fxn,...)
  }
)

## ExpressionSet
#' @export
#' @rdname fitVarPartModel-method
#' @aliases fitVarPartModel,ExpressionSet-method
setMethod("fitVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, showWarnings=TRUE,fxn=identity, ...)
  {
    .fitVarPartModel( as.matrix(exprs(exprObj)), formula, data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings, fxn=fxn, ...)
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
#' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
#' @param adjust remove variation from specified variables from the denominator.  This computes the adjusted ICC with respect to the specified variables
#' @param adjustAll adjust for all variables.  This computes the adjusted ICC with respect to all variables.  This overrides the previous argument, so all variables are include in adjust.
#' @param showWarnings show warnings about model fit (default TRUE)
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
#' The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omited by the underlying call to lm/lmer.
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
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
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
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
      standardGeneric("fitExtractVarPartModel")
)

# internal driver function
.fitExtractVarPartModel <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, colinearityCutoff=.999,...){ 

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

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fit, showWarnings=showWarnings, colinearityCutoff )

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lm( eval(parse(text=form)), data=data, weights=gene$weights,...)

			calcVarPart( fit, adjust, adjustAll, showWarnings, colinearityCutoff )
		}
		modelType = "anova"

	}else{

		# fit first model to initialize other model fits
		# this make the other models converge faster
		gene = nextElem(exprIter(exprObj, weightsMatrix, useWeights))
		fitInit <- lmer( eval(parse(text=form)), data=data,..., REML=REML, control=lmerControl(calc.derivs=FALSE))

		# check that model fit is valid, and throw warning if not
		checkModelStatus( fitInit, showWarnings=showWarnings, colinearityCutoff )

		varPart <- foreach(gene=exprIter(exprObj, weightsMatrix, useWeights) ) %dopar% {
			# fit linear mixed model
			fit = lmer( eval(parse(text=form)), data=data, ..., REML=REML, weights=gene$weights, control=lmerControl(calc.derivs=FALSE), start=fitInit@theta)

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
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
  {
    .fitExtractVarPartModel(exprObj, formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, ...)
  }
)

# data.frame
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,data.frame-method
setMethod("fitExtractVarPartModel", "data.frame",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, ...)
  }
)

# EList
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,EList-method
setMethod("fitExtractVarPartModel", "EList",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprObj$E), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=exprObj$weights, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings, ...)
  }
)

# ExpressionSet
#' @export
#' @rdname fitExtractVarPartModel-method
#' @aliases fitExtractVarPartModel,ExpressionSet-method
setMethod("fitExtractVarPartModel", "ExpressionSet",
  function(exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
  {
    .fitExtractVarPartModel( as.matrix(exprs(exprObj)), formula, data,
                     REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, adjust=adjust, adjustAll=adjustAll, showWarnings=showWarnings,...)
  }
)



# Check that lm/lmer model is valid
# Throw warning if
#	1) Intercept is ommited
#	2) Any coefficient is NA
#	3) a categorical variable is modeled as a fixed effect
setGeneric("checkModelStatus", signature="fit",
  function( fit, showWarnings=TRUE, colinearityCutoff=.999 )
      standardGeneric("checkModelStatus")
)

setMethod("checkModelStatus", "lm",
  function( fit, showWarnings=TRUE, colinearityCutoff=.999 )
	{
		# if no intercept is specified, give warning
		if( showWarnings && length(which(names(coef(fit)) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if any coefficient is NA
		if( showWarnings && any(is.na(coef(fit))) ){
			stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
		}

		# check colinearity
		score = colinearityScore(fit)
		if( score > colinearityCutoff ){
			stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
		}
	}
)

setMethod("checkModelStatus", "lmerMod",
  function( fit, showWarnings=TRUE, colinearityCutoff=.999 )
	{
		# if no intercept is specified, give warning
		if( showWarnings && length(which(colnames(fit@pp$X) == "(Intercept)")) == 0 ){
			warning("No Intercept term was specified in the formula:\nThe results will not behave as expected and may be very wrong!!")
		}

		# if any coefficient is NA
		if( showWarnings && any(is.na(coef(fit))) ){
			stop("The variables specified in this model are redundant,\nso the design matrix is not full rank")
		}

		# check colinearity
		score = colinearityScore(fit)
		if( score > colinearityCutoff ){
			stop(paste("Colinear score =", format(score, digits=4), ">", colinearityCutoff,"\nCovariates in the formula are so strongly correlated that the\nparameter estimates from this model are not meaningful.\nDropping one or more of the covariates will fix this problem"))
		}

		# if a factor|character is modeled as a fixed effect
		fixedFactors = sapply( attr(terms(fit), "term.labels"), function(key) is.factor(fit@frame[[key]]) | is.character(fit@frame[[key]]))
		if( showWarnings && length(fixedFactors) > 0 && sum(fixedFactors) > 0 ){
			warning(paste("Categorical variables modeled as fixed effect:", paste(names(which(fixedFactors)), collapse=', '), "\nThe results will not behave as expected and may be very wrong!!"))
		}

		# Cannot model continuous variable as random effect
		variables = sapply(fit@frame, class)[-1]

		for( i in 1:length(variables) ){

			if( variables[i] == "numeric" && is.na(match( names(variables)[i], names(fixedFactors))) && showWarnings && names(variables)[i] != "(weights)"){
				stop(paste("Continuous variable cannot be modeled as a random effect:", names(variables)[i]))
			}
		}

		isMultipleVaryingCoefficientTerms( fit )
	}
)





#' Violin plot of variance fractions
#' 
#' Violin plot of variance fraction for each gene and each variable
#'
#' @param obj varParFrac object returned by fitExtractVarPart or extractVarPart
#' @param col vector of colors
#' @param label.angle angle of labels on x-axis
#' @param ylim limits of y-axis
#' @param main title of plot
#' @param ... additional arguments
#' 
#' @return
#' Makes violin plots of variance components model.  This function uses the graphics interface from ggplot2.  Warnings produced by this function usually ggplot2 warning that the window is too small.  
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
#' # Specify variables to consider
#' # Age is continuous so we model it as a fixed effect
#' # Individual and Tissue are both categorical, so we model them as random effects
#' form <- ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' plotVarPart( sortCols( varPart ) )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
#' @rdname plotVarPart-method
setGeneric("plotVarPart", signature="obj",
	function( obj, col=ggColorHue(ncol(obj)), label.angle=20, ylim=c(0,100), main="",...)
      standardGeneric("plotVarPart")
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,matrix-method
setMethod("plotVarPart", "matrix",
	function( obj, col=ggColorHue(ncol(obj)), label.angle=20, ylim=c(0,100), main="", ...){
 		.plotVarPart( obj, col, label.angle, ylim, main )
 	}
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,varPartResults-method
setMethod("plotVarPart", "data.frame",
	function( obj, col=ggColorHue(ncol(obj)), label.angle=20, ylim=c(0,100), main="",...){
 		.plotVarPart( obj, col, label.angle, ylim, main )
 	}
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,matrix-method
setMethod("plotVarPart", "varPartResults",
	function( obj, col=ggColorHue(ncol(obj)), label.angle=20, ylim=c(0,100), main="", ...){
		
 		.plotVarPart( data.frame(obj, check.names=FALSE), col, label.angle, ylim, main, ylab=obj@method)
 	}
)

# internal driver function
.plotVarPart <- function( obj, col=ggColorHue(ncol(obj)), label.angle=20, ylim=c(0,100), main="", ylab=''){

	# convert to data.frame
	obj = as.data.frame(obj, check.names=FALSE)

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
	fig = ggplot(data=data, aes(x=variable, y=value)) + 
		geom_violin( scale="width", aes(fill = factor(variable))) + 
		ylab(ylab) + xlab('') + ylim(ylim) + theme_bw() + 
		geom_boxplot(width=0.07, fill="grey", outlier.colour='black') + 
		scale_fill_manual(values=col) +
		theme(legend.position="none") +
		theme(axis.text.x = element_text(size  = 13,
	                            angle = label.angle,
	                            hjust = 1,
	                            vjust = 1)) 

	if( main != ""){
		fig = fig + ggtitle( main ) + theme(plot.title = element_text(lineheight=.8, face="bold"))
	}

	return( fig )
}


#' Compute variance statistics
#' 
#' Compute fraction of variation attributable to each variable in regression model.  Also interpretable as the intra-class correlation after correcting for all other variables in the model.
#'
#' @param fit model fit from lm() or lmer()
#' @param adjust remove variation from specified variables from the denominator.  This computes the adjusted ICC with respect to the specified variables
#' @param adjustAll adjust for all variables.  This computes the adjusted ICC with respect to all variables
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param ... additional arguments (not currently used)
#' 
#' @return
#' fraction of variance explained / ICC for each variable in the model
#' 
#' @examples
#' library(lme4)
#' data(varPartData)
#'
#' # Linear mixed model
#' fit <- lmer( geneExpr[1,] ~ (1|Tissue) + Age, info)
#' calcVarPart( fit )
#'
#' # Linear model
#' # Note that the two models produce slightly different results
#' # This is expected: they are different statistical estimates 
#' # of the same underlying value
#' fit <- lm( geneExpr[1,] ~ Tissue + Age, info)
#' calcVarPart( fit )
#' 
#' @export
#' @docType methods
#' @rdname calcVarPart-method
setGeneric("calcVarPart", signature="fit",
  function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
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
function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, ...)
{

	# check validity of model fit
	checkModelStatus( fit, showWarnings,...)

	# Get ANOVA
	a = anova(fit) 

	# get variables to remove from denominator
	# also, check variables
	adjust = getAdjustVariables( rownames(a), adjust, adjustAll)

	# get Sum of Squares
	if( is.null(adjust) ){
		varFrac = a[['Sum Sq']] / sum( a[['Sum Sq']] )
	}else{

		v = a[['Sum Sq']]
		names(v) = rownames(a)

		varFrac = c()
		for( i in 1:(length(v)-1) ){
			# total variance minus components to remove, but add back the current variable
			varFrac[i] = v[i] / get_denom( v, adjust, i)
		}
		# Residual variance
		varFrac[length(v)] = v[["Residuals"]] / sum(v)

	}

	# set names 
	names(varFrac) = rownames(a)

	return( varFrac )
}
)

get_denom = function( v, adjust, i){

	currentSet = setdiff(adjust, names(v)[i])

	if( length(currentSet) > 0){
		denom = sum(v) - sum(sapply(currentSet, function(x) v[x]))
	}else{
		denom = sum(v)
	}
	return(denom)
}

isMultipleVaryingCoefficientTerms = function( fit ){

	# get variance components values
	varComp = getVarianceComponents( fit )

	# get varying coefficient terms
	res = which(sapply(varComp, length) > 1)

	# if there are more than 1
	if( length(res) > 1){
		warning(paste("Cannot have more than one varying coefficient term:", paste(names(res), collapse=", "), "\nThe results will not behave as expected and may be very wrong!!"))
	}
}

#' @export
#' @rdname calcVarPart-method
#' @aliases calcVarPart,lmerMod-method
setMethod("calcVarPart", "lmerMod",
function(fit, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE,...)
{
	# check validity of model fit
	checkModelStatus( fit, showWarnings, ...)

	# get variance components values
	varComp = getVarianceComponents( fit )

	# get variables to remove from denominator
	# also, check variables
	adjust = getAdjustVariables( names(varComp), adjust, adjustAll)

	if( max(sapply(varComp, length)) > 1 && !is.null(adjust) ){
		stop("The adjust and adjustAll arguments are not currently supported for varying coefficient models")
	}

	variableLevels = list()
	for(key in colnames(fit@frame)){
		key2 = paste(paste(key, levels(fit@frame[[key]]), sep=''), collapse=',')
		variableLevels[[key2]] = key
	}

 	# Get variance terms for each variable
 	# for standard terms, return the variance
 	# for varying coefficient terms return weighted sum of variance
 	# 	weighted by the sample sizes corresponding to the subsets of the data
	get_total_variance = function(varComp){
		sapply(varComp, function(x){
		if( length(x) == 1){
			return( x )
		}else{			
			key = variableLevels[[paste(names(x), collapse=',')]]

			# get total variance from varying coefficient model
			weights = (table(fit@frame[[key]]) / nrow(fit@frame))
			x %*% weights / sum(weights)
		}
		})
	}

	varPart = c()

	for( key in names(varComp) ){
		for( i in 1:length(varComp[[key]]) ){

			# get variance contributed by other variables
			varOtherVariables = varComp[-which(names(varComp) == key)]

			# remove this variance from the denominator
			currentSet = setdiff(adjust, key)

			if( length(currentSet) > 0 && key != "Residuals"){
				adjustVariance = sum(sapply(currentSet, function(x) varComp[[x]]))
			}else{
				adjustVariance = 0
			}

			# demoninator: remove variables in this class (i.e. same key)
			# from the variables in the current class, only consider 1 variable
			totVar = tryCatch({
			     get_total_variance(varOtherVariables)
			}, error = function(e) {
			    stop("Problem with varying coefficient model in formula: should have form (A+0|B)")
			}, finally = {
			})
			denom = sum(totVar) + varComp[[key]][i] - adjustVariance

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


#' Extract variance terms
#' 
#' Extract variance terms from a model fit with lm() or lmer()
#'
#' @param fit list of lmer() model fits
#'  
#' @return 
#'  variance explained by each variable
# @details
#' @examples
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
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' modelList <- fitVarPartModel( geneExpr, form, info )
#' 
#' fit <- modelList[[1]]
#' getVarianceComponents( fit )
#' 
#' # stop cluster
#' stopCluster(cl)
#' 
#' @export
getVarianceComponents = function( fit ){
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

	return( varComp )
}

# Construct list of variable names to remove from the denominator
getAdjustVariables = function( variables, adjust, adjustAll){

	if( adjustAll ){
		adjust = variables
	}

	# Check that variables in adjust array are present
	idx = match(adjust, variables)

	if( any(is.na(idx)) ){
		stop(paste("The following variables cannot be used in the 'adjust' argument\nbecause they are not present in the model fit:", paste(adjust[which(is.na(idx))], collapse='; ')))
	}

		# remove Residuals from the adjust list
	if( "Residuals" %in% adjust){
		adjust = adjust[-match("Residuals", adjust)]
	}

	return( adjust )
}


#' Extract variance statistics
#' 
#' Extract variance statistics from list of models fit with lm() or lmer()
#'
#' @param modelList list of lmer() model fits
#' @param adjust remove variation from specified variables from the denominator.  This computes the adjusted ICC with respect to the specified variables
#' @param adjustAll adjust for all variables.  This computes the adjusted ICC with respect to all variables. This overrides the previous argument, so all variables are include in adjust.
#' @param showWarnings show warnings about model fit (default TRUE)
#' @param ... other arguments
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
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
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
#' # Advanced: 
#' # Fit model and extract variance in two separate steps
#' # Step 1: fit model for each gene, store model fit for each gene in a list
#' results <- fitVarPartModel( geneExpr, form, info )
#' 
#' # Step 2: extract variance fractions
#' varPart <- extractVarPart( results )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
extractVarPart <- function( modelList, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE,... ){

	# get results from first model to enumerate all variables present
	singleResult = calcVarPart( modelList[[1]], adjust, adjustAll, showWarnings=showWarnings,... )

	# get variables to remove from denominator
	# also, check variables
	adjust = getAdjustVariables( names(singleResult), adjust, adjustAll)

	# for each model fit, get R^2 values
	entry <- 1
	varPart <- lapply( modelList, function( entry ) 
		calcVarPart( entry, adjust, adjustAll, showWarnings=showWarnings,... )
	)

	varPartMat <- data.frame(matrix(unlist(varPart), nrow=length(varPart), byrow=TRUE))
	colnames(varPartMat) <- names(varPart[[1]])
	rownames(varPartMat) <- names(modelList)

	modelType = ifelse(class(modelList[[1]])[1] == "lm", "anova", "linear mixed model")

	if( is.null(adjust) ) adjust = NA

	if( any(!is.na(adjust)) ){
		method = "adjusted intra-class correlation"
	}else{
		method = "Variance explained (%)"
	}	

	res <- new("varPartResults", varPartMat, type=modelType, adjustedFor=array(adjust), method=method)
	
	return( res )
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
#' @details
#' If model is fit with missing data, residuals returns NA for entries that were 
#' missing in the original data
#'
#' @examples
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
#' # Fit model
#' modelFit <- fitVarPartModel( geneExpr, form, info )
#' 
#' # Extract residuals of model fit
#' res <- residuals( modelFit )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
#' @aliases residuals
setMethod("residuals", "VarParFitList",
  function(object, ...) {

		# extract residuals
		res <- lapply( object, function(fit)
			residuals( fit )
		)

		# identify which samples were omitted
		excludeList <- sapply( object, function(fit)
			getOmitted( fit )
		)

		# get total number of samples
		n_samples = sapply(res, length) + sapply(excludeList, length)

		if( max(n_samples) - min(n_samples)  > 0 ){
			stop("Samples were dropped from model fit.  Either expression data or metadata contained NA values")
		}

		# create matrix of total size, including omitted samples
		resMatrix <- matrix(NA, nrow=length(res), ncol=n_samples[1])

		# fill non-missing entries with residuals
		for( j in 1:nrow(resMatrix)){
			excl = excludeList[[j]]

			if( is.null(excl) ){
				resMatrix[j,] = res[[j]] 
			}else{				
				resMatrix[j,-excl] = res[[j]] 
			}
		}

		# gene names along rows
		rownames(resMatrix) <- names(object)

		# get columns names, filling NA values from excludeList
		sampleNames = rep(NA, ncol(resMatrix))
		excl = excludeList[[1]]

		if( is.null(excl) ){
			sampleNames = names(fitted.values(object[[1]])) 
		}else{				
			sampleNames[-excl] = names(fitted.values(object[[1]]))
			sampleNames[excl] = names(excl)
		}
		colnames(resMatrix) = sampleNames
		
		return( as.matrix( resMatrix ) )
	}
)

setGeneric("getOmitted", signature="fit",
  function(fit, ...)
      standardGeneric("getOmitted")
)

setMethod("getOmitted", "lm",
  function(fit, ...)
  {
    fit$na.action
  }
)

setMethod("getOmitted", "lmerMod",
  function(fit, ...)
  {
    attr(fit@frame, "na.action")
  }
)




#' Colinearity score
#' 
#' Colinearity score for a regression model indicating if variables are too highly correlated to give meaningful results
#'
#' @param fit regression model fit from lm() or lmer()
#' 
#' @return
#' Returns the colinearity score between 0 and 1, where a score > 0.999 means the degree of colinearity is too high.  This function reports the correlation matrix between coefficient estimates for fixed effects.  The colinearity score is the maximum absolute correlation value of this matrix. Note that the values are the correlation between the parameter estimates, and not between the variables themselves.
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
#' form <- ~ Age + (1|Individual) + (1|Tissue) 
#' 
#' res <- fitVarPartModel( geneExpr[1:10,], form, info )
#'  
#' # evaluate the colinearity score on the first model fit
#' # this reports the correlation matrix between coefficients estimates
#' # for fixed effects
#' # the colinearity score is the maximum absolute correlation value
#' # If the colinearity score > .999 then the variance parition 
#' # estimates may be problematic
#' # In that case, a least one variable should be omitted
#' colinearityScore(res[[1]])
#' 
colinearityScore = function(fit){
	 # get correlation matrix
	 V = vcov(fit)
	 if( any(is.na(vcov(fit))) ){
	 	C = NA
	 }else{
	 	C = cov2cor(as.matrix(V))
	 } 

	 if( any(is.na(vcov(fit))) ){
	 	score = 1
	 }else if( nrow(C) > 1 ){
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
#' @param legend show legend
#' @param x.labels show x axis labels
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
plotStratifyBy = function( geneExpr, xval, yval, xlab=xval, ylab=yval, main=NULL, sortBy=xval, colorBy=xval, sort=TRUE, text=NULL, text.y=1, text.size=5, pts.cex=1, ylim=NULL, legend=TRUE, x.labels=FALSE ){

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
   
    pOut = ggplot( geneExpr, aes_string(x=ord, y=yval)) + theme_bw() + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ylab(ylab) + xlab(xlab) 

    if(  is.null(colorBy) || is.na(colorBy) ){
        pOut = pOut + geom_boxplot(color="grey", fill="grey", outlier.colour='black',outlier.shape = 20)
    }else{

    	# if colors are specified and all levels of xval are represented
    	if( sum(levels(geneExpr[[xval]]) %in% names(colorBy)) == nlevels(geneExpr[[xval]]) ){

    		i = match(levels(geneExpr[[ord]]), levels(geneExpr[[xval]]) )
 			pOut = pOut + geom_boxplot(aes_string(fill=xval), color=colorBy[i], outlier.colour='black',outlier.shape = 20) + scale_fill_manual( values=array(colorBy))
 		}else{

	        # color boxes by colorBy variable in geneExpr
	        pOut = pOut + geom_boxplot( aes_string(color=colorBy, fill=colorBy), outlier.colour='black', outlier.shape = 20)
	    }

	    # add legend
	    if( legend ){
	    	pOut = pOut + theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key = element_rect(color = 'white'), axis.text.x=element_text(angle=30))
	    }else{
	    	pOut = pOut + theme(legend.position="none", axis.text.x=element_text(angle=30))
	    }
    }

    if( ! x.labels ){
		pOut = pOut + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
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
        pOut = pOut + annotate("text", label = text, x = xpos, y=ypos, size = text.size, hjust=0)
    }

    #pOut = pOut + geom_jitter(size=pts.cex,height=0, width=0, col="black")

    return( pOut )
}



#' Sort variance parition statistics
#' 
#' Sort columns returned by extractVarPart() or fitExtractVarPartModel()
#'
#' @param x object returned by extractVarPart() or fitExtractVarPartModel()
#' @param FUN function giving summary statistic to sort by.  Defaults to median
#' @param decreasing  logical.  Should the sorting be increasing or decreasing?  
#' @param ... other arguments to sort 
#'
#' @return
#' data.frame with columns sorted by mean value, with Residuals in last column
#' 
#' @examples
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
#' # Step 1: fit linear mixed model on gene expresson
#' # If categoritical variables are specified, a linear mixed model is used
#' # If all variables are modeled as continuous, a linear model is used
#' # each entry in results is a regression model fit on a single gene
#' # Step 2: extract variance fractions from each model fit
#' # for each gene, returns fraction of variation attributable to each variable 
#' # Interpretation: the variance explained by each variable
#' # after correction for all other variables
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'  
#' # violin plot of contribution of each variable to total variance
#' # sort columns by median value
#' plotVarPart( sortCols( varPart ) )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
#' @docType methods
#' @rdname sortCols-method
setGeneric("sortCols", signature="x",
	function( x, FUN=median, decreasing = TRUE, ... )
      standardGeneric("sortCols")
)

#' @export
#' @rdname sortCols-method
#' @aliases sortCols,matrix-method
setMethod("sortCols", "matrix",
	function( x, FUN=median, decreasing = TRUE, ... ){
 		.sortCols(x, FUN, decreasing, ... )
 	}
)

#' @export
#' @rdname sortCols-method
#' @aliases sortCols,data.frame-method
setMethod("sortCols", "data.frame",
	function( x, FUN=median, decreasing = TRUE, ... ){
 		.sortCols(x, FUN, decreasing, ... )
 	}
)

#' @export
#' @rdname sortCols-method
#' @aliases sortCols,varPartResults-method
setMethod("sortCols", "varPartResults",
	function( x, FUN=median, decreasing = TRUE, ... ){
 		res = .sortCols( data.frame(x, check.names=FALSE), FUN, decreasing, ... )

 		vp = new( "varPartResults", res, type=x@type, adjustedFor=x@adjustedFor, method=x@method)

 		return( vp )
 	}
)

# internal driver function
.sortCols = function( x, FUN=median, decreasing = TRUE, ... ){

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




# #' Effective sample size
# #' 
# #' Compute effective sample size based on correlation structure in linear mixed model
# #'
# #' @param fit model fit from lmer()
# #' @param method "full" uses the full correlation structure of the model. The "approximate" method makes the simplifying assuption that the study has a mean of m samples in each of k groups, and computes m based on the study design.  When the study design is evenly balanced (i.e. the assumption is met), this gives the same results as the "full" method.  
# #' 
# #' @return
# #' effective sample size for each random effect in the model
# #' 
# #' @details
# #'
# #' Effective sample size calculations are based on:

# #' Liu, G., and Liang, K. Y. (1997). Sample size calculations for studies with correlated observations. Biometrics, 53(3), 937-47.
# #'
# #' "full" method: if V_x = var(Y;x) is the variance-covariance matrix of Y, the respose, based on the covariate x, then the effective sample size corresponding to this covariate is \\Sigma_\{i,j\} (V_x^\{-1\})_\{i,j\}.  In R notation, this is: sum(solve(V_x)).
# #'
# #' "approximate" method: Letting m be the mean number of samples per group, k be the number of groups, and rho be the intraclass correlation, the effective sample size is m*k / (1+rho*(m-1))
# #'
# #' Note that these values are equal when there are exactly m samples in each group.  If m is only an average then this an approximation.
# #'
# #' @examples
# #' library(lme4)
# #' data(varPartData)
# #'
# #' # Linear mixed model
# #' fit <- lmer( geneExpr[1,] ~ (1|Individual) + (1|Tissue) + Age, info)
# #'
# #' # Effective sample size
# #' ESS( fit )
# #' 
# #' @export
# #' @docType methods
# #' @rdname ESS-method
# setGeneric("ESS", signature="fit",
#   function(fit, method="full")
#       standardGeneric("ESS")
# )

# #' @export
# #' @rdname ESS-method
# #' @aliases ESS,lmerMod-method
# setMethod("ESS", "lmerMod",
# 	function( fit, method="full" ){

# 		if( !(method %in% c("full", "approximate")) ){
# 			stop(paste("method is not valid:", method))
# 		}

# 		# get correlation terms
# 		vp = calcVarPart( fit )

# 		n_eff = c()

# 		if( method == 'full'){

# 			# get structure of study design
# 			sigG = get_SigmaG( fit )

# 			ids = names(coef(fit))
# 			for( key in ids){
# 				i = which( key == ids)
# 				C = as.matrix(sigG$G[[i]]) * vp[[key]]
# 				diag(C) = 1 # set diagonals to 1
# 				n_eff[i] = sum(ginv(as.matrix(C)))
# 			}
# 			names(n_eff) = ids
# 		}else{

# 			ids = names(coef(fit))
# 			for( key in ids){
# 				i = which( key == ids)
# 				rho = vp[[key]]
# 				k = nlevels(fit@frame[[key]])
# 				m = nrow(fit@frame) / k 
# 				n_eff[i] = m*k / (1+rho*(m-1))
# 			}
# 			names(n_eff) = ids
# 		}
# 		return( n_eff )
# 	}
# )



#' Bar plot of variance fractions
#'
#' Bar plot of variance fractions for a subset of genes
#'
#' @param varPart object returned by extractVarPart() or fitExtractVarPartModel()
#' @param col color of bars for each variable
#'
#' @return Returns ggplot2 barplot
#' @examples
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
#' form <- ~ Age + (1|Individual) + (1|Tissue)
#'
#' # Fit model
#' varPart <- fitExtractVarPartModel( geneExpr, form, info )
#'
#' # Bar plot for a subset of genes showing variance fractions
#' plotPercentBars( varPart[1:5,] )
#'
#' # Move the legend to the top
#' plotPercentBars( varPart[1:5,] ) + theme(legend.position="top") 
#' 
#' # stop cluster
#' stopCluster(cl)
#' @export
plotPercentBars = function( varPart, col = ggColorHue(ncol(varPart)) ){

	if( !is.matrix(varPart) && !is.data.frame(varPart)){
		stop("Argument must be a matrix or data.frame")
	} 

	if( length(col) < ncol(varPart) ){
		stop("Number of colors is less than number of variables")
	}

	# convert matrix to tall data.frame
	df = melt(varPart, id.vars=NULL)

	# assign gene names
	df$gene = rep(rownames(varPart), ncol(varPart))

	# convert gene names to factors sorted so first gene is 
	# plotted on top
	df$gene = factor(df$gene, rev(rownames(varPart)))

	# convert values from [0-1] to [0-100]
	df$value = 100*df$value

	# Initialize variables to satisfy R CMD check
	gene = value = variable = 0

	# plot
	fig = ggplot(df, aes(x = gene, y = value, fill = variable)) + 
		geom_bar(stat = "identity") + theme_bw() + 
		theme(panel.grid.major = element_blank(), 
		panel.grid.minor = element_blank()) + coord_flip() + 
		ylab("Variance explained (%)") + xlab("")

	fig = fig + theme(axis.line = element_line(colour = "white"),
		axis.line.x = element_line(colour = "black"),
		panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(),
		panel.border = element_blank(),
		panel.background = element_blank(), 
		axis.ticks.y = element_blank(), 
		legend.key = element_blank()) +
		guides(fill=guide_legend(title=NULL)) +
		scale_fill_manual( values = col) + scale_y_continuous(expand=c(0,0.03))

	fig		
}


