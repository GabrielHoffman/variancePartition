# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai


# define classes
setClass("varParFrac")
setClass("VarParFitList", representation(method="character"), contains="list")


# Iterator over genes
exprIter = function( exprObj, weights, useWeights = TRUE, scale=TRUE){
    xit <- icountn( nrow(exprObj) )
    nextEl <- function() {
    	j <- nextElem(xit)

    	if( useWeights && !is.null(weights) ){    		
			# scale weights to have mean of 1, otherwise it affects the residual variance too much
    		if(scale){
    			w = weights[j,] /  mean(weights[j,])
    		}else{
    			w = weights[j,]
    		}
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




#' Collinearity score
#' 
#' Collinearity score for a regression model indicating if variables are too highly correlated to give meaningful results
#'
#' @param fit regression model fit from lm() or lmer()
#' 
#' @return
#' Returns the collinearity score between 0 and 1, where a score > 0.999 means the degree of collinearity is too high.  This function reports the correlation matrix between coefficient estimates for fixed effects.  The collinearity score is the maximum absolute correlation value of this matrix. Note that the values are the correlation between the parameter estimates, and not between the variables themselves.
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
#' # evaluate the collinearity score on the first model fit
#' # this reports the correlation matrix between coefficients estimates
#' # for fixed effects
#' # the collinearity score is the maximum absolute correlation value
#' # If the collinearity score > .999 then the variance partition 
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




# Check if model contains a random effect
# 
# Check if model contains a random effect
#
# @param formula model formula
# @param data data.frame
# 
.isMixedModelFormula = function(formula, data ){

    # don't throw an error if the LHS is missing
    control = lme4::lmerControl(check.formula.LHS = "ignore")

    possibleError <- tryCatch( lFormula( formula, data, control=control), error = function(e) e)

    mesg <- "No random effects terms specified in formula"
    result = inherits(possibleError, "error") && identical(possibleError$message, mesg)

    return( ! result )
}