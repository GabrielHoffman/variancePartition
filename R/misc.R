# (C) 2016 Gabriel E. Hoffman
# Icahn School of Medicine at Mount Sinai


#' Class varParFrac
#'
#' Class \code{varParFrac} 
#'
#' @name varParFrac-class
#' @rdname varParFrac-class
#' @exportClass varParFrac
setClass("varParFrac")

#' Class VarParFitList
#'
#' Class \code{VarParFitList} 
#'
#' @name VarParFitList-class
#' @rdname VarParFitList-class
#' @exportClass VarParFitList
setClass("VarParFitList", representation(method="character"), contains="list")

# requared to that iterator return NULL after the last element
icount2 = function (count){
    if (missing(count))
        count <- NULL
    else if (!is.numeric(count) || length(count) != 1)
        stop("count must be a numeric value")
    i <- 0L
    nextEl <- function(){
    	if( is.null(i) )    		
        	(i <<- NULL)
        else if (is.null(count) || i < count)
            (i <<- i + 1L)
        else 
        	(i <<- NULL)
    }
    it <- list(nextElem = nextEl)
    class(it) <- c("abstractiter", "iter")
    it
}


# Iterator over genes
exprIter = function( exprObj, weights, useWeights = TRUE, scale=TRUE, iterCount = "icount"){

	n_features = nrow(exprObj)

	if( iterCount == 'icount2'){
		xit <- icount2( n_features )
	}else{
		xit <- icount( n_features )
	}

    nextEl <- function() {
    	j <- nextElem(xit)

    	if( is.null(j) || j > n_features){
    		res = NULL
    	}else{
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

	       	res = list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
	     }
	     res
    }
    it <- list(nextElem = nextEl)
    class(it) <- c("abstractiter", "iter")
    it
}


# exprIterOrig = function( exprObj, weights, useWeights = TRUE, scale=TRUE){

# 	n_features = nrow(exprObj)
# 	xit <- icountn( n_features )

#     nextEl <- function() {
#     	j <- nextElem(xit)

#     	if( useWeights && !is.null(weights) ){    		
# 			# scale weights to have mean of 1, otherwise it affects the residual variance too much
#     		if(scale){
#     			w = weights[j,] /  mean(weights[j,])
#     		}else{
#     			w = weights[j,]
#     		}
#     	}else{
#     		w = NULL		
# 		}

#        	list(E = exprObj[j,], weights = w, n_iter = j, max_iter = n_features)
#     }
#     it <- list(nextElem = nextEl)
#     class(it) <- c("abstractiter", "iter")
#     it
# }


iterBatch <- function(exprObj, weights, useWeights = TRUE, scale=TRUE, n_chunks = nrow(exprObj) / 500 ) {

	# if there are fewer rows than chuncks, set n_chunks=1
	n_chunks = ifelse( nrow(exprObj) >= n_chunks, n_chunks, 1)

	# specify chunks
    idx <- parallel::splitIndices(nrow(exprObj), min(nrow(exprObj), n_chunks))
    i <- 0L

    f = function() {
    	if (i == length(idx)){
            return(NULL)
        }
        i <<- i + 1L
        E = exprObj[ idx[[i]],, drop = FALSE ]

        if( useWeights && !is.null(weights) ){    		
			# scale weights to have mean of 1, otherwise it affects the residual variance too much
    		if(scale){
    			# w = weights[j,] /  mean(weights[j,])
    			w = weights[idx[[i]],,drop=FALSE]
    			# for each row, devide by row mean
    			w = w / rowMeans(w)
    		}else{
    			# w = weights[j,]
    			w = weights[idx[[i]],,drop=FALSE]
    		}
    	}else{
    		w = matrix(1, nrow(E), ncol(E))		
		}
        
       	list(E = E, weights = w )
    }

    # get number of chunks
    attr( f, "n_chunks") = length(idx)
    f
}




# results: the variancePartition statistics
# type: indicate variance fractions or adjusted ICC
# adjustedFor: variables whose variance were removed from the denominator 
# type: lmm or anova
# setClass("varPartResults", representation(results = "data.frame", type = "character", adjustedFor="array", method="character"))

#' Class varPartResults
#'
#' Class \code{varPartResults} 
#'
#' @name varPartResults-class
#' @rdname varPartResults-class
#' @exportClass varPartResults
setClass("varPartResults", representation(type = "character", method="character"), contains="data.frame")

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

#' Convert to matrix
#' 
#' Convert varPartResults to matrix
#'
#' @param x varPartResults
#' @param ... other arguments.
#' 
#' @return
#' matrix
#' @examples
#' # load library
#' # library(variancePartition)
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
#' varPart <- fitExtractVarPartModel( geneExpr[1:5,], form, info )
#' 
#' # convert to matrix
#' as.matrix(varPart)
#'
#' @export
# @docType methods
#' @aliases as.matrix
setMethod("as.matrix", "varPartResults",
function(x, ...)
{
	df = as.data.frame( x@.Data )
	colnames(df) = names(x)
	rownames(df) = x@row.names
	# colnames(df) = colnames(x)
	# rownames(df) = rownames(x)

	return( as.matrix(df) )
})

#' Convert to data.frame
#' 
#' Convert varPartResults to data.frame
#'
#' @param x varPartResults
#' @param row.names pass thru to generic
#' @param optional pass thru to generic
#' @param ... other arguments.
#' 
#' @return
#' data.frame
#' @examples
#' # load library
#' # library(variancePartition)
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
#' varPart <- fitExtractVarPartModel( geneExpr[1:5,], form, info )
#' 
#' # convert to matrix
#' as.data.frame(varPart)
#'
#' @export
# @docType methods
#' @aliases as.data.frame
setMethod("as.data.frame", "varPartResults",
function(x, row.names = NULL, optional = FALSE,...)
{
	df = as.data.frame( x@.Data )
	colnames(df) = names(x)
	rownames(df) = x@row.names
	# colnames(df) = colnames(x)
	# rownames(df) = rownames(x)

	return( df )
})

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
#' # Intialize parallel backend with 4 cores
#' library(BiocParallel)
#' register(SnowParam(4))
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
# # stop cluster
# stopCluster(cl)
#'
#' @export 
# @docType methods
#' @aliases residuals,VarParFitList-method
setMethod("residuals", "VarParFitList",
  function(object, ...) {

		# extract residuals
		res <- lapply( object, function(fit)
			residuals( fit )
		)

		resMatrix <- matrix(NA, nrow=length(res), ncol=length(res[[1]]))
		colnames(resMatrix) = names(res[[1]])

		# fill non-missing entries with residuals
		for( j in 1:nrow(resMatrix)){
			resMatrix[j,] = res[[j]] 
		}

		# get total number of samples
		n_samples = sapply(res, length)

		if( max(n_samples) - min(n_samples) > 0 ){
			stop("Samples were dropped from model fit.  Either expression data or metadata contained NA values")
		}

		# identify which samples were omitted
		# excludeList <- sapply( object, function(fit)
		# 	getOmitted( fit )
		# )

		# # get total number of samples
		# n_samples = sapply(res, length)# + sapply(excludeList, length)

		# if( max(n_samples) - min(n_samples) > 0 ){
		# 	stop("Samples were dropped from model fit.  Either expression data or metadata contained NA values")
		# }

		# # create matrix of total size, including omitted samples
		# resMatrix <- matrix(NA, nrow=length(res), ncol=n_samples[1])

		# # fill non-missing entries with residuals
		# for( j in 1:nrow(resMatrix)){
		# 	excl = excludeList[[j]]

		# 	if( is.null(excl) ){
		# 		resMatrix[j,] = res[[j]] 
		# 	}else{				
		# 		resMatrix[j,-excl] = res[[j]] 
		# 	}
		# }

		# # gene names along rows
		# rownames(resMatrix) <- names(object)

		# # get columns names, filling NA values from excludeList
		# sampleNames = rep(NA, ncol(resMatrix))
		# excl = excludeList[[1]]

		# if( is.null(excl) ){
		# 	sampleNames = names(fitted.values(object[[1]])) 
		# }else{				
		# 	sampleNames[-excl] = names(fitted.values(object[[1]]))
		# 	sampleNames[excl] = names(excl)
		# }
		# colnames(resMatrix) = sampleNames
		
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
	 if( any(is.na(vcov(fit))) || nrow(V) == 0 ){
	 	C = NA
	 }else{
		 C = cov2cor(as.matrix(V))		 
	 } 

	 score = 0
	 if( any(is.na(vcov(fit))) ){
	 	score = 1
	 }else if( length(C) > 1 &&  nrow(C) > 1 ){
	 	 # return largest off-diagonal absolute correlation
	 	score = max(abs(C[lower.tri(C)]))
	 }

	 attr( score, "vcor") = C
	 return(score)
}




#' Check if model contains a random effect
#' 
#' Check if model contains a random effect
#'
#' @param formula model formula
#' 
#' @importFrom stats as.formula
# @importFrom lme4 lmerControl 
#' @importFrom lme4 findbars 
#' @keywords internal
.isMixedModelFormula = function(formula ){

	! is.null( findbars( as.formula( formula ) ) )
}



# Get process IDs
# 
# Get process IDs for parallel processing
#
# @param vmax minumum number of processes
#' @keywords internal
.get_pids = function( vmax = 2){

	# use do instead of dopar, since only need single pid
	foreach(i = 1:vmax, .combine=c) %do% { 
		Sys.getpid()
	}
}



# Try evaluating foreach dopar loop
#
# if it fails, connection is disconnected, so return TRUE, else FALSE
#
# This detects the error thrown by foreach dopar, after stopCluster() is run 
# By no new clusters is registered
.isDisconnected = function(){
	i = NULL
	possibleError <- tryCatch( suppressWarnings(foreach(i = seq_len(2)) %dopar% {i}), error = function(e) e)
	return( inherits(possibleError, "error") && identical(possibleError$message, "invalid connection") )
}


