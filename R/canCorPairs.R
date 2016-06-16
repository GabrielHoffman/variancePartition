# Gabriel Hoffman
# April 18, 2016

#' canCorPairs
#'
#' Assess correlation between all pairs of variables in a formula 
#'
#' @param formula standard linear model formula (doesn't support random effects currently, so just change the syntax)
#' @param data data.frame with the data for the variables in the formula
#'
#' @details
#' Canonical Correlation Analysis (CCA) is similar to correlation between two vectors, except that CCA can accommodate matricies as well.  For a pair of variables, canCorPairs assesses the degree to which they co-vary and contain the same information.  Variables in the formula can be a continuous variable or a discrete variable expanded to a matrix (which is done in the backend of a regression model).  For a pair of variables, canCorPairs uses CCA to compute the correlation between these variables and returns the pairwise correlation matrix.
#' 
#' Statistically, let rho be the array of correlation values returned by the standard R function cancor to compute CCA.  canCorPairs returns rho / sum(rho) which is the fraction of the maximum possible correlation.   
#' 
#' Note that CCA returns correlations values between 0 and 1
#'
#' @return 
#' Matrix of correlation values between all pairs of variables.
#'
#' @examples
#'
#' # load library
#' # library(variancePartition)
#' 
#' # load simulated data:
#' data(varPartData)
#' 
#' # specify formula
#' form <- ~ Individual + Tissue + Batch + Age + Height
#'
#' # Compute Canonical Correlation Analysis (CCA)
#' # between all pairs of variables
#' # returns absolute corelation value  
#' C = canCorPairs( form, info)
#' 
#' # Plot correlation matrix
#' plotCorrMatrix( C )
#'
#' @export
canCorPairs = function(formula, data){

	if(.isMixedModelFormula(formula, data)) {
		stop("Cannot handle mixed effects models with (1|x) terms.\nUse ~ x instead of ~ (1|x) in formula")
	}

	X = model.matrix(formula, data)

	varLabels = attr(terms(formula),"term.labels")

	# get grouping
	idx = attr(X, "assign")

	variableList = list()

	# extract data
	for( i in unique(idx)){
		if(i==0) next
		variableList[[varLabels[i]]] = X[,idx==i,drop=FALSE]
	}
	
	C = matrix(NA, length(varLabels), length(varLabels))
	colnames(C) = varLabels
	rownames(C) = varLabels
	diag(C) = 1

	pairs = combn(varLabels, 2)

	for(i in 1:ncol(pairs)){
		key1 = pairs[1,i]
		key2 = pairs[2,i]

		fit = cancor( variableList[[key1]], variableList[[key2]] ) 
		C[key1, key2] = sum(fit$cor) / length(fit$cor)
		C[key2, key1] = C[key1, key2]
	}

	return( C )
}
