
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


