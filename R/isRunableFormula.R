# Gabriel Hoffman
# April 6, 2021
 
#' Test if formula is full rank on this dataset
#'
#' Test if formula is full rank on this dataset
#' 
#' @param exprObj expression object
#' @param formula formula
#' @param data data
#'
#' @export
isRunableFormula = function( exprObj, formula, data){ 

	isRunable = TRUE

	possibleError <- tryCatch( 
		# variancePartition:::.getContrastInit(res$geneExpr$E, form_mod, data[colnames(res$geneExpr),])
		.getContrastInit(exprObj, formula, data)
		, error = function(e) e)

	mesg <- "the fixed-effects model matrix is column rank deficient"
	if( isTRUE(inherits(possibleError, "error") && grep(mesg, possibleError$message) == 1) ){
		isRunable = FALSE
	}

	isRunable
}