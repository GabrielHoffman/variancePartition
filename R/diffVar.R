# Gabriel Hoffman
# Nov 17, 2022
# 
# A test of differential variance

#' Test differential variance
#' 
#' Test the association between a covariate of interest and the response's deviation from expectation. 
#' 
#' @param fit model fit from \code{dream()}
#' @param method transform the residuals using absolute deviation ("AD") or squared deviation "SQ".
#' @param BPPARAM parameters for parallel evaluation
#' @param ... other parameters passed to \code{dream()}
#'  
#' @details
#' This method performs a test of differential variance between two subsets of the data, in a way that generalizes to multiple categories, continuous variables and metrics of spread beyond variance.  For the two category test, this method is simular to Levene's test.  This model was adapted from Phipson, et al (2014), extended to linear mixed models, and adapted to be compatible with \code{dream()}.
#' 
#' This method is composed of multiple steps where 1) a typical linear (mixed) model is fit with \code{dream()}, 2) residuals are computed and transformed based on an absolute value or squaring transform, 3) a second regression is performed with \code{dream()} to test if a variable is associated with increased deviation from expectation.  Both regression take advantage of the \code{dream()} linear (mixed) modelling framework followed by empirical Bayes shrinkage that extends the \code{limma::voom()} framework.
#'
#' Note that \code{diffVar()} takes the results of the first regression as a parameter to use as a starting point.
#' 
#' @references{
#'   \insertRef{phipson2014diffvar}{variancePartition}
#' }
#' 
#' @examples
#' # library(variancePartition)
#' library(edgeR)
#' data(varPartDEdata)
#' 
#' # filter genes by number of counts
#' isexpr = rowSums(cpm(countMatrix)>0.1) >= 5
#' 
#' # Standard usage of limma/voom
#' geneExpr = DGEList( countMatrix[isexpr,] )
#' geneExpr = calcNormFactors( geneExpr )
#' 
#' # make this vignette faster by analyzing a subset of genes
#' geneExpr = geneExpr[1:1000,]
#' 
#' # regression formula
#' form <- ~ Disease 
#' 
#' # estimate precision weights
#' vobj = voomWithDreamWeights( geneExpr, form, metadata )
#' 
#' # fit dream model
#' fit = dream( vobj, form, metadata )
#' fit = eBayes(fit)
#' 
#' # fit differential variance model
#' res = diffVar( fit )
#' 
#' # extract results for differential variance based on Disease
#' topTable(res, coef = "Disease1", number=3)
#' 
#' # Box plot of top hit
#' # Since ASCL3 has a negative logFC, 
#' # the deviation from expectation is *smaller* in 
#' # Disease==1 compared to baseline.
#' gene = "ENST00000325884.1 gene=ASCL3"
#' boxplot(vobj$E[gene,] ~ metadata$Disease, main=gene)
#' 
#' @seealso \code{missMethyl::diffVar()}, \code{car::leveneTest()}
#' @export
#' 
#' @docType methods
#' @rdname diffVar-method
setGeneric("diffVar", signature="fit",
	function( fit, method = c("AD", "SQ"), BPPARAM = SerialParam(), ... )
      standardGeneric("diffVar")
)


#' @export
#' @rdname diffVar-method
#' @aliases diffVar,MArrayLM2-method
setMethod("diffVar", "MArrayLM",
	function( fit, method = c("AD", "SQ"),
		BPPARAM = SerialParam(), ...){

	method = match.arg(method)

	# 1) extract residuals
	z = residuals(fit)

	# extract leverage
	hv = fit$hatvalues

	# fixed effects models have the same hatvalues for each response
	# convert array into matrix with the same rows
	if( ! is.matrix(hv) ){
		hv = t(matrix(hv, ncol(z), nrow(z)))
	}

	# scale residuals by sqrt(1-h)
	z = z / sqrt(1-hv)

	# 2) transform residuals
	if( method == "AD"){
		z = abs(z)
	}else if( method == "SQ"){
		z = z^2
	}

	# 3) fit regression on deviations
	fit2 = dream(z, fit$formula, fit$data, BPPARAM=BPPARAM,..., quiet=TRUE)
	fit2 = eBayes(fit2)

	fit2
})


