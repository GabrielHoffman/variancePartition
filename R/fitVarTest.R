
# # Define return objects
# setClass("varPartTest", representation(method="character"), contains="list")
# setClass("varPartCompositeTest", representation(method="character"), contains="list")


# #' Fit linear (mixed) model for hypothesis testing
# #' 
# #' Fit linear (mixed) model before extracting results of hypothesis tests
# #'
# #' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
# #' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)
# #' @param data data.frame with columns corresponding to formula 
# #' @param useWeights if TRUE, analysis uses heteroskadatistic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
# #' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
# #' @param showWarnings show warnings about model fit (default TRUE)
# #' @param colinearityCutoff cutoff used to determine if model is computationally singular
# #' @param ddf method to compute degrees of freedom for hpothesis testing in linear mixed model 
# #' @param ... Additional arguments for lmer() or lm()
# #'
# #' @return
# #' list with summary(fit) evaluated on each model fit
# #' 
# #' @examples
# #' # library(variancePartition)
# #'
# #' # optional step to run analysis in parallel on multicore machines
# #' # Here, we used 4 threads
# #' library(doParallel)
# #' cl <- makeCluster(4)
# #' registerDoParallel(cl)
# #' # or by using the doSNOW package
# #'
# #' # load simulated data:
# #' # geneExpr: matrix of gene expression values
# #' # info: information/metadata about each sample
# #' data(varPartData)
# #' 
# #' # Specify variables to consider
# #' # Age is continuous so we model it as a fixed effect
# #' # Individual and Tissue are both categorical, so we model them as random effects
# #' form <- ~ Age + (1|Individual) + (1|Tissue) 
# #' 
# #' # Fit model for standard hypothesis test
# #' fitTest <- fitVarTest( geneExpr, form, info )
# #' 
# #' # Extract hypothesis test results
# #' res1 = extractTestResults( fitTest, "Age")
# #' 
# #' # stop cluster
# #' stopCluster(cl)
# #'
# ##### @export
# fitVarTest <- function( exprObj, formula, data, useWeights=TRUE, weightsMatrix=NULL, showWarnings=FALSE, colinearityCutoff=.999, ddf="Kenward-Roger",...){ 

# 	f = eval(parse(text=paste("function(fit){
# 		summary(fit, ddf='", ddf, "')
# 	}", sep='')))

# 	res = fitVarPartModel(exprObj, formula, data, REML=TRUE, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings,fxn=f)

#  	res2 = new( "varPartTest", lapply(res, identity), method=res@method)
# 	names(res2) = names(res)

# 	res2
# }

# #' Fit linear (mixed) model for composite hypothesis testing
# #' 
# #' Fit linear (mixed) model before extracting results of composite hypothesis tests
# #'
# #' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
# #' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used a a response. e.g.: ~ a + b + (1|c)
# #' @param data data.frame with columns corresponding to formula 
# #' @param useWeights if TRUE, analysis uses heteroskadatistic error estimates from voom().  Value is ignored unless exprObj is an EList() from voom() or weightsMatrix is specified
# #' @param weightsMatrix matrix the same dimension as exprObj with observation-level weights from voom().  Used only if useWeights is TRUE 
# #' @param showWarnings show warnings about model fit (default TRUE)
# #' @param colinearityCutoff cutoff used to determine if model is computationally singular
# #' @param ddf method to compute degrees of freedom for hpothesis testing in linear mixed model 
# #' @param ... Additional arguments for lmer() or lm()
# #'
# #' @return
# #' list with anova(fit) evaluated on each model fit
# #' 
# #' @examples
# #' # library(variancePartition)
# #'
# #' # optional step to run analysis in parallel on multicore machines
# #' # Here, we used 4 threads
# #' library(doParallel)
# #' cl <- makeCluster(4)
# #' registerDoParallel(cl)
# #' # or by using the doSNOW package
# #'
# #' # load simulated data:
# #' # geneExpr: matrix of gene expression values
# #' # info: information/metadata about each sample
# #' data(varPartData)
# #' 
# #' # Specify variables to consider
# #' # Age is continuous so we model it as a fixed effect
# #' # Individual and Tissue are both categorical, so we model them as random effects
# #' form <- ~ Age + (1|Individual) + (1|Tissue) 
# #' 
# #' # Fit model for composite hypothesis test
# #' fitCompTest <- fitVarCompositeTest( geneExpr, form, info )
# #'
# #' # Extract hypothesis test results from composite test
# #' res2 = extractTestResults( fitCompTest, "Age")
# #' 
# #' # stop cluster
# #' stopCluster(cl)
# #'
# ##### @export
# fitVarCompositeTest <- function( exprObj, formula, data, useWeights=TRUE, weightsMatrix=NULL, showWarnings=FALSE, colinearityCutoff=.999, ddf="Kenward-Roger",...){ 

# 	fit = fitVarPartModel(exprObj, formula, data, showWarnings=showWarnings)

# 	if( fit@method == 'lm' ){

# 		f = function(fit){
# 			anova( fit)
# 		}
# 	}else{		
# 		f = eval(parse(text=paste("function(fit){
# 			anova(fit, ddf='", ddf, "')
# 		}", sep='')))
# 	}

# 	res = fitVarPartModel(exprObj, formula, data, REML=TRUE, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings,fxn=f)

# 	res2 = new( "varPartCompositeTest", lapply(res, identity), method=res@method)
# 	names(res2) = names(res)

# 	res2
# }



# #' Extract hypothesis test results
# #' 
# #' Extract hypothesis test results from model fit with fitVarTest() or fitVarCompositeTest()
# #'
# #' @param x list of models fit with fitVarTest() or fitVarCompositeTest()
# #' @param coef coefficient to be tested
# #' @param ... other arguments 
# #'
# #' @return
# #' data.frame with results for each gene including p-values (P.Value) and p-values corrected for multiple testing using Benjamini-Hochberg method (adj.P.Value).  Also returns coefficient estimate, standard error and t-statistic for standard hypothesis tests.  For composite hypothesis tests, the p-values are based on multiple coefficient and so are not directly interpretable and are not reported for technical reasons. 
# #' 
# #' @examples
# #' # library(variancePartition)
# #'
# #' # optional step to run analysis in parallel on multicore machines
# #' # Here, we used 4 threads
# #' library(doParallel)
# #' cl <- makeCluster(4)
# #' registerDoParallel(cl)
# #' # or by using the doSNOW package
# #'
# #' # load simulated data:
# #' # geneExpr: matrix of gene expression values
# #' # info: information/metadata about each sample
# #' data(varPartData)
# #' 
# #' # Specify variables to consider
# #' # Age is continuous so we model it as a fixed effect
# #' # Individual and Tissue are both categorical, so we model them as random effects
# #' form <- ~ Age + (1|Individual) + (1|Tissue) 
# #' 
# #' # Fit model for standard hypothesis test
# #' fitTest <- fitVarTest( geneExpr, form, info )
# #' 
# #' # Fit model for composite hypothesis test
# #' fitCompTest <- fitVarCompositeTest( geneExpr, form, info )
# #'
# #' # Extract hypothesis test results
# #' res1 = extractTestResults( fitTest, "Age")
# #'
# #' # Extract hypothesis test results from composite test
# #' res2 = extractTestResults( fitCompTest, "Age")
# #' 
# #' # stop cluster
# #' stopCluster(cl)
# #'
# ##### @export
# #' @docType methods
# #' @rdname extractTestResults-method
# setGeneric("extractTestResults", signature="x",
# 	function( x, coef, ... )
#       standardGeneric("extractTestResults")
# )

# ##### @export
# #' @rdname extractTestResults-method
# #' @aliases extractTestResults,varPartTest-method
# setMethod("extractTestResults", "varPartTest",
# 	function( x, coef, ... ){

# 		# check if coef is missing
# 		if( missing(coef) ){
# 			stop("Must specify coef")
# 		}

# 		# check of coef is valid
# 		validCoef = rownames(coef(x[[1]]))

# 		if( ! (coef %in% validCoef) ){
# 			stop("coef was not found in model fit. Valid coefficients to test are:\n", paste(validCoef, collapse=', '), "Note that only fixed effects can be tested, so discrete variables most be modeled as fixed effects for hypothesis testing purposes")
# 		}

# 		# initialize to pass R CMD check
# 		fit = 1

# 		# extract hypothesis test results for each fit
# 		res = foreach(fit = x, .combine=rbind) %dopar% {
# 			values = coef( fit )
# 			values[coef,]
# 		}
# 		rownames(res) = names(x)

# 		# drop df from lmer summary
# 		if( x@method == "lmer" ){
# 			res = res[,colnames(res) != "df"]
# 		}

# 		# change column names and convert to data.frame
# 		colnames(res) = c("estimate", "se",  "t", "P.Value")
# 		res = data.frame(res)

# 		# BH multiple testing correction
# 		res$adj.P.Value = p.adjust( res$P.Value, method="BH")

# 		res
#  	}
# )

# ##### @export
# #' @rdname extractTestResults-method
# #' @aliases extractTestResults,varPartCompositeTest-method
# setMethod("extractTestResults", "varPartCompositeTest",
# 	function( x, coef, ... ){

# 		# check if coef is missing
# 		if( missing(coef) ){
# 			stop("Must specify coef")
# 		}

# 		# check of coef is valid
# 		validCoef = rownames( x[[1]] )

# 		if( ! (coef %in% validCoef) ){			
# 			stop("coef was not found in model fit. Valid coefficients to test are:\n", paste(validCoef, collapse=', '), "Note that only fixed effects can be tested, so discrete variables most be modeled as fixed effects for hypothesis testing purposes")
# 		}

# 		# initialize to pass R CMD check
# 		fit = 1

# 		# extract hypothesis test results for each fit
# 		res = foreach(fit = x, .combine=rbind) %dopar% {
# 			 fit[coef,'Pr(>F)']
# 		}
# 		rownames(res) = names(x)

# 		# change column names and convert to data.frame
# 		colnames(res) = c("P.Value")
# 		res = data.frame(res)

# 		# BH multiple testing correction
# 		res$adj.P.Value = p.adjust( res$P.Value, method="BH")

# 		res
#  	}
# )





	


# # fit = lmer( geneExpr[1,] ~ Age + (1|Individual) + (1|Tissue), info)

# # test = "Age == 0"

# # # apply generalized linear hypothesis test
# # # to get linear combination of coefficients
# # glRes = glht( fit, test)

# # # Get KR degrees of freedom for this test
# # df = get_Lb_ddf( fit, glRes$linfct)
# # df = as.integer(df)

# # # apply generalized linear hypothesis test
# # # with KR degrees of freedon
# # summary(glht( fit, test, df=df, type=adjusted("none")))


# # fitFix = lm( geneExpr[1,] ~ Age + Individual + Tissue, info)
# # head(coef(summary(fitFix)))



# # fitFix = lm( geneExpr[1,] ~ Age , info)
# # head(coef(summary(fitFix)))