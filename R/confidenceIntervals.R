


setClass("VarParCIList", representation(method="character"), contains="list")


#' Linear mixed model confidence intervals
#' 
#' Fit linear mixed model to estimate contribution of multiple sources of variation while simultaneously correcting for all other variables. Then perform parametric bootstrap sampling to get a 95\% confidence intervals for each variable for each gene.
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
#' @param colinearityCutoff cutoff used to determine if model is computationally singular
#' @param control control settings for lmer()
#' @param nsim number of bootstrap datasets
#' @param ... Additional arguments for lmer() or lm()
#' 
#' @return
#' list() of where each entry is the result for a gene.  Each entry is a matrix of the 95% confidence interval of the variance fraction for each variable 
#'
#' @details 
#' A linear mixed model is fit for each gene, and bootMer() is used to generate parametric boostrap confidence intervals.  use.u=TRUE is used so that the \\hat(u) values from the random effects are used as estimated and are not resampled.  This gives confindence intervals as if additional data were generated from these same current samples.  Conversely, use.u=FALSE assumes that  this dataset is a sample from a larger population.   Thus it simulates \\hat(u) based on the estimated variance parameter.  This approach gives confidence intervals as if additional data were collected from the larger population from which this dataset is sampled.  Overall, use.u=TRUE gives smaller confidence intervals that are appropriate in this case.
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
#' # Compute bootstrap confidence intervals for each variable for each gene
#' resCI <- varPartConfInf( geneExpr[1:5,], form, info, nsim=100 )
#' 
#' # stop cluster
#' stopCluster(cl)
#'
#' @export
varPartConfInf <- function( exprObj, formula, data, REML=FALSE, useWeights=TRUE, weightsMatrix=NULL, adjust=NULL, adjustAll=FALSE, showWarnings=TRUE, colinearityCutoff=.999, control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" ),nsim=1000,...){ 

	if( !is.numeric(nsim) || nsim <= 0 ){
		stop("Must specify nsim as positive number")
	}

	# define bootstrap function
	bootStrapFxn = local(function( fit ){

		if( "lmerMod" %in% class(fit)){
			# use bootMer from lme4 to sample to bootstraps and refit the model
			# use.u=TRUE says that the \hat{u} values from the random effects are used
			bootRes = lme4::bootMer( fit, calcVarPart, use.u=TRUE, nsim=nsim, parallel="no", ncpus=1)
			apply( bootRes$t, 2, function(x) quantile(x, c(0.025, 0.975)))
		}else if( "lm" %in% class(fit)){
			stop("Bootstrap of fixed effect ANOVA model not currently supported")
		}
	})

	# fit the model and run the bootstrap function for each gene
	res = fitVarPartModel( exprObj=exprObj, formula=formula, data=data, REML=REML, useWeights=useWeights, weightsMatrix=weightsMatrix, showWarnings=showWarnings,fxn=bootStrapFxn, colinearityCutoff=colinearityCutoff, control = control,...) #  adjust=adjust, adjustAll=adjustAll, 

	res@method="bootstrap"

	return(res)
}




# form = ~ Age + Individual + Tissue

#  fit = lm(geneExpr[1,]~ Age + Individual + Tissue, info, weights=1:nrow(info))



# bs2 <- function(formula, data, indices) {
# 	d <- data[indices,]
# 	fit2 <- lm(formula, data=d, weights=fit$weights[indices])
# 	return( calcVarPart(fit2) ) 
# } 

# bootRes <- boot(data=data.frame(info, gene=geneExpr[1,]), statistic=bs2, R=1000, formula= paste("gene ~", as.character(form))



# apply( bootRes$t, 2, function(x) quantile(x, c(0.025, 0.975)))

# calcVarPart( fit )




