# Gabriel Hoffman
# March 4, 2015

#' plotCorrStructure
#'
#' Plot correlation structure of a gene based on random effects
#'
#' @param fit linear mixed model fit of a gene produced by lmer() or fitVarPartModel()
#' @param varNames variables in the metadata for which the correlation structure should be shown.  Variables must be random effects
#' @param reorder how to reorder the rows/columns of the correlation matrix.  reorder=FALSE gives no reorder.  reorder=TRUE reorders based on hclust.  reorder can also be an array of indeces to reorder the samples manually
#' @param pal color palette  
#' @param hclust.method clustering methods for hclust
#'
#' @return
#' Image of correlation structure beteen each pair of experiments for a single gene
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
#' data(varPartData)
#' 
#' # specify formula
#' form <- ~ Age + (1|Individual) + (1|Tissue)
#' 
#' # fit and return linear mixed models for each gene
#' fitList <- fitVarPartModel( geneExpr[1:10,], form, info )
#' 
#' # Focus on the first gene
#' fit = fitList[[1]]
#' 
#' # plot correlation sturcture based on Individual, reordering samples with hclust
#' plotCorrStructure( fit, "Individual" )
#' 
#' # don't reorder
#' plotCorrStructure( fit, "Individual", reorder=FALSE )
#' 
#' # plot correlation sturcture based on Tissue, reordering samples with hclust
#' plotCorrStructure( fit, "Tissue" )
#' 
#' # don't reorder
#' plotCorrStructure( fit, "Tissue", FALSE )
#' 
#' # plot correlation structure based on all random effects
#' # reorder manually by Tissue and Individual
#' idx = order(info$Tissue, info$Individual)
#' plotCorrStructure( fit, reorder=idx )
#' 
#' # plot correlation structure based on all random effects
#' # reorder manually by Individual, then Tissue 
#' idx = order(info$Individual, info$Tissue)
#' plotCorrStructure( fit, reorder=idx )
#' 
#' # stop cluster
#' stopCluster(cl)
#' 
#' @export
plotCorrStructure = function( fit, varNames = names(coef(fit)), reorder=TRUE, 
	pal = colorRampPalette(c("white", "red", "darkred")), hclust.method = "complete" ){

	if( ! inherits(fit, "lmerMod") ){
		stop("fit must be a linear mixed model fit of class lmerMod")
	}

	# If any variable names are not in fit
	if( any(!varNames %in% names(coef(fit))) ){

		# indeces of variables not found in model fit
		idx = which(!varNames %in% names(coef(fit)))

		stop(paste("Variables are not found in model fit:", paste(varNames[idx], collapse=", "))) 
	}

	# get correlation terms
	vp = calcVarPart( fit )

	# get structure of study design
	sigG = pbkrtest::get_SigmaG( fit )

	# get variable names
	ids = names(coef(fit))

	# store correlation matrix for each variable
	C = list()
	SideColors = c()
	Sigma = 0

	# get the metadata from the model fit
	data = fit@frame

	# for each desired variable compute the correlation matrix
	for( key in varNames){
		# get the variable index
		i = which( key == ids)

		# multiply the sparse matrix by the correlation value
		C[[i]] = sigG$G[[i]] * vp[[key]]

		# sum up with previous correlation matrix
		Sigma = Sigma + C[[i]]

		# convert the relevant variable to a factor
		variable = factor(data[,ids[i]])

		# color the variable levels and cbind them
		sideColset = ggColorHue(nlevels(variable))
		SideColors = cbind(SideColors, sideColset[variable])
	}
	colnames(SideColors) = varNames
	diag(Sigma) = 1 # set diagonals to 1

	main = paste(varNames, collapse=" + ")


	if( is.logical(reorder) && reorder ){
		# reorder the samples by hclust
		hcr <- hclust(as.dist(1-as.matrix(Sigma)), method=hclust.method)
		ddr <- as.dendrogram(hcr)
		ddr <- stats::reorder(ddr, 1)
		idx <- order.dendrogram(ddr)

	}else if( length(reorder) == nrow(Sigma) ){
		# reorder based on pre-specified order
		idx = reorder
	}else{
		# don't reorder
		idx = 1:ncol(Sigma)
	}

	heatmap.3( as.matrix(Sigma)[idx,idx], dendrogram="none", col=pal(100), trace="none", tracecol=FALSE, symm=TRUE, ColSideColors=SideColors[idx,,drop=FALSE], RowSideColors=t(SideColors[idx,ncol(SideColors):1]), keysize = 1, key.title="", density.info="none", KeyValueName = "correlation", main=main, labRow='', labCol='', margins=c(2,2), Rowv=FALSE, Colv=FALSE, NumColSideColors=ncol(SideColors)*.5, NumRowSideColors=ncol(SideColors)*.5 )
}





# library(variancePartition)


# source("~/workspace/Bioconductor/variancePartition/R/heatmap.R")
