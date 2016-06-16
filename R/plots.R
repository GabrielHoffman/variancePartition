

#' Violin plot of variance fractions
#' 
#' Violin plot of variance fraction for each gene and each variable
#'
#' @param obj varParFrac object returned by fitExtractVarPart or extractVarPart
#' @param col vector of colors
#' @param label.angle angle of labels on x-axis
#' @param main title of plot
#' @param convertToPercent multiply fractions by 100 to convert to percent values
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
	function( obj, col=c(ggColorHue(ncol(obj)-1), "#bebebe99"), label.angle=20, main="", convertToPercent=TRUE,...)
      standardGeneric("plotVarPart")
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,matrix-method
setMethod("plotVarPart", "matrix",
	function( obj, col=c(ggColorHue(ncol(obj)-1), "#bebebe99"), label.angle=20, main="", convertToPercent=TRUE, ...){
 		.plotVarPart( obj, col, label.angle, main, convertToPercent,...)
 	}
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,varPartResults-method
setMethod("plotVarPart", "data.frame",
	function( obj, col=c(ggColorHue(ncol(obj)-1), "#bebebe99"), label.angle=20, main="", convertToPercent=TRUE,...){
 		.plotVarPart( obj, col, label.angle, main, convertToPercent,... )
 	}
)

#' @export
#' @rdname plotVarPart-method
#' @aliases plotVarPart,matrix-method
setMethod("plotVarPart", "varPartResults",
	function( obj, col=c(ggColorHue(ncol(obj)-1), "#bebebe99"), label.angle=20, main="", convertToPercent=TRUE, ...){
		
		# don't convert if values are actual variances
		convertToPercent = !(obj@method == "Variance (log2 CPM scale)")

 		.plotVarPart( data.frame(obj, check.names=FALSE), col, label.angle, main, ylab=obj@method, convertToPercent,...)
 	}
)

# internal driver function
.plotVarPart <- function( obj, col=c(ggColorHue(ncol(obj)-1), "#bebebe99"), label.angle=20, main="", ylab='', convertToPercent=TRUE, ylim,...){

	# convert to data.frame
	obj = as.data.frame(obj, check.names=FALSE)

	if( length(col) < ncol(obj) ){
		stop("Not enough colors specified by col")
	}

	# get gene name of each row
	obj$gene <- rownames(obj)

	# convert to data.frame for ggplot
	data <- melt(obj, id="gene")

	if( convertToPercent ){
		data$value <- data$value * 100

		if( missing(ylim) ){
			ylim = c(0, 100)
		}

	}else{
		if( missing(ylim)){
			ylim = c(0, max(data$value))
		}
	}

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
plotPercentBars = function( varPart, col = c(ggColorHue(ncol(varPart)-1), "#bebebe99") ){

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


