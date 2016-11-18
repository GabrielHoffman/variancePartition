

#' plotStratify
#'
#' Plot gene expression stratified by another variable
#'
#' @param formula specify variables shown in the x- and y-axes.  Y-axis should be continuous variable, x-axis should be discrete.
#' @param data data.frame storing continuous and discrete variables specified in formula  
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
#' # Note: This is a newer, more convient interface to plotStratifyBy()
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
#' plotStratify( Expression ~ Tissue, GE )
#'
#' # Omit legend and color boxes grey
#' plotStratify( Expression ~ Tissue, GE, colorBy = NULL)
#'
#' # Specify colors
#' col = c( B="green", A="red", C="yellow")
#' plotStratify( Expression ~ Tissue, GE, colorBy=col, sort=FALSE)
#'
#' @export
plotStratify = function( formula, data, xlab, ylab, main, sortBy, colorBy, sort=TRUE, text=NULL, text.y=1, text.size=5, pts.cex=1, ylim=NULL, legend=TRUE, x.labels=FALSE ){

	mc <- match.call()
	m <- match(c("formula","data"), names(mc), 0L)
	mf <- mc[c(1L, m)]
	mf[[1L]] <- as.name("model.frame")
	mf <- eval(mf, parent.frame())
	data.st <- data.frame(mf)

	if( ncol(data.st) != 2){
		stop("formula must have exactly 2 entries")
	}

	xval = colnames(data.st)[attr(attr(mf, "terms"), "response")+1]
	yval = colnames(data.st)[attr(attr(mf, "terms"), "response")] 

	if( missing(xlab) ){ 
		xlab = xval
	}
	if( missing(sortBy) || is.null(sortBy) ){ 
		sortBy = xval
	}
	if( missing(colorBy) ){
		colorBy = xval
	}	
	if( missing(ylab) ){
		ylab = yval
	}

    # check that sortBy exist in data.st
    if( !(sortBy %in% colnames(data.st)) ){
    	stop(paste("sortBy is not found in colnames(data): sortBy =", sortBy))
    }

    data.st[[yval]] = as.numeric( data.st[[yval]] )

    xpos = 0.5 #text.x * nlevels(data.st[[xval]])
    ypos = text.y * (max(data.st[[yval]]) - min(data.st[[yval]])) + min(data.st[[yval]])

    if( sort ){    	
	    # sort categories by median expression
	    data.st[['reorder']] = reorder(data.st[[sortBy]],data.st[[yval]], FUN=median)
	    ord = "reorder"
    }else{
    	ord = xval
    }
   
    pOut = ggplot( data.st, aes_string(x=ord, y=yval)) + theme_bw() + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ylab(ylab) + xlab(xlab) + theme(plot.title=element_text(hjust=0.5)) 

    if(  is.null(colorBy) || is.na(colorBy) ){
        pOut = pOut + geom_boxplot(color="grey", fill="grey", outlier.colour='black',outlier.shape = 20)
    }else{

    	# if colors are specified and all levels of xval are represented
    	if( sum(levels(data.st[[xval]]) %in% names(colorBy)) == nlevels(data.st[[xval]]) ){

    		i = match(levels(data.st[[ord]]), levels(data.st[[xval]]) )
 			pOut = pOut + geom_boxplot(aes_string(fill=xval), color=colorBy[i], outlier.colour='black',outlier.shape = 20) + scale_fill_manual( values=array(colorBy))
 		}else{

	        # color boxes by colorBy variable in data.st
	        pOut = pOut + geom_boxplot( aes_string(color=colorBy, fill=colorBy), outlier.colour='black', outlier.shape = 20)
	    }

	    # add legend
	    if( legend ){
	    	pOut = pOut + theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key = element_rect(fill="transparent"), axis.text.x=element_text(angle=30), legend.background = element_rect(fill="transparent"))
	    }else{
	    	pOut = pOut + theme(legend.position="none", axis.text.x=element_text(angle=30))
	    }
    }

    if( ! x.labels ){
		pOut = pOut + theme(axis.ticks.x = element_blank(), axis.text.x = element_blank())
	}

    # add median bar
    pOut = pOut + stat_summary(geom = "crossbar", width=0.65, fatten=0, color="black", fun.data = function(x){ return(c(y=median(x), ymin=median(x), ymax=median(x))) })

    if( ! missing(ylim) ){
    	pOut = pOut + ylim(ylim)
    }

    if( ! missing(main) ){
    	 pOut = pOut + ggtitle(main)
   	}

    if( ! missing(text) ){
        pOut = pOut + annotate("text", label = text, x = xpos, y=ypos, size = text.size, hjust=0)
    }

    #pOut = pOut + geom_jitter(size=pts.cex,height=0, width=0, col="black")

    return( pOut )
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
   
    pOut = ggplot( geneExpr, aes_string(x=ord, y=yval)) + theme_bw() + theme( plot.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + ylab(ylab) + xlab(xlab) + theme(plot.title=element_text(hjust=0.5))

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
	    	pOut = pOut + theme(legend.justification=c(1,0), legend.position=c(1,0), legend.key = element_rect(fill="transparent"), axis.text.x=element_text(angle=30), legend.background = element_rect(fill="transparent"))
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
