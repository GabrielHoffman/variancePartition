
#' Plot representation of contrast matrix
#' 
#' Plot contrast matrix to clarify interpretation of hypothesis tests with linear contrasts
#'
#' @param L contrast matrix
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
#' # geneExpr: matrix of gene expression values
#' # info: information/metadata about each sample
#' data(varPartData)
#' 
#' # get contrast matrix testing if the coefficient for Batch2 is zero 
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L1 = getContrast( geneExpr, form, info, "Batch3")
#' 
#' # get contrast matrix testing if the coefficient for Batch2 is different from Batch3 
#' form <- ~ Batch + (1|Individual) + (1|Tissue) 
#' L2 = getContrast( geneExpr, form, info, c("Batch2", "Batch3"))
#' 
#' # combine contrasts into single matrix
#' L_combined = cbind(L1, L2)
#' 
#' # plot contrasts
#' plotContrasts( L_combined )
#' 
#' @importFrom reshape2 melt
#' @import ggplot2
#' @export
plotContrasts = function( L ){
	df = melt(t(L))
	df$Var1 = factor(df$Var1)

	if( identical(levels(df$Var1), "1") ){
		df$Var1 = factor(rep('', nrow(df)))
	}

	Var1 = Var2 = value = NULL

	h = length(unique(df$Var1))
	w = length(unique(df$Var2))

	ggplot(df, aes(Var2, y=Var1, fill=value)) + geom_tile(color="black") + theme_minimal() + theme(aspect.ratio=h/w,
		panel.grid.major =element_blank(),
		panel.grid.minor =element_blank(),plot.title = element_text(hjust = 0.5)) + scale_fill_gradient2(name="Contrast coef", limits=c(-1,1), low=alpha("blue",.8), mid="white", high=alpha("red", .8)) + xlab("Variable") + ylab("Contrasts") + ggtitle("Graphical representation of linear contrasts") + geom_text(aes(label=value), fontface = "bold")
}

