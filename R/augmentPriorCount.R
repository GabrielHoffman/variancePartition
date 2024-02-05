#' Augment observed read counts with prior count
#' 
#' Augment observed read counts with prior count since log of zero counts is undefined.  The prior count added to each sample is scaled so that no variance is introduced
#' 
#' @param counts matrix of read counts with genes as rows and samples as columns
#' @param lib.size library sizes, the sum of all ready for each sample 
#' @param prior.count average prior count added to each sample. 
#' @param scaledByLib if \code{TRUE}, scale pseudocount by \code{lib.size}.  Else to standard constant pseudocount addition 
#' 
#' @return
#' matrix with augmented counts
#'
#' @seealso \code{edgeR::cpm()}
#' @details Adding prior counts removes the issue of evaluating the log of zero counts, and stabilizes the log transform when counts are very small. However, adding a constant prior count to all samples can introduced an artifact. Consider two samples each with zero counts for a given gene, but one as a library size of 1k and the other of 50k. After applying the prior count values become pc / 1k and pc / 50k. It appears that there is variance in the expression of this gene, even though no counts are observed. This is driven only by variation in the library size, which does not reflect biology. This issue is most problematic for small counts.
#'
#' Instead, we make the reasonable assumption that a gene does not have expression variance unless supported sufficiently by counts in the numerator.  Consider adding a different prior count to each sample so that genes with zero counts end up woth zero variance. This corresponds to adding \code{prior.count * lib.size[i] / mean(lib.size)} to sample \code{i}.  
#'
# Not that \code{prior.count} should be small compared to the library size.  Adding a prior count of 2 in the case of bulk RNA-seq with a library size of 30M reads is reasonable.  But for single cell data with 5k reads per cell, adding 2 counts to each gene may be too much.
#'
#' This is done in the backend of \code{edgeR::cpm()}, but this function allows users to apply it more generally.
#'
#' @examples
# library(variancePartition)
#' library(edgeR)
#' 
#' data(varPartDEdata)
#' 
#' # normalize RNA-seq counts
#' dge <- DGEList(counts = countMatrix)
#' dge <- calcNormFactors(dge)
#' 
#' countsAugmented <- augmentPriorCount( dge$counts, dge$samples$lib.size, 1)
# 
#' @importFrom matrixStats colSums2 
#' @export
augmentPriorCount = function(counts, 
					lib.size = colSums2(counts), 
					prior.count = 0.5, scaledByLib = FALSE){

	if( ! scaledByLib ){
		# do standard, constant pseudocount addition
		lib.size = rep(1, ncol(counts))
	}
	stopifnot(length(lib.size) == ncol(counts))

	# prior count scaled by library size deviation
	# constant for all genes within a sample
	lambda = prior.count * lib.size / mean(lib.size)

	# matrix of lambda values
	pcMat = matrix(lambda, nrow=nrow(counts), ncol=ncol(counts), byrow=TRUE)

	# augment count matrix with prior counts
	counts + pcMat
}