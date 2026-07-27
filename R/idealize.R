

#' Idealized Variance Fractions
#' 
#' Compute idealized variance fractions removing \code{CountNoise} or other variables
#' 
#' @param x \code{matrix} or \code{data.frame} of variance fracdtions
#' @param remove colnames of \code{x} to remove from the variance fractions
#' 
#' @details Variance partitioning analysis includes all components of variance in the denominator when computing the variance fractions.  This function recomputes the variance fractions in idealized case where \code{CountNoise} or any set of variances given in the \code{remove} argument is zero.
#'
#' @examples
#' # Simulate counts
#' set.seed(1)
#' countMatrix <- matrix(rnbinom(n=100000, mu=20, size=3), ncol=10)
#' rownames(countMatrix) <- paste0("gene_", seq(nrow(countMatrix)))
#' colnames(countMatrix) <- paste0("sample_", seq(ncol(countMatrix)))
#' 
#' condition <- factor(rep(1:2, each=5))
#' 
#' # DESeq2 model #
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(countMatrix, 
#'   DataFrame(condition), 
#'   ~ condition)
#' dds <- DESeq(dds)
#' res <- results(dds)
#' 
#' # Variance partition analysis
#' vp1 <- varpart(dds)
#'
#' # Compute variance fractions in the idealized
#' #  case where there is no count noise
#' head(idealized(vp1))
#
#' @export
idealized <- function(x, remove = "CountNoise"){

  exclude <- match(remove, colnames(x))

  if( any(is.na(exclude)) ){
    stop("Entry in remove not found in colnames of x: ", paste(remove[is.na(exclude)], collapse=', '))
  }

  x.exclude <- x[,-exclude,drop=FALSE]

  x.exclude / rowSums(x.exclude)
}