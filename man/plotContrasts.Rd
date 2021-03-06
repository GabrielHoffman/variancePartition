% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotContrasts.R
\name{plotContrasts}
\alias{plotContrasts}
\title{Plot representation of contrast matrix}
\usage{
plotContrasts(L)
}
\arguments{
\item{L}{contrast matrix}
}
\value{
ggplot2 object
}
\description{
Plot contrast matrix to clarify interpretation of hypothesis tests with linear contrasts
}
\examples{

# load library
# library(variancePartition)

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# get contrast matrix testing if the coefficient for Batch2 is zero 
form <- ~ Batch + (1|Individual) + (1|Tissue) 
L1 = getContrast( geneExpr, form, info, "Batch3")

# get contrast matrix testing if the coefficient for Batch2 is different from Batch3 
form <- ~ Batch + (1|Individual) + (1|Tissue) 
L2 = getContrast( geneExpr, form, info, c("Batch2", "Batch3"))

# combine contrasts into single matrix
L_combined = cbind(L1, L2)

# plot contrasts
plotContrasts( L_combined )

}
