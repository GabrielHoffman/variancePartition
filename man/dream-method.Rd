% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dream.R
\name{dream}
\alias{dream}
\title{Differential expression with linear mixed model}
\usage{
dream(
  exprObj,
  formula,
  data,
  L,
  ddf = c("adaptive", "Satterthwaite", "Kenward-Roger"),
  useWeights = TRUE,
  control = vpcontrol,
  hideErrorsInBackend = FALSE,
  REML = TRUE,
  BPPARAM = SerialParam(),
  ...
)
}
\arguments{
\item{exprObj}{matrix of expression data (g genes x n samples), or \code{ExpressionSet}, or \code{EList} returned by voom() from the limma package}

\item{formula}{specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)}  Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.}

\item{data}{data.frame with columns corresponding to formula}

\item{L}{contrast matrix specifying a linear combination of fixed effects to test}

\item{ddf}{Specifiy "Satterthwaite" or "Kenward-Roger" method to estimate effective degress of freedom for hypothesis testing in the linear mixed model.  Note that Kenward-Roger is more accurate, but is *much* slower.  Satterthwaite is a good enough approximation for most datasets. "adaptive" (Default) uses KR for <= 20 samples.}

\item{useWeights}{if TRUE, analysis uses heteroskedastic error estimates from \code{voom()}.  Value is ignored unless exprObj is an \code{EList()} from \code{voom()} or \code{weightsMatrix} is specified}

\item{control}{control settings for \code{lmer()}}

\item{hideErrorsInBackend}{default FALSE.  If TRUE, hide errors in \code{attr(.,"errors")} and \code{attr(.,"error.initial")}}

\item{REML}{use restricted maximum likelihood to fit linear mixed model. default is TRUE.  See Details.}

\item{BPPARAM}{parameters for parallel evaluation}

\item{...}{Additional arguments for \code{lmer()} or \code{lm()}}
}
\value{
MArrayLM2 object (just like MArrayLM from limma), and the directly estimated p-value (without eBayes)
}
\description{
Fit linear mixed model for differential expression and preform hypothesis test on fixed effects as specified in the contrast matrix \code{L}
}
\details{
A linear (mixed) model is fit for each gene in exprObj, using formula to specify variables in the regression (Hoffman and Roussos, 2021).  If categorical variables are modeled as random effects (as is recommended), then a linear mixed model us used.  For example if formula is \code{~ a + b + (1|c)}, then the model is

\code{fit <- lmer( exprObj[j,] ~ a + b + (1|c), data=data)}

\code{useWeights=TRUE} causes \code{weightsMatrix[j,]} to be included as weights in the regression model.

Note: Fitting the model for 20,000 genes can be computationally intensive.  To accelerate computation, models can be fit in parallel using \code{BiocParallel} to run code in parallel.  Parallel processing must be enabled before calling this function.  See below.

The regression model is fit for each gene separately. Samples with missing values in either gene expression or metadata are omitted by the underlying call to lmer.

Hypothesis tests and degrees of freedom are producted by \code{lmerTest} and \code{pbkrtest} pacakges

While \code{REML=TRUE} is required by \code{lmerTest} when ddf='Kenward-Roger', ddf='Satterthwaite' can be used with \code{REML} as \code{TRUE} or \code{FALSE}.  Since the Kenward-Roger method gave the best power with an accurate control of false positive rate in our simulations, and since the Satterthwaite method with REML=TRUE gives p-values that are slightly closer to the Kenward-Roger p-values, \code{REML=TRUE} is the default.  See Vignette "3) Theory and practice of random effects and REML"
}
\examples{
# library(variancePartition)

# load simulated data:
# geneExpr: matrix of *normalized* gene expression values
# info: information/metadata about each sample
data(varPartData)

form <- ~ Batch + (1 | Individual) + (1 | Tissue)

# Fit linear mixed model for each gene
# run on just 10 genes for time
# NOTE: dream() runs on *normalized* data
fit <- dream(geneExpr[1:10, ], form, info)
fit <- eBayes(fit)

# view top genes
topTable(fit, coef = "Batch2", number = 3)

# get contrast matrix testing if the coefficient for Batch3 is
# different from coefficient for Batch2
# Name this comparison as 'compare_3_2'
# The variable of interest must be a fixed effect
L <- makeContrastsDream(form, info, contrasts = c(compare_3_2 = "Batch3 - Batch2"))

# plot contrasts
plotContrasts(L)

# Fit linear mixed model for each gene
# run on just 10 genes for time
fit2 <- dream(geneExpr[1:10, ], form, info, L)
fit2 <- eBayes(fit2)

# view top genes for this contrast
topTable(fit2, coef = "compare_3_2", number = 3)

# Parallel processing using multiple cores with reduced memory usage
param <- SnowParam(4, "SOCK", progressbar = TRUE)
fit3 <- dream(geneExpr[1:10, ], form, info, L, BPPARAM = param)
fit3 <- eBayes(fit3)

# Fit fixed effect model for each gene
# Use lmFit in the backend
form <- ~Batch
fit4 <- dream(geneExpr[1:10, ], form, info, L)
fit4 <- eBayes(fit4)

# view top genes
topTable(fit4, coef = "compare_3_2", number = 3)

# Compute residuals using dream
residuals(fit4)[1:4, 1:4]

}
\references{
{
  \insertRef{hoffman2021dream}{variancePartition}
}
}
