# Linear mixed model confidence intervals

Fit linear mixed model to estimate contribution of multiple sources of
variation while simultaneously correcting for all other variables. Then
perform parametric bootstrap sampling to get a 95% confidence intervals
for each variable for each gene.

## Usage

``` r
varPartConfInf(
  exprObj,
  formula,
  data,
  REML = FALSE,
  useWeights = TRUE,
  control = vpcontrol,
  nsim = 1000,
  ...
)
```

## Arguments

- exprObj:

  matrix of expression data (g genes x n samples), or `ExpressionSet`,
  or `EList` returned by
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html) from the `limma`
  package

- formula:

  specifies variables for the linear (mixed) model. Must only specify
  covariates, since the rows of exprObj are automatically used as a
  response. e.g.: `~ a + b + (1|c)`

- data:

  `data.frame` with columns corresponding to formula

- REML:

  use restricted maximum likelihood to fit linear mixed model. default
  is FALSE. Strongly discourage against changing this option, but here
  for compatibility.

- useWeights:

  if TRUE, analysis uses heteroskedastic error estimates from
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html). Value is ignored
  unless exprObj is an `EList` from
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html) or `weightsMatrix`
  is specified

- control:

  control settings for
  [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

- nsim:

  number of bootstrap datasets

- ...:

  Additional arguments for
  [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html) or l`m()`

## Value

[`list()`](https://rdrr.io/r/base/list.html) of where each entry is the
result for a gene. Each entry is a matrix of the 95% confidence interval
of the variance fraction for each variable

## Details

A linear mixed model is fit for each gene, and
[`bootMer()`](https://rdrr.io/pkg/lme4/man/bootMer.html) is used to
generate parametric bootstrap confidence intervals. `use.u=TRUE` is used
so that the \\\hat{u}\\ values from the random effects are used as
estimated and are not re-sampled. This gives confidence intervals as if
additional data were generated from these same current samples.
Conversely, `use.u=FALSE` assumes that this dataset is a sample from a
larger population. Thus it simulates \\\hat{u}\\ based on the estimated
variance parameter. This approach gives confidence intervals as if
additional data were collected from the larger population from which
this dataset is sampled. Overall, `use.u=TRUE` gives smaller confidence
intervals that are appropriate in this case.

## Examples

``` r

# load library
# library(variancePartition)

library(BiocParallel)

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# Specify variables to consider
# Age is continuous so we model it as a fixed effect
# Individual and Tissue are both categorical, so we model them as random effects
form <- ~ Age + (1 | Individual) + (1 | Tissue)

# Compute bootstrap confidence intervals for each variable for each gene
resCI <- varPartConfInf(geneExpr[1:5, ], form, info, nsim = 100)
```
