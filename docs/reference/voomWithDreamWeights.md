# Transform RNA-Seq Data Ready for Linear Mixed Modelling with `dream()`

Transform count data to log2-counts per million (logCPM), estimate the
mean-variance relationship and use this to compute appropriate
observation-level weights. The data are then ready for linear mixed
modelling with
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md).
This method is the same as
[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html), except that
it allows random effects in the formula

## Usage

``` r
voomWithDreamWeights(
  counts,
  formula,
  data,
  lib.size = NULL,
  normalize.method = "none",
  span = 0.5,
  weights = NULL,
  prior.count = 0.5,
  prior.count.for.weights = prior.count,
  plot = FALSE,
  save.plot = TRUE,
  rescaleWeightsAfter = FALSE,
  scaledByLib = FALSE,
  priorWeightsAsCounts = FALSE,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- counts:

  a numeric `matrix` containing raw counts, or an `ExpressionSet`
  containing raw counts, or a `DGEList` object. Counts must be
  non-negative and NAs are not permitted.

- formula:

  specifies variables for the linear (mixed) model. Must only specify
  covariates, since the rows of exprObj are automatically used as a
  response. e.g.: `~ a + b + (1|c)` Formulas with only fixed effects
  also work, and [`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html)
  followed by contrasts.fit() are run.

- data:

  `data.frame` with columns corresponding to formula

- lib.size:

  numeric vector containing total library sizes for each sample.
  Defaults to the normalized (effective) library sizes in `counts` if
  `counts` is a `DGEList` or to the columnwise count totals if `counts`
  is a matrix.

- normalize.method:

  the microarray-style normalization method to be applied to the logCPM
  values (if any). Choices are as for the `method` argument of
  `normalizeBetweenArrays` when the data is single-channel. Any
  normalization factors found in `counts` will still be used even if
  `normalize.method="none"`.

- span:

  width of the lowess smoothing window as a proportion. Setting
  `span="auto"` uses
  [`fANCOVA::loess.as()`](https://rdrr.io/pkg/fANCOVA/man/loess.as.html)
  to estimate the tuning parameter from the data

- weights:

  Can be a numeric matrix of individual weights of same dimensions as
  the `counts`, or a numeric vector of sample weights with length equal
  to `ncol(counts)`

- prior.count:

  average count to be added to each observation to avoid taking log of
  zero. The count applied to each sample is normalized by library size
  so given equal log CPM for a gene with zero counts across multiple
  samples

- prior.count.for.weights:

  count added to regularize weights

- plot:

  logical, should a plot of the mean-variance trend be displayed?

- save.plot:

  logical, should the coordinates and line of the plot be saved in the
  output?

- rescaleWeightsAfter:

  default = FALSE, should the output weights be scaled by the input
  weights

- scaledByLib:

  if `TRUE`, scale pseudocount by `lib.size`. Else to standard constant
  pseudocount addition

- priorWeightsAsCounts:

  if `weights` is `NULL`, set weights to be equal to counts, following
  delta method for log2 CPM

- BPPARAM:

  parameters for parallel evaluation

- ...:

  other arguments are passed to `lmer`.

## Value

An `EList` object just like the result of
[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html)

## Details

Adapted from [`voom()`](https://rdrr.io/pkg/limma/man/voom.html) in
`limma` v3.40.2

## See also

[`limma::voom()`](https://rdrr.io/pkg/limma/man/voom.html)

## Examples

``` r
# library(variancePartition)
library(edgeR)
library(BiocParallel)

data(varPartDEdata)

# normalize RNA-seq counts
dge <- DGEList(counts = countMatrix)
dge <- calcNormFactors(dge)

# specify formula with random effect for Individual
form <- ~ Disease + (1 | Individual)

# compute observation weights
vobj <- voomWithDreamWeights(dge[1:20, ], form, metadata)

# fit dream model
res <- dream(vobj, form, metadata)
res <- eBayes(res)

# extract results
topTable(res, coef = "Disease1", number = 3)
#>                                    logFC  AveExpr        t      P.Value
#> ENST00000456159.1 gene=MET     1.0182945 2.458926 6.241638 7.270470e-07
#> ENST00000418210.2 gene=TMEM64  1.0375652 4.715367 6.424903 3.175466e-06
#> ENST00000555834.1 gene=RPS6KL1 0.9355651 5.272063 5.653604 3.749850e-06
#>                                   adj.P.Val        B    z.std
#> ENST00000456159.1 gene=MET     1.454094e-05 5.834316 4.953996
#> ENST00000418210.2 gene=TMEM64  2.499900e-05 5.811378 4.659131
#> ENST00000555834.1 gene=RPS6KL1 2.499900e-05 4.223847 4.624786
```
