# eBayes generic for for MArrayLM and MArrayLM2

eBayes for result of linear (mixed) model

eBayes for result of linear model for with uses
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)

eBayes for result of linear mixed model for with
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
using residual degrees of freedom approximated with
[`rdf.merMod()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/rdf.merMod.md)

## Usage

``` r
eBayes(fit, ...)

# S4 method for class 'MArrayLM'
eBayes(fit, ...)

# S4 method for class 'MArrayLM2'
eBayes(
  fit,
  proportion = 0.01,
  stdev.coef.lim = c(0.1, 4),
  trend = FALSE,
  span = NULL,
  robust = FALSE,
  winsor.tail.p = c(0.05, 0.1),
  legacy = NULL
)
```

## Arguments

- fit:

  fit

- ...:

  all other args

- proportion:

  proportion

- stdev.coef.lim:

  stdev.coef.lim

- trend:

  trend

- span:

  span

- robust:

  robust

- winsor.tail.p:

  winsor.tail.p

- legacy:

  legacy

## Value

results of eBayes using approximated residual degrees of freedom

## See also

[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md),
[`rdf.merMod()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/rdf.merMod.md),
[`limma::eBayes()`](https://rdrr.io/pkg/limma/man/ebayes.html)
