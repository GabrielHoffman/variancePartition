# Sqrt of co-variance matrix for `dream()` fit

Define generic `vcovSqrt()` for result of
[`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) and
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

## Usage

``` r
vcovSqrt(object, vobj, coef, approx = TRUE)

# S4 method for class 'MArrayLM'
vcovSqrt(object, vobj, coef, approx = TRUE)

# S4 method for class 'MArrayLM2'
vcovSqrt(object, vobj, coef, approx = TRUE)
```

## Arguments

- object:

  `MArrayLM` object return by
  [`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) or
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- vobj:

  `EList` object returned by
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html)

- coef:

  name of coefficient to be extracted

- approx:

  use fast approximation

## Value

Computes factor of covariance matrix so that `vcov(object)` is the same
as `crossprod(vcovSqrt(object))`

## Examples

``` r
# load simulated data:
# geneExpr: matrix of *normalized* gene expression values
# info: information/metadata about each sample
data(varPartData)

form <- ~Batch

fit <- dream(geneExpr[1:2, ], form, info)
fit <- eBayes(fit)

# Compute covariance directly
Sigma <- vcov(fit, geneExpr[1:2, ])

# Compute factor of covariance
S <- crossprod(vcovSqrt(fit, geneExpr[1:2, ]))
```
