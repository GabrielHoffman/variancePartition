# Co-variance matrix for `dream()` fit

Define generic [`vcov()`](https://rdrr.io/r/stats/vcov.html) for result
of [`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) and
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

## Usage

``` r
# S4 method for class 'MArrayLM2'
vcov(object, vobj, coef)
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

## Value

variance-covariance matrix
