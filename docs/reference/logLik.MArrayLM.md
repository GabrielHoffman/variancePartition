# Log-likelihood from model fit

Log-likelihood from model fit

## Usage

``` r
# S3 method for class 'MArrayLM'
logLik(object, vobj, ...)
```

## Arguments

- object:

  result of [`lmFit()`](https://rdrr.io/pkg/limma/man/lmFit.html) or
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- vobj:

  `EList` used to fit model

- ...:

  See [`?stats::logLik`](https://rdrr.io/r/stats/logLik.html)
