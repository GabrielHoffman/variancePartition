# Multiple Testing Genewise Across Contrasts

For each gene, classify a series of related t-statistics as up, down or
not significant.

## Usage

``` r
classifyTestsF(object, ...)
```

## Arguments

- object:

  numeric matrix of t-statistics or an 'MArrayLM2' object from which the
  t-statistics may be extracted.

- ...:

  additional arguments

## Details

Works like limma::classifyTestsF, except object can have a list of
covariance matrices object\$cov.coefficients.list, instead of just one
in object\$cov.coefficients

## See also

[`limma::classifyTestsF`](https://rdrr.io/pkg/limma/man/classifytestsF.html)
