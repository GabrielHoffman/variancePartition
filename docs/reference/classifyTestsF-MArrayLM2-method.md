# Multiple Testing Genewise Across Contrasts

For each gene, classify a series of related t-statistics as up, down or
not significant.

## Usage

``` r
# S4 method for class 'MArrayLM2'
classifyTestsF(
  object,
  cor.matrix = NULL,
  df = Inf,
  p.value = 0.01,
  fstat.only = FALSE
)
```

## Arguments

- object:

  numeric matrix of t-statistics or an 'MArrayLM2' object from which the
  t-statistics may be extracted.

- cor.matrix:

  covariance matrix of each row of t-statistics. Defaults to the
  identity matrix.

- df:

  numeric vector giving the degrees of freedom for the t-statistics. May
  have length 1 or length equal to the number of rows of tstat.

- p.value:

  numeric value between 0 and 1 giving the desired size of the test

- fstat.only:

  logical, if 'TRUE' then return the overall F-statistic as for 'FStat'
  instead of classifying the test results

## Details

Works like limma::classifyTestsF, except object can have a list of
covariance matrices object\$cov.coefficients.list, instead of just one
in object\$cov.coefficients

## See also

[`limma::classifyTestsF`](https://rdrr.io/pkg/limma/man/classifytestsF.html)
