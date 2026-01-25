# Multivariate tests on results from `dream()`

Evaluate multivariate tests on results from
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
using [`vcov()`](https://rdrr.io/r/stats/vcov.html) to compute the
covariance between estimated regression coefficients across multiple
responses. A joint test to see if the coefficients are jointly different
from zero is performed using meta-analysis methods that account for the
covariance.

## Usage

``` r
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'MArrayLM,EList,vector'
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'MArrayLM,EList,missing'
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'MArrayLM,EList,list'
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'mvTest_input,ANY,ANY'
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)

# S4 method for class 'MArrayLM,matrix,ANY'
mvTest(
  fit,
  vobj,
  features,
  coef,
  method = c("FE.empirical", "FE", "RE2C", "tstat", "hotelling", "sidak", "fisher"),
  shrink.cov = TRUE,
  BPPARAM = SerialParam(),
  ...
)
```

## Arguments

- fit:

  `MArrayLM` or `MArrayLM2` returned by
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- vobj:

  matrix or `EList` object returned by
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html)

- features:

  a\) indeces or names of features to perform multivariate test on, b)
  list of indeces or names. If missing, perform joint test on all
  features.

- coef:

  name of coefficient or contrast to be tested

- method:

  statistical method used to perform multivariate test. See details.
  `'FE'` is a fixed effect test that models the covariance between
  coefficients. `'FE.empirical'` use compute empirical p-values by
  sampling from the null distribution and fitting with a gamma. `'RE2C'`
  is a random effect test of heterogeneity of the estimated coefficients
  that models the covariance between coefficients, and also incorporates
  a fixed effects test too. `'tstat'` combines the t-statistics and
  models the covariance between coefficients. `'hotelling'` performs the
  Hotelling T2 test. `'sidak'` returns the smallest p-value and
  accounting for the number of tests. `'fisher'` combines the p-value
  using Fisher's method assuming independent tests.

- shrink.cov:

  shrink the covariance matrix between coefficients using the
  Schafer-Strimmer method

- BPPARAM:

  parameters for parallel evaluation

- ...:

  other arugments

## Value

Returns a `data.frame` with the statistics from each test, the `pvalue`
from the test, `n_features`, `method`, and `lambda` from the
Schafer-Strimmer method to shrink the estimated covariance. When
`shrink.cov=FALSE`, `lambda = 0`.

## Details

See package `remaCor` for details about the
[`remaCor::RE2C()`](http://gabrielhoffman.github.io/remaCor/reference/RE2C.md)
test, and see
[`remaCor::LS()`](http://gabrielhoffman.github.io/remaCor/reference/LS.md)
for details about the fixed effect test. When only 1 feature is
selected, the original p-value is returned and the test statistic is set
to `NA`.

For the `"RE2C"` test, the final test statistic is the sum of a test
statistic for the mean effect (`stat.FE`) and heterogeneity across
effects (`stat.het`). `mvTest()` returns 0 if `stat.het` is negative in
extremely rare cases.

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
fit <- dream(vobj, form, metadata)
fit <- eBayes(fit)

# Multivariate test of features 1 and 2
mvTest(fit, vobj, 1:2, coef = "Disease1")
#>        beta        se     stat       pvalue n_features lambda       method
#> 1 0.9301722 0.1325913 7.015336 7.022133e-10          2   0.01 FE.empirical

# Test multiple sets of features
lst <- list(a = 1:2, b = 3:4)
mvTest(fit, vobj, lst, coef = "Disease1", BPPARAM = SnowParam(2))
#>   ID      beta        se     stat       pvalue n_features lambda       method
#> 1  a 0.9301722 0.1325913 7.015336 7.420512e-10          2   0.01 FE.empirical
#> 2  b 0.9005393 0.1284976 7.008220 6.457419e-09          2   0.01 FE.empirical
```
