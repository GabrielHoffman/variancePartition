# Collinearity score

Collinearity score for a regression model indicating if variables are
too highly correlated to give meaningful results

## Usage

``` r
colinearityScore(fit)
```

## Arguments

- fit:

  regression model fit from lm() or lmer()

## Value

Returns the collinearity score between 0 and 1, where a score \> 0.999
means the degree of collinearity is too high. This function reports the
correlation matrix between coefficient estimates for fixed effects. The
collinearity score is the maximum absolute correlation value of this
matrix. Note that the values are the correlation between the parameter
estimates, and not between the variables themselves.

## Examples

``` r

# load library
# library(variancePartition)

# load simulated data:
data(varPartData)
#
form <- ~ Age + (1 | Individual) + (1 | Tissue)

res <- fitVarPartModel(geneExpr[1:10, ], form, info)

# evaluate the collinearity score on the first model fit
# this reports the correlation matrix between coefficients estimates
# for fixed effects
# the collinearity score is the maximum absolute correlation value
# If the collinearity score > .999 then the variance partition
# estimates may be problematic
# In that case, a least one variable should be omitted
colinearityScore(res[[1]])
#> [1] 0.7397006
#> attr(,"vcor")
#>             (Intercept)        Age
#> (Intercept)   1.0000000 -0.7397006
#> Age          -0.7397006  1.0000000
```
