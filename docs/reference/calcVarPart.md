# Compute variance statistics

Compute fraction of variation attributable to each variable in
regression model. Also interpretable as the intra-class correlation
after correcting for all other variables in the model.

## Usage

``` r
calcVarPart(fit, ...)
```

## Arguments

- fit:

  regression model

- ...:

  other arguments passed to
  [`fastglmm::varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)

## Details

this is now an iterface to
[`fastglmm::varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)

## See also

[`fastglmm::varpart()`](http://gabrielhoffman.github.io/fastglmm/reference/varpart.md)

## Examples

``` r
library(lme4)
data(varPartData)

# Linear mixed model
fit <- lmer(geneExpr[1, ] ~ (1 | Tissue) + Age, info)
calcVarPart(fit)
#>          Age       Tissue    Residuals 
#> 0.0005195716 0.0478594193 0.9516210091 

# Linear model
# Note that the two models produce slightly different results
# This is expected: they are different statistical estimates
# of the same underlying value
fit <- lm(geneExpr[1, ] ~ Tissue + Age, info)
calcVarPart(fit)
#>       Tissue          Age    Residuals 
#> 0.0784058883 0.0009845912 0.9206095205 
```
