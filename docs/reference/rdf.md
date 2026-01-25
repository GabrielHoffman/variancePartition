# Residual degrees of freedom

Residual degrees of freedom

## Usage

``` r
rdf(fit)
```

## Arguments

- fit:

  model fit from [`lm()`](https://rdrr.io/r/stats/lm.html),
  [`glm()`](https://rdrr.io/r/stats/glm.html),
  [`lmer()`](https://rdrr.io/pkg/lme4/man/lmer.html)

## See also

`rdf.merMod`

## Examples

``` r
library(lme4)

fit <- lm(Reaction ~ Days, sleepstudy)
rdf(fit)
#> [1] 178
```
