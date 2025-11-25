# Effective sample size

Compute effective sample size based on correlation structure in linear
mixed model

## Usage

``` r
ESS(fit, method = "full")

# S4 method for class 'lmerMod'
ESS(fit, method = "full")
```

## Arguments

- fit:

  model fit from lmer()

- method:

  "full" uses the full correlation structure of the model. The
  "approximate" method makes the simplifying assumption that the study
  has a mean of m samples in each of k groups, and computes m based on
  the study design. When the study design is evenly balanced (i.e. the
  assumption is met), this gives the same results as the "full" method.

## Value

effective sample size for each random effect in the model

## Details

Effective sample size calculations are based on:

Liu, G., and Liang, K. Y. (1997). Sample size calculations for studies
with correlated observations. Biometrics, 53(3), 937-47.

"full" method: if \$\$V_x = var(Y;x)\$\$ is the variance-covariance
matrix of Y, the response, based on the covariate x, then the effective
sample size corresponding to this covariate is \$\$\Sigma\_{i,j}
(V_x^{-1})\_{i,j}\$\$. In R notation, this is: `sum(solve(V_x))`. In
practice, this can be evaluted as sum(w), where R

"approximate" method: Letting m be the mean number of samples per group,
\$\$k\$\$ be the number of groups, and \$\$\rho\$\$ be the intraclass
correlation, the effective sample size is \$\$mk / (1+\rho(m-1))\$\$

Note that these values are equal when there are exactly m samples in
each group. If m is only an average then this an approximation.

## Examples

``` r
library(lme4)
#> Loading required package: Matrix
data(varPartData)

# Linear mixed model
fit <- lmer(geneExpr[1, ] ~ (1 | Individual) + (1 | Tissue) + Age, info)

# Effective sample size
ESS(fit)
#> Individual     Tissue 
#>   27.24628   53.67295 
```
