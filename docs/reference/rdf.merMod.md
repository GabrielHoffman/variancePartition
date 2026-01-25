# Approximate residual degrees of freedom

For a linear model with \\n\\ samples and \\p\\ covariates,
\\RSS/\sigma^2 \sim \chi^2\_{\nu}\\ where \\\nu = n-p\\ is the residual
degrees of freedom. In the case of a linear mixed model, the
distribution is no longer exactly a chi-square distribution, but can be
approximated with a chi-square distribution.

Given the hat matrix, \\H\\, that maps between observed and fitted
responses, the approximate residual degrees of freedom is \\\nu =
tr((I-H)^T(I-H))\\. For a linear model, this simplifies to the well
known form \\\nu = n - p\\. In the more general case, such as a linear
mixed model, the original form simplifies only to \\n - 2tr(H) +
tr(HH)\\ and is an approximation rather than being exact. The third term
here is quadratic time in the number of samples, \\n\\, and can be
computationally expensive to evaluate for larger datasets. Here we
develop a linear time algorithm that takes advantage of the fact that
\\H\\ is low rank.

\\H\\ is computed as \\A^TA + B^TB\\ for `A=CL` and `B=CR` defined in
the code. Since \\A\\ and \\B\\ are low rank, there is no need to
compute \\H\\ directly. Instead, the terms \\tr(H)\\ and \\tr(HH)\\ can
be computed using the eigen decompositions of \\AA^T\\ and \\BB^T\\
which is linear time in the number of samples.

## Usage

``` r
rdf.merMod(model, method = c("linear", "quadratic"))
```

## Arguments

- model:

  An object of class `merMod`

- method:

  Use algorithm that is "linear" (default) or quadratic time in the
  number of samples

## Value

residual degrees of freedom

## Details

Compute the approximate residual degrees of freedom from a linear mixed
model.

## See also

rdf_from_matrices

## Examples

``` r
library(lme4)

# Fit linear mixed model
fit <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)

# Evaluate the approximate residual degrees of freedom
rdf.merMod(fit)
#> [1] 146.3644
```
