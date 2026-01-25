# Extract contrast matrix for linear mixed model

Extract contrast matrix, L, testing a single variable. Contrasts
involving more than one variable can be constructed by modifying L
directly

## Usage

``` r
getContrast(exprObj, formula, data, coefficient)
```

## Arguments

- exprObj:

  matrix of expression data (g genes x n samples), or `ExpressionSet`,
  or `EList` returned by
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html) from the `limma`
  package

- formula:

  specifies variables for the linear (mixed) model. Must only specify
  covariates, since the rows of exprObj are automatically used as a
  response. e.g.: `~ a + b + (1|c)` Formulas with only fixed effects
  also work

- data:

  data.frame with columns corresponding to formula

- coefficient:

  the coefficient to use in the hypothesis test

## Value

Contrast matrix testing one variable

## Examples

``` r

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# get contrast matrix testing if the coefficient for Batch2 is zero
# The variable of interest must be a fixed effect
form <- ~ Batch + (1 | Individual) + (1 | Tissue)
L <- getContrast(geneExpr, form, info, "Batch3")

# get contrast matrix testing if Batch3 - Batch2 = 0
form <- ~ Batch + (1 | Individual) + (1 | Tissue)
L <- getContrast(geneExpr, form, info, c("Batch3", "Batch2"))

# To test against Batch1 use the formula:
#   ~ 0 + Batch + (1|Individual) + (1|Tissue)
# to estimate Batch1 directly instead of using it as the baseline
```
