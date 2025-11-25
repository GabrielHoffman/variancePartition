# Get all univariate contrasts

Get all univariate contrasts

## Usage

``` r
.getAllUniContrasts(formula, data)
```

## Arguments

- formula:

  specifies variables for the linear (mixed) model. Must only specify
  covariates, since the rows of exprObj are automatically used as a
  response. e.g.: `~ a + b + (1|c)` Formulas with only fixed effects
  also work

- data:

  data.frame with columns corresponding to formula

## Value

Matrix testing each variable one at a time. Contrasts are on rows
