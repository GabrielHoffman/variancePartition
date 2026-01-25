# Table of Top Genes from Linear Model Fit

topTable generic

topTable generic MArrayLM

topTable generic MArrayLM2

## Usage

``` r
topTable(
  fit,
  coef = NULL,
  number = 10,
  genelist = fit$genes,
  adjust.method = "BH",
  sort.by = "B",
  resort.by = NULL,
  p.value = 1,
  lfc = 0,
  confint = FALSE
)

# S4 method for class 'MArrayLM'
topTable(
  fit,
  coef = NULL,
  number = 10,
  genelist = fit$genes,
  adjust.method = "BH",
  sort.by = "p",
  resort.by = NULL,
  p.value = 1,
  lfc = 0,
  confint = FALSE
)

# S4 method for class 'MArrayLM2'
topTable(
  fit,
  coef = NULL,
  number = 10,
  genelist = fit$genes,
  adjust.method = "BH",
  sort.by = "p",
  resort.by = NULL,
  p.value = 1,
  lfc = 0,
  confint = FALSE
)
```

## Arguments

- fit:

  fit

- coef:

  coef

- number:

  number

- genelist:

  genelist

- adjust.method:

  adjust.method

- sort.by:

  sort.by

- resort.by:

  resort.by

- p.value:

  p.value

- lfc:

  lfc

- confint:

  confint

## Value

results of toptable

results of toptable

results of toptable
