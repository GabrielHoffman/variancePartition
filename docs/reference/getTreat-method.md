# Test if coefficient is different from a specified value

Test if coefficient is different from a specified value

## Usage

``` r
getTreat(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p")

# S4 method for class 'MArrayLM'
getTreat(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p")

# S4 method for class 'MArrayLM2'
getTreat(fit, lfc = log2(1.2), coef = 1, number = 10, sort.by = "p")
```

## Arguments

- fit:

  fit

- lfc:

  a minimum log2-fold-change below which changes not considered
  scientifically meaningful

- coef:

  which coefficient to test

- number:

  number of genes to return

- sort.by:

  column to sort by

## Value

results of getTreat

## Examples

``` r

data(varPartData)

form <- ~ Age + Batch + (1 | Individual) + (1 | Tissue)

fit <- dream(geneExpr, form, info)
fit <- eBayes(fit)

coef <- "Age"

# Evaluate treat()/topTreat() in a way that works seamlessly for dream()
getTreat(fit, lfc = log2(1.03), coef, sort.by = "none", number = 3)
#>             logFC     AveExpr t   P.Value adj.P.Val         B z.std
#> gene1 0.005551165 -10.4664549 0 0.9882946 0.9997854 -6.323485     0
#> gene2 0.013121135  -1.1281610 0 0.9586703 0.9997854 -6.084858     0
#> gene3 0.020151494   0.1702122 0 0.9211006 0.9997854 -5.586873     0
```
