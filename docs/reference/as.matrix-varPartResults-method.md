# Convert to matrix

Convert varPartResults to matrix

## Usage

``` r
# S4 method for class 'varPartResults'
as.matrix(x, ...)
```

## Arguments

- x:

  varPartResults

- ...:

  other arguments.

## Value

matrix

## Examples

``` r
# load library
# library(variancePartition)

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# Specify variables to consider
# Age is continuous so we model it as a fixed effect
# Individual and Tissue are both categorical, so we model them as random effects
form <- ~ Age + (1 | Individual) + (1 | Tissue)

# Fit model
varPart <- fitExtractVarPartModel(geneExpr[1:5, ], form, info)

# convert to matrix
as.matrix(varPart)
#>       Individual     Tissue          Age  Residuals
#> gene1  0.8903138 0.02468003 4.354754e-05 0.08496264
#> gene2  0.8060304 0.10102037 3.336677e-04 0.09261554
#> gene3  0.8899201 0.03630060 1.374661e-03 0.07240464
#> gene4  0.7688265 0.12531473 1.014416e-03 0.10484437
#> gene5  0.6997239 0.20910172 3.871483e-05 0.09113566
```
