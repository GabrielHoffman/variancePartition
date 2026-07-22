# Convert to data.frame

Convert varPartResults to data.frame

## Usage

``` r
# S3 method for class 'varPartResults'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  varPartResults

- row.names:

  pass thru to generic

- optional:

  pass thru to generic

- ...:

  other arguments.

## Value

data.frame

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
as.data.frame(varPart)
#>                Age Individual     Tissue  Residuals
#> gene1 4.403094e-05  0.8911946 0.02470444 0.08405689
#> gene2 3.342172e-04  0.7992842 0.10017486 0.10020673
#> gene3 1.378790e-03  0.8836673 0.03604554 0.07890839
#> gene4 1.022546e-03  0.7672383 0.12505585 0.10668333
#> gene5 3.876328e-05  0.6935935 0.20726975 0.09909796
```
