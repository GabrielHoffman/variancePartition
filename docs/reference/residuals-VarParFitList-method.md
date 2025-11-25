# Residuals from model fit

Extract residuals for each gene from model fit with fitVarPartModel()

## Usage

``` r
# S4 method for class 'VarParFitList'
residuals(object, ...)
```

## Arguments

- object:

  object produced by fitVarPartModel()

- ...:

  other arguments.

## Value

Residuals extracted from model fits stored in object

## Details

If model is fit with missing data, residuals returns NA for entries that
were missing in the original data

## Examples

``` r
# load library
# library(variancePartition)

library(BiocParallel)

# load simulated data:
# geneExpr: matrix of gene expression values
# info: information/metadata about each sample
data(varPartData)

# Specify variables to consider
# Age is continuous so we model it as a fixed effect
# Individual and Tissue are both categorical, so we model them as random effects
form <- ~ Age + (1 | Individual) + (1 | Tissue)

# Fit model
modelFit <- fitVarPartModel(geneExpr, form, info)

# Extract residuals of model fit
res <- residuals(modelFit)
```
