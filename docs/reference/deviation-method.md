# Deviation from expectation for each observation

Given a model fit for each features, residuals are computed and
transformed based on an absolute value or squaring transform.

## Usage

``` r
deviation(fit, method = c("AD", "SQ"), scale = c("leverage", "none"))

# S4 method for class 'MArrayLM'
deviation(fit, method = c("AD", "SQ"), scale = c("leverage", "none"))
```

## Arguments

- fit:

  model fit from
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- method:

  transform the residuals using absolute deviation ("AD") or squared
  deviation ("SQ").

- scale:

  scale each observation by "leverage", or no scaling ("none")

## Value

matrix of deviations from expection for each observation

## See also

[`diffVar()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/diffVar-method.md)

## Examples

``` r
# library(variancePartition)
library(edgeR)
data(varPartDEdata)

# filter genes by number of counts
isexpr <- rowSums(cpm(countMatrix) > 0.1) >= 5

# Standard usage of limma/voom
dge <- DGEList(countMatrix[isexpr, ])
dge <- calcNormFactors(dge)

# make this vignette faster by analyzing a subset of genes
dge <- dge[1:1000, ]

# regression formula
form <- ~Disease

# estimate precision weights
vobj <- voomWithDreamWeights(dge, form, metadata)

# fit dream model
fit <- dream(vobj, form, metadata)
fit <- eBayes(fit)

# Compute deviation from expection for each observation
# using model residuals
z <- deviation(fit)
z[1:4, 1:4]
#>                               sample_01 sample_02 sample_03 sample_04
#> ENST00000570099.1 gene=YPEL3 0.04040445 0.1318337 0.3916526 0.4099356
#> ENST00000589123.1 gene=NFIC  0.19896787 0.1255554 0.5928201 0.8993084
#> ENST00000360314.3 gene=CASS4 0.02940427 0.8317505 0.2532487 0.3961843
#> ENST00000456159.1 gene=MET   0.11576522 0.2188826 0.0129174 0.3682492
```
