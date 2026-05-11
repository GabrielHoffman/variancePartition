# Variance Partitioning Analysis

Variance partitioning analysis on each gene

## Usage

``` r
varpart(
  x,
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  ...
)

# S4 method for class 'DESeqDataSet'
varpart(
  x,
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  ...
)

# S4 method for class 'DGELRT'
varpart(
  x,
  method = c("exact", "approximate"),
  pseudocount = 1,
  p.tail = 1e-04,
  nthreads = parallelly::availableCores(),
  dispObj,
  formula,
  ...
)
```

## Arguments

- x:

  regression model fit

- method:

  select method for count models: `"exact"` or `"approximate"` for
  faster approximation

- pseudocount:

  pseudocount used for `"exact"` and `"approximate"` methods for count
  models

- p.tail:

  probability threashold for evaluating expectations for `"exact"`
  methods for count models

- nthreads:

  number of threads used for count models

- ...:

  other arguments

- dispObj:

  result of
  [`estimateDisp()`](https://rdrr.io/pkg/edgeR/man/estimateDisp.html)

- formula:

  formula used for the design matrix

## Examples

``` r
# Simulate counts
set.seed(1)
countMatrix <- matrix(rnbinom(n=100000, mu=20, size=3), ncol=10)
rownames(countMatrix) <- paste0("gene_", seq(nrow(countMatrix)))
colnames(countMatrix) <- paste0("sample_", seq(ncol(countMatrix)))

condition <- factor(rep(1:2, each=5))

# DESeq2 model #
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countMatrix, 
  DataFrame(condition), 
  ~ condition)
#> converting counts to integer mode
dds <- DESeq(dds)
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
res <- results(dds)

# Variance partition analysis
vp1 <- varpart(dds)

# Plot contribution of each component
plotVarPart(vp1, main="DESeq2")


# Plot count noise vs expression magnitude
plotTrendVP( dds, vp1, "CountNoise" )



# edgeR model #
library(edgeR)
design <- model.matrix( ~ condition, data.frame(condition))
d <- DGEList(countMatrix)
d <- normLibSizes(d)
d <- estimateDisp(d, design)
fit <- glmQLFit(d, design)
fit <- glmQLFTest(fit)

vp2 <- varpart(fit, dispObj = d, formula = ~ cond)

# Plot contribution of each component
plotVarPart(vp2, main="edgeR")


# Plot count noise vs expression magnitude
plotTrendVP( dds, vp2, "CountNoise" )
#> Warning: Failed to fit group -1.
#> Caused by error in `nls()`:
#> ! step factor 0.000488281 reduced below 'minFactor' of 0.000976562

```
