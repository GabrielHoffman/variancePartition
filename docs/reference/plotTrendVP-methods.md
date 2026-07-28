# Plot Trend of Variance Fractions Vs Count Magnitude

Plot trend of variance fractions for a specified component versus count
magnitude for each gene and cell cluster

## Usage

``` r
plotTrendVP(x, vp, component, ...)

# S4 method for class 'DESeqDataSet,data.frame'
plotTrendVP(x, vp, component, ...)

# S4 method for class 'DGELRT,data.frame'
plotTrendVP(x, vp, component, ...)
```

## Arguments

- x:

  regression fits

- vp:

  `data.frame` from `fitVarPart()`

- component:

  variance component to extract from `vp`

- ...:

  additional arguments

## Value

Plot of variance fraction vs count magnitude

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
dds <- DESeq(dds, quiet=TRUE)
res <- results(dds)

# Variance partition analysis
vp1 <- varpart(dds)

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

# Plot count noise vs expression magnitude
plotTrendVP( fit, vp2, "CountNoise" )

```
