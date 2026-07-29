# Plot Trend of Variance Fractions Vs Count Magnitude

Plot trend of variance fractions for a specified component versus count
magnitude for each gene and cell cluster

## Usage

``` r
plotTrendVP(x, vp, component, ...)

# S4 method for class 'DESeqDataSet,data.frame'
plotTrendVP(x, vp, component, ...)

# S4 method for class 'DGEGLM,data.frame'
plotTrendVP(x, vp, component, dispObj, ...)

# S4 method for class 'DGELRT,data.frame'
plotTrendVP(x, vp, component, dispObj, ...)
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

- dispObj:

  dispersion object if `edgeR` is used

## Value

Plot of variance fraction vs count magnitude

## Examples

``` r
# Simulate counts
set.seed(1)
eta <- rnorm(10, 3, 1)
mu <- exp(eta)

countMatrix <- matrix(rnbinom(n=100000, mu=mu, size=3), ncol=10)
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

vp2 <- varpart(fit, dispObj = d, formula = ~ cond)
#> Error in .local(fit, method, pseudocount, p.tail, ...): object 'x' not found

# Plot count noise vs expression magnitude
plotTrendVP( fit, vp2, "CountNoise", dispObj = d )
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'vp' in selecting a method for function 'plotTrendVP': object 'vp2' not found
```
