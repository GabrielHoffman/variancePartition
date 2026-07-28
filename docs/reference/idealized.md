# Idealized Variance Fractions

Compute idealized variance fractions removing `CountNoise` or other
variables

## Usage

``` r
idealized(x, remove = "CountNoise")
```

## Arguments

- x:

  `matrix` or `data.frame` of variance fracdtions

- remove:

  colnames of `x` to remove from the variance fractions

## Details

Variance partitioning analysis includes all components of variance in
the denominator when computing the variance fractions. This function
recomputes the variance fractions in idealized case where `CountNoise`
or any set of variances given in the `remove` argument is zero.

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
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following object is masked from ‘package:lme4’:
#> 
#>     refit
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following object is masked from ‘package:limma’:
#> 
#>     plotMA
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> 
#> Attaching package: ‘S4Vectors’
#> The following objects are masked from ‘package:Matrix’:
#> 
#>     expand, unname
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: GenomicRanges
#> Loading required package: Seqinfo
#> Loading required package: SummarizedExperiment
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
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

# Compute variance fractions in the idealized
#  case where there is no count noise
head(idealized(vp1))
#>           condition Residuals
#> gene_1 0.3545538734 0.6454461
#> gene_2 0.0001032731 0.9998967
#> gene_3 0.4544569666 0.5455430
#> gene_4 0.0264197616 0.9735802
#> gene_5 0.1288757222 0.8711243
#> gene_6 0.0010682706 0.9989317
```
