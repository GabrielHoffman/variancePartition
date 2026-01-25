# Augment observed read counts with prior count

Augment observed read counts with prior count since log of zero counts
is undefined. The prior count added to each sample is scaled so that no
variance is introduced

## Usage

``` r
augmentPriorCount(
  counts,
  lib.size = colSums2(counts),
  prior.count = 0.5,
  scaledByLib = FALSE
)
```

## Arguments

- counts:

  matrix of read counts with genes as rows and samples as columns

- lib.size:

  library sizes, the sum of all ready for each sample

- prior.count:

  average prior count added to each sample.

- scaledByLib:

  if `TRUE`, scale pseudocount by `lib.size`. Else to standard constant
  pseudocount addition

## Value

matrix with augmented counts

## Details

Adding prior counts removes the issue of evaluating the log of zero
counts, and stabilizes the log transform when counts are very small.
However, adding a constant prior count to all samples can introduced an
artifact. Consider two samples each with zero counts for a given gene,
but one as a library size of 1k and the other of 50k. After applying the
prior count values become pc / 1k and pc / 50k. It appears that there is
variance in the expression of this gene, even though no counts are
observed. This is driven only by variation in the library size, which
does not reflect biology. This issue is most problematic for small
counts.

Instead, we make the reasonable assumption that a gene does not have
expression variance unless supported sufficiently by counts in the
numerator. Consider adding a different prior count to each sample so
that genes with zero counts end up woth zero variance. This corresponds
to adding `prior.count * lib.size[i] / mean(lib.size)` to sample `i`.

This is done in the backend of
[`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html), but this
function allows users to apply it more generally.

## See also

[`edgeR::cpm()`](https://rdrr.io/pkg/edgeR/man/cpm.html)

## Examples

``` r
library(edgeR)

data(varPartDEdata)

# normalize RNA-seq counts
dge <- DGEList(counts = countMatrix)
dge <- calcNormFactors(dge)

countsAugmented <- augmentPriorCount( dge$counts, dge$samples$lib.size, 1)
```
