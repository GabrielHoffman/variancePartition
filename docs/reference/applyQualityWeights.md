# Apply pre-specified sample weights

Apply pre-specified sample weights by scaling existing precision weights

## Usage

``` r
applyQualityWeights(vobj, weights)
```

## Arguments

- vobj:

  `EList` from `voom` or `voomWithDreamWeights`.

- weights:

  sample level weights

## Details

Apply pre-specified sample-level weights to the existing precision
weights estimated from the data. While the
[`limma::voomWithQualityWeights`](https://rdrr.io/pkg/limma/man/voomWithQualityWeights.html)
function of Lui et al. (2015) estimates the sample-level weights from
`voom` fit, here the weights are fixed beforehand.

## References

Liu R, Holik AZ, Su S, Jansz N, Chen K, Leong HS, Blewitt ME,
Asselin-Labat M, Smyth GK, Ritchie ME (2015). “Why weight? Modelling
sample and observational level variability improves power in RNA-seq
analyses.” *Nucleic acids research*, **43**(15), e97–e97.

## See also

[`limma::voomWithQualityWeights`](https://rdrr.io/pkg/limma/man/voomWithQualityWeights.html)
