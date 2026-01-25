# Plot Variance Estimates

Plot Variance Estimates

## Usage

``` r
plotVarianceEstimates(
  fit,
  fitEB,
  var_true = NULL,
  xmax = quantile(fit$sigma^2, 0.999)
)
```

## Arguments

- fit:

  model fit from
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- fitEB:

  model fit from
  [`eBayes()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/eBayes.md)

- var_true:

  array of true variance values from simulation (optional)

- xmax:

  maximum value on the x-axis
