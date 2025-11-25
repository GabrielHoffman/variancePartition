# Compute standard post-processing values

These values are typically computed by eBayes

## Usage

``` r
.standard_transform(fit, sigma = fit$sigma)
```

## Arguments

- fit:

  result of dream (MArrayLM2)

- sigma:

  vector of standard errors used to compute t-statistic. Can be maximum
  likelihood estimates, or posterior means

## Value

MArrayLM2 object with values computed
