# Shrinkage metric for eBayes

Evaluates the coefficient from the linear regression of
`s2.post ~ sigmaSq`. When there is no shrinkage, this value is 1. Values
less than 1 indicate the amount of shrinkage.

## Usage

``` r
shrinkageMetric(sigmaSq, s2.post)
```

## Arguments

- sigmaSq:

  maximum likelihood residual variance for every gene

- s2.post:

  empirical Bayes posterior estimate of residual variance for every gene

## Details

Shrinkage metric for eBayes quantifying the amount of shrinkage that is
applied to shrink the maximum likelihood residual variance to the
empirical Bayes posterior estimate
