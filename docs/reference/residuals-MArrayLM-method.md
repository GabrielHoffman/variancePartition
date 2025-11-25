# residuals for MArrayLM

residuals for MArrayLM

## Usage

``` r
# S4 method for class 'MArrayLM'
residuals(object, y, ..., type = c("response", "pearson"))
```

## Arguments

- object:

  MArrayLM object from dream

- y:

  `EList` object used in
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- ...:

  other arguments, currently ignored

- type:

  compute either response or pearson residuals

## Value

results of residuals
