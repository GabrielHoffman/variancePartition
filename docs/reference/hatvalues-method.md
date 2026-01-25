# Compute hatvalues

Compute hatvalues from dream fit

## Usage

``` r
# S4 method for class 'MArrayLM'
hatvalues(model, vobj, ...)

# S4 method for class 'MArrayLM2'
hatvalues(model, ...)
```

## Arguments

- model:

  model fit from
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- vobj:

  `EList` returned by
  [`voom()`](https://rdrr.io/pkg/limma/man/voom.html) or
  [`voomWithDreamWeights()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/voomWithDreamWeights.md).

- ...:

  other arguments, currently ignored
