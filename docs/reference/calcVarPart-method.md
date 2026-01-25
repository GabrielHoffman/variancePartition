# Compute variance statistics

Compute fraction of variation attributable to each variable in
regression model. Also interpretable as the intra-class correlation
after correcting for all other variables in the model.

## Usage

``` r
calcVarPart(fit, returnFractions = TRUE, ...)

# S4 method for class 'lm'
calcVarPart(fit, returnFractions = TRUE, ...)

# S4 method for class 'lmerMod'
calcVarPart(fit, returnFractions = TRUE, ...)

# S4 method for class 'glm'
calcVarPart(fit, returnFractions = TRUE, ...)

# S4 method for class 'negbin'
calcVarPart(fit, returnFractions = TRUE, ...)

# S4 method for class 'glmerMod'
calcVarPart(fit, returnFractions = TRUE, ...)
```

## Arguments

- fit:

  model fit from lm() or lmer()

- returnFractions:

  default: TRUE. If TRUE return fractions that sum to 1. Else return
  unscaled variance components.

- ...:

  additional arguments (not currently used)

## Value

fraction of variance explained / ICC for each variable in the regression
model

## Details

For linear model, variance fractions are computed based on the sum of
squares explained by each component. For the linear mixed model, the
variance fractions are computed by variance component estimates for
random effects and sum of squares for fixed effects.

For a generalized linear model, the variance fraction also includes the
contribution of the link function so that fractions are reported on the
linear (i.e. link) scale rather than the observed (i.e. response) scale.
For linear regression with an identity link, fractions are the same on
both scales. But for logit or probit links, the fractions are not well
defined on the observed scale due to the transformation imposed by the
link function.

The variance implied by the link function is the variance of the
corresponding distribution:

logit -\> logistic distribution -\> variance is \\\pi^\frac{2}{3}\\

probit -\> standard normal distribution -\> variance is 1

For the Poisson distribution with rate \\\lambda\\, the variance is
\\log(1 + 1/\lambda)\\.

For the negative binomial distribution with rate \\\lambda\\ and shape
\\\theta\\, the variance is \\log(1 + 1/\lambda + 1/\theta)\\.

Variance decomposition is reviewed by Nakagawa and Schielzeth (2012),
and expanded to other GLMs by Nakagawa, et al. (2017). See McKelvey and
Zavoina (1975) for early work on applying to GLMs. Also see DeMaris
(2002).

We note that Nagelkerke's pseudo R^2 evaluates the variance explained by
the full model. Instead, a variance partitioning approach evaluates the
variance explained by each term in the model, so that the sum of each
systematic plus random term sums to 1 (Hoffman and Schadt, 2016;
Nakagawa and Schielzeth, 2012).

## References

Nakagawa S, Johnson PC, Schielzeth H (2017). “The coefficient of
determination R 2 and intra-class correlation coefficient from
generalized linear mixed-effects models revisited and expanded.”
*Journal of the Royal Society Interface*, **14**(134), 20170213.  
  
Nakagawa S, Schielzeth H (2013). “A general and simple method for
obtaining R2 from generalized linear mixed-effects models.” *Methods in
ecology and evolution*, **4**(2), 133–142.  
  
McKelvey RD, Zavoina W (1975). “A statistical model for the analysis of
ordinal level dependent variables.” *Journal of mathematical sociology*,
**4**(1), 103–120.  
  
DeMaris A (2002). “Explained variance in logistic regression: A Monte
Carlo study of proposed measures.” *Sociological Methods & Research*,
**31**(1), 27–74.  
  
Hoffman GE, Schadt EE (2016). “variancePartition: interpreting drivers
of variation in complex gene expression studies.” *BMC bioinformatics*,
**17**(1), 1–13.  
  

## Examples

``` r
library(lme4)
data(varPartData)

# Linear mixed model
fit <- lmer(geneExpr[1, ] ~ (1 | Tissue) + Age, info)
calcVarPart(fit)
#>     Tissue        Age  Residuals 
#> 0.08672786 0.00093212 0.91234002 

# Linear model
# Note that the two models produce slightly different results
# This is expected: they are different statistical estimates
# of the same underlying value
fit <- lm(geneExpr[1, ] ~ Tissue + Age, info)
calcVarPart(fit)
#>     Tissue        Age  Residuals 
#> 0.08080435 0.00101471 0.91818094 
```
