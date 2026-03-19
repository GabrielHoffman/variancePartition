# Theory and practice of random effects

Abstract

The distinction between modeling a variable as a fixed versus a random
effect depends on the goal of the statistical analysis. While some
theory and software make a strong distinction, `variancePartition` and
`dream` take different approaches based on the goal of each type of
analysis. Here we consider the distinction between fixed and random
effects, and the usage of REML in `variancePartition` and `dream`.

### `variancePartition`

#### Estimating contributions to expression variation

In traditional statistics and biostatistics, there is a strong
distinction between modeling categorical variants as fixed and random
effects. Random effects correspond to a sample of units from a larger
population, while fixed effects correspond to properties of specific
individuals. Random effects are typically treated as nuisance variables
and integrated out, and hypothesis testing is performed on the fixed
effect.

The `r2glmm` package fits into this traditional framework, by computing
the variance fractions for a given fixed effect as: \begin{eqnarray}
\sigma^2\_{fixed}/ \left(\sigma^2\_{fixed} + \sigma^2\_{error}\right)
\end{eqnarray}

Importantly, the random effects are not in the denominator. The fraction
is only determined by fixed effects and residuals.

In my experience in bioinformatics, this was a problem. Making such
distinctions between fixed and random effects seemed arbitrary. Variance
in a phenotype could be due to age (fixed) or to variation across
subject (random). Including all of the variables in the denominator
produced more intuitive results so that 1) the variance fractions sum to
one across all components and 2) fixed and random effects could be
interpreted on the same scale 3) fractions could be compared across
studies with different designs, 4) estimates of variance fractions were
most accurate. So in variancePartition the fractions are defined as:
\begin{eqnarray} \sigma^2\_{X}/ \left(\sigma^2\_{fixed} +
\sigma^2\_{random} + \sigma^2\_{error}\right) \end{eqnarray}

just plugging the each variable in the numerator.

Thus the faction evaluated by variancePartition is different than
`r2glmm` by definition.

Here is some code explicitly demonstrating this difference:

``` r

library("variancePartition")
library("lme4")
library("r2glmm")

set.seed(1)

N <- 1000
beta <- 3
alpha <- c(1, 5, 7)

# generate 1 fixed variable and 1 random variable with 3 levels
data <- data.frame(X = rnorm(N), Subject = sample(c("A", "B", "C"), 100, replace = TRUE))

# simulate variable
# y = X\beta + Subject\alpha + \sigma^2
data$y <- data$X * beta + model.matrix(~ data$Subject) %*% alpha + rnorm(N, 0, 1)

# fit model
fit <- lmer(y ~ X + (1 | Subject), data, REML = FALSE)

# calculate variance fraction using variancePartition
# include the total sum in the denominator
frac <- calcVarPart(fit)
frac
```

      Subject         X Residuals 
       0.4505    0.4952    0.0543 

``` r

# the variance fraction excluding the random effect from the denominator
# is the same as from r2glmm
frac[["X"]] / (frac[["X"]] + frac[["Residuals"]])
```

    [1] 0.901

``` r

# using r2glmm
r2beta(fit)
```

      Effect   Rsq upper.CL lower.CL
    1  Model 0.896    0.904    0.886
    2      X 0.896    0.904    0.886

So the formulas are different. But why require categorical variables as
random effects?

At practical level, categorical variables with too many levels are
problematic. Using a categorical variable with 200 categories as a fixed
effect is statistically unstable. There are so many degrees of freedom
that that variable will absorb a lot of variance even under the null.
Statistically, estimating the variance fraction for a variable with many
categories can be biased if that variable is a fixed effect. Therefore,
`variancePartition` requires all categorical variables to be random
effects. Modeling this variable as a random effect produces unbiased
estimates of variance fractions in practice. See simulations in the
Supplement (section 1.5) of [Hoffman and Schadt
(2016)](https://doi.org/10.1186/s12859-016-1323-z).

The distinction between fixed and random effects is important in the
formulation because it affects which variables are put in the
denominator. So choosing to model a variable as a fixed versus random
effect will definitely change the estimated fraction.

Yet for the `variancePartition` formulation, all variables are in the
denominator and it isn\`t affected by the fixed/random decision.
Moreover, using a random effect empirically reduces the bias of the
estimated fraction.

Finally, why use maximum likelihood to estimate the paramters instead of
the default REML ()? Maximum likelihood fits all parameters jointly so
that it estimates the fixed and random effects together. This is
essential if we want to compare fixed and random effects later.
Conversely, REML estimates the random effects by removing the fixed
effects from the response before estimation. This implicitly removes the
fixed effects from the denominator when evaluating the variance
fraction. REML treats fixed effects as nuisance variables, while
`variancePartition` considers fixed effects to be a core part of the
analysis.

While REML produced unbiased estimates of the variance components, the
goal of `variancePartition` is to estimate the variance fractions for
fixed and random effects jointly. In simulations from the Supplement
(section 1.5) of [Hoffman and Schadt
(2016)](https://doi.org/10.1186/s12859-016-1323-z), REML produced biased
estimates of the variance fractions while maximum likelihood estimates
are unbiased.

### `dream`

#### Hypothesis testing

While `dream` is also based on a linear mixed model, the goal of this
analysis is to perform hypothesis testing on fixed effects. Random
effects are treated as nuisance variables to be integrated out, and the
approximate null distribution of a t- or F-statistic is constructed from
the model fit.

Since the goal of the analysis is different, the consideration of using
REML versus ML is different than above. While is required by called by
when , can be used with as either or . Since the Kenward-Roger method
gave the best power with an accurate control of false positive rate in
our simulations, and since the Satterthwaite method with gives p-values
that are slightly closer to the Kenward-Roger p-values, is set as the
default.

## Session Info

    R version 4.5.1 (2025-06-13)
    Platform: aarch64-apple-darwin23.6.0
    Running under: macOS Sonoma 14.7.1

    Matrix products: default
    BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.31_1/lib/libopenblasp-r0.3.31.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] knitr_1.51

    loaded via a namespace (and not attached):
     [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0    
     [5] xfun_0.56         cachem_1.1.0      htmltools_0.5.9   rmarkdown_2.30   
     [9] lifecycle_1.0.5   cli_3.6.5         sass_0.4.10       pkgdown_2.2.0    
    [13] textshaping_1.0.5 jquerylib_0.1.4   systemfonts_1.3.2 compiler_4.5.1   
    [17] tools_4.5.1       ragg_1.5.1        bslib_0.10.0      evaluate_1.0.5   
    [21] yaml_2.3.12       otel_0.2.0        jsonlite_2.0.0    rlang_1.1.7      
    [25] fs_1.6.7          htmlwidgets_1.6.4
