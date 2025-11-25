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
    BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/New_York
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
    [1] r2glmm_0.1.3             lme4_1.1-38              Matrix_1.7-4            
    [4] variancePartition_1.39.5 BiocParallel_1.42.2      limma_3.64.3            
    [7] ggplot2_4.0.1            knitr_1.50              

    loaded via a namespace (and not attached):
     [1] tidyselect_1.2.1    dplyr_1.1.4         farver_2.1.2       
     [4] S7_0.2.1            bitops_1.0-9        fastmap_1.2.0      
     [7] digest_0.6.39       lifecycle_1.0.4     statmod_1.5.1      
    [10] magrittr_2.0.4      compiler_4.5.1      rlang_1.1.6        
    [13] sass_0.4.10         tools_4.5.1         yaml_2.3.10        
    [16] htmlwidgets_1.6.4   plyr_1.8.9          RColorBrewer_1.1-3 
    [19] KernSmooth_2.23-26  withr_3.0.2         purrr_1.2.0        
    [22] numDeriv_2016.8-1.1 BiocGenerics_0.54.1 desc_1.4.3         
    [25] grid_4.5.1          aod_1.3.3           caTools_1.18.3     
    [28] scales_1.4.0        gtools_3.9.5        iterators_1.0.14   
    [31] MASS_7.3-65         dichromat_2.0-0.1   cli_3.6.5          
    [34] mvtnorm_1.3-3       rmarkdown_2.30      ragg_1.5.0         
    [37] reformulas_0.4.2    generics_0.1.4      reshape2_1.4.5     
    [40] minqa_1.2.8         cachem_1.1.0        stringr_1.6.0      
    [43] splines_4.5.1       parallel_4.5.1      matrixStats_1.5.0  
    [46] vctrs_0.6.5         boot_1.3-32         jsonlite_2.0.0     
    [49] pbkrtest_0.5.5      systemfonts_1.3.1   jquerylib_0.1.4    
    [52] tidyr_1.3.1         glue_1.8.0          nloptr_2.2.1       
    [55] pkgdown_2.2.0       codetools_0.2-20    stringi_1.8.7      
    [58] gtable_0.3.6        EnvStats_3.1.0      lmerTest_3.1-3     
    [61] tibble_3.3.0        remaCor_0.0.20      pillar_1.11.1      
    [64] htmltools_0.5.8.1   gplots_3.2.0        R6_2.6.1           
    [67] textshaping_1.0.4   Rdpack_2.6.4        evaluate_1.0.5     
    [70] lattice_0.22-7      Biobase_2.68.0      rbibutils_2.4      
    [73] backports_1.5.0     RhpcBLASctl_0.23-42 broom_1.0.10       
    [76] fANCOVA_0.6-1       corpcor_1.6.10      bslib_0.9.0        
    [79] Rcpp_1.1.0          nlme_3.1-168        xfun_0.54          
    [82] fs_1.6.6            pkgconfig_2.0.3    
