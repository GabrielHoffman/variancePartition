# Error handling

Errors and warnings in `variancePartition` are mostly designed to let
the user know that there is an isssue with the model. Note that some of
these warnings and errors can be overridden by specifying
`hideErrorsInBackend=TRUE` for
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
and `showWarnings=FALSE` for
[`fitExtractVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitExtractVarPartModel-method.md)
and
[`fitVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitVarPartModel-method.md).

## Errors in `dream()`

The linear mixed model used by
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
can be a little fragile for small sample sizes and correlated
covariates.

- `Initial model failed: the fixed-effects model matrix is column rank deficient (rank(X) = 3 < 4 = p); the fixed effects will be jointly unidentifiable`

  The design matrix has redundant variables, so the model is singular
  and coefficients can’t be estimated. Fix by dropping one or more
  variables. Use
  [`canCorPairs()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/canCorPairs.md)
  to examine correlation betweeen variables.

### Gene-level errors

The most common issue is when
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
analysis succeeds for most genes, but a handful of genes fail. These
genes can fail if the iterative process of fitting the linear mixed
model does not converge, or if the estimated covariance matrix that is
supposed be positive definite has an eigen-value that is negative or too
close to zero due to rounding errors in floating point arithmetic.

In these cases,
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
gives a warning that the model has failed for a subset of genes, and
also provides the gene-level errors. All **successful** model fits are
returned to be used for downstream analysis.

Here we demonstrate how
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
handles model fits:

``` r

library(variancePartition)
data(varPartData)

# Redundant formula
# This example is an extreme example of redundancy
# but more subtle cases often show up in real data
form <- ~ Tissue + (1 | Tissue)

fit <- dream(geneExpr[1:30, ], form, info)

## Warning in dream(geneExpr[1:30, ], form, info): Model failed for 29 responses.
##   See errors with attr(., 'errors')

# Extract gene-level errors
attr(fit, "errors")[1:2]

## gene1
## "Error in lmerTest:::as_lmerModLT(model, devfun, tol = tol): (converted from warning)
## Model may not have converged with 1 eigenvalue close to zero: -2.0e-09\n"

## gene2
## "Error: (converted from warning) Model failed to converge
##   with 1 negative eigenvalue: -1.5e-08\n"
```

## Shared by multiple functions

These are shared by
[`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md),
[`fitVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitVarPartModel-method.md)
and
[`fitExtractVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitExtractVarPartModel-method.md).
Note that some of the these can be found in “1) Tutorial on using
variancePartition”.

### Warnings

- `No Intercept term was specified in the formula: The results will not behave as expected and may be very wrong!!`

  An intercept (i.e. mean term) must be specified order for the results
  to be statistically valid. Otherwise, the variance percentages will be
  *very* overestimated.

- `Categorical variables modeled as fixed effect: The results will not behave as expected and may be very wrong!!`

  If a linear mixed model is used, all categorical variables must be
  modeled as a random effect. Alternatively, a fixed effect model can be
  used by modeling all variables as fixed.

- `Cannot have more than one varying coefficient term:\newline The results will not behave as expected and may be very wrong!!`

  Only one varying coefficient term can be specified. For example, the
  formula `~(Tissue+0|Individual) + (Batch+0|Individual)` contains two
  varying coefficient terms and the results from this analysis are not
  easily interpretable. Only a formula with one term like
  `(Tissue+0|Individual)` is allowed.

### Errors

- `Colinear score > .99: Covariates in the formula are so strongly correlated that the parameter estimates from this model are not meaningful. Dropping one or more of the covariates will fix this problem`

- `Error in asMethod(object) : not a positive definite matrix`

- `In vcov.merMod(fit) : Computed variance-covariance matrix problem: not a positive definite matrix; returning NA matrix`

- `fixed-effect model matrix is rank deficient so dropping 26 columns / coefficients`

  Including variables that are highly correlated can produce misleading
  results (see Section “Detecting problems caused by collinearity of
  variables”). In this case, parameter estimates from this model are not
  meaningful. Dropping one or more of the covariates will fix this
  problem.

- `Error in checkNlevels(reTrms$flist, n = n, control): number of levels of each grouping factor must be < number of observations`

  This arises when using a varying coefficient model that examines the
  effect of one variable inside subsets of the data defined by another:
  `~(A+0|B)`. See Section “Variation within multiple subsets of the
  data”. There must be enough observations of each level of the variable
  B with each level of variable A. Consider an example with samples from
  multiple tissues from a set of individual where we are interested in
  the variation across individuals within each tissue using the formula:
  `~(Tissue+0|Individual)`. This analysis will only work if there are
  multiple samples from the same individual in at least one tissue. If
  all tissues only have one sample per individual, the analysis will
  fail and `variancePartition` will give this error.

- `Problem with varying coefficient model in formula: should have form (A+0|B)`

  When analyzing the variation of one variable inside another (see
  Section “Variation within multiple subsets of the data”.), the formula
  most be specified as `(Tissue+0|Individual)`. This error occurs when
  the formula contains `(Tissue|Individual)` instead.

- `fatal error in wrapper code`

- `Error in mcfork() : unable to fork, possible reason: Cannot allocate memory`

- `Error: cannot allocate buffer`

  This error occurs when `fitVarPartModel` uses too many threads and
  takes up too much memory. The easiest solution is to use
  `fitExtractVarPartModel` instead. Occasionally there is an issue in
  the parallel backend that is out of my control. Using fewer threads or
  restarting R will solve the problem.

#### Errors: Problems removing samples with NA/NaN/Inf values

`variancePartition` fits a regression model for each gene and drops
samples that have NA/NaN/Inf values in each model fit. This is generally
seamless but can cause an issue when a variable specified in the formula
no longer varies within the subset of samples that are retained.
Consider an example with variables for sex and age where age is NA for
all males samples. Dropping samples with invalid values for variables
included in the formula will retain only female samples. This will cause
`variancePartition` to throw an error because there is now no variation
in sex in the retained subset of the data. This can be resolved by
removing either age or sex from the formula.

This situtation is indicated by the following errors:

- `Error: grouping factors must have > 1 sampled level`

- `Error: Invalid grouping factor specification, Individual`

- `Error in contrasts<-(*tmp*, value = contr.funs[1 + isOF[nn]]): contrasts can be applied only to factors with 2 or more levels`

- `Error in checkNlevels(reTrms\$flist, n = n, control): grouping factors must have > 1 sampled level`

#### Errors with BiocParallel multithreading backend

- `Error: 'bpiterate' receive data failed: error reading from connection`

- `Error in serialize(data, node$con, xdr = FALSE) : ignoring SIGPIPE signal`

  `variancePartition` uses the `BiocParallel` package to run analysis in
  parallel across multiple cores. If there is an issue with the parallel
  backend you might see these errors. This often occurs in long
  interactive sessions, or if you manually kill a function running in
  parallel. There are two ways to address this issue.

  - **Global**: set the number of threads to be a smaller number. I have
    found that reducing the number of threads reduces the chance of
    random failures like this.

    ``` r

    library(BiocParallel)

    # globally specify that all multithreading using bpiterate from BiocParallel
    # should use 8 cores
    register(SnowParam(8))
    ```

  - **Local**: set the number of theads at each function call. This
    re-initializes the parallel backend and should address the error

    ``` r

    fitExtractVarPartModel(..., BPPARAM = SnowParam(8))

    fitVarPartModel(..., BPPARAM = SnowParam(8))

    dream(..., BPPARAM = SnowParam(8))

    voomWithDreamWeights(..., BPPARAM = SnowParam(8))
    ```

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.30/lib/libopenblasp-r0.3.30.dylib;  LAPACK version 3.12.0
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## time zone: America/New_York
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] digest_0.6.39     desc_1.4.3        R6_2.6.1          fastmap_1.2.0     xfun_0.54        
    ##  [6] cachem_1.1.0      knitr_1.50        htmltools_0.5.8.1 rmarkdown_2.30    lifecycle_1.0.4  
    ## [11] cli_3.6.5         sass_0.4.10       pkgdown_2.2.0     textshaping_1.0.4 jquerylib_0.1.4  
    ## [16] systemfonts_1.3.1 compiler_4.5.1    tools_4.5.1       ragg_1.5.0        bslib_0.9.0      
    ## [21] evaluate_1.0.5    yaml_2.3.10       jsonlite_2.0.0    rlang_1.1.6       fs_1.6.6         
    ## [26] htmlwidgets_1.6.4
