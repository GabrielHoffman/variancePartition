
# variancePartition 1.39.5
  - Nov 25, 2025
  - move calls to `nobars()` and `findbars()` to `reformulas` package

# variancePartition 1.39.4
  - handle `NA` values for `fitExtractVarPartModel()`

# variancePartition 1.39.2
  - update checks to handle `voom(...,adaptive.span)`

# variancePartition 1.39.1
  - version bump

# variancePartition 1.38.1
 - Aug 11, 2025
 - Handle change in function prototype of `limma::eBayes()`
  - `variancePartition::eBayes()` now uses arguments `...`

# variancePartition 1.37.4
 - Aug 4, 2025
 - Handle missing expression values like in `limma::lmFit()`

# variancePartition 1.37.3
 - June 22, 2025
 - update docs

# variancePartition 1.36.3
 - Jan 2, 2025
 - examine adding `REDUCE = c` to `bpiterate()`
 - fix issue `mvTest()` for test of single coefficient

# variancePartition 1.36.2
 - Nov 7, 2024
 - fix docs

# variancePartition 1.36.1
 - Nov 5, 2024
 - Upgrade to Bioc 3.20
 - update generic for `eBayes()` to be comptabile with `legacy` parameter added in `limma` v3.62.0 

# variancePartition 1.35.5
 - June 17, 2024
 - fix bug in `mvTest()` for `se` value with 1 feature

# variancePartition 1.35.4
 - June 13, 2024
 - fix bug in `[.MArrayLM2` dropping residuals to a vector

# variancePartition 1.35.3
 - June 11, 2024
 - in `mvTest()` with `shrink.cov = TRUE` uses `lambda = 0.01`

# variancePartition 1.35.2
 - June 6, 2024
 - in `mvTest()` return `beta` and `se`

# variancePartition 1.35.1
 - May 31, 2024
 - version bump for BioC

# variancePartition 1.33.15
 - May 16, 2024
 - in `voomWithDreamWeights()`, default changed to `rescaleWeightsAfter = FALSE`

# variancePartition 1.33.14
 - April 19, 2024
 - in `topTable()` resolve issue when specifying multiple coeffs

# variancePartition 1.33.13
 - April 6, 2024
 - in `voomWithDreamWeights()`, add argument  `prior.count.for.weights`

# variancePartition 1.33.12
 - March 22, 2024
 - Fix use with `voomWithDreamWeights()` taking raw counts
  - https://github.com/GabrielHoffman/variancePartition/issues/97
 - in `voomWithDreamWeights()`, add argument `priorWeightsAsCounts=FALSE` and `prior.count.for.weights`

# variancePartition 1.33.11
 - Feb 7, 2024
 - fix bug in `dream(...,ddf="Kenward-Roger")` that gave false positives and negatives
   - scaled weights properly to fix this issue, since `df` in `lmerTest::contest()` assumes mean of weights is 1.
   - old code used unscaled weights, so df values were too large

# variancePartition 1.33.10
 - Feb 5, 2024
 - in `augmentPriorCount()` and `voomWithDreamWeights()`, add argument `scaledByLib=FALSE`

# variancePartition 1.33.9
 - Jan 25, 2024
 - bug fix in `BIC()` and `.fitExtractVarPartModel()`

# variancePartition 1.33.8
 - Jan 22, 2024
 - add `BIC()` for result of `dream()` and `lmFit()`

# variancePartition 1.33.7
 - Jan 19, 2024
 - add `logLik()` for result of `dream()` and `lmFit()`

# variancePartition 1.33.6
 - fix issue with subsetting using logicals

# variancePartition 1.33.5
 - Dec 21, 2023
 - fix issue with ddf always calling `"Satterthwaite"`
 - enforce subsetting of residuals in `assign( "[.MArrayLM2",)`
   - this comes with when `rdf < 1`

# variancePartition 1.33.4
 - Dec 13, 2023
 - fix bug with BiocParallel in Windows
   - https://github.com/GabrielHoffman/variancePartition/issues/91
 - handle case where a single contrasts fails in `makeContrastsDream()`

# variancePartition 1.33.3
 - Dec 5, 2023
 - in `residuals()` when dividing by `sqrt(1-hatvalues)` add small offset to make sure the value is positive

# variancePartition 1.33.2
 - Nov 13, 2023
 - add `augmentPriorCount()`
 - add `prior.count` argument to `voomWithDreamWeights()` and feed it to `augmentPriorCount()`

# variancePartition 1.33.1
 - Nov 2, 2023
 - move `omp_set_num_threads()` deeper in nesting

# variancePartition 1.31.22
 - Oct 19, 2023
 - fix handling of variables with missing data
 - return `fit$genes` properly

# variancePartition 1.31.21
 - Oct 16, 2023
 - handle weights properly when the linear mixed model fails for some genes
 - `lmFit()` and 
  - in `iterRows()` set `scale = FALSE` as default
  - in `voomWithDreamWeights()`, scale in input weights and weights in side `fitVarPartModel()`
    - get weights estimated most similar to `voomLmFit()`
  - in `dream()` use `rescaleWeights = FALSE` to get `sigma` estimates compatable with `lmFit()`

# variancePartition 1.31.20
 - Sept 26, 2023
 - allow `weights` to be a matrix in `voomWithDreamWeights()`

# variancePartition 1.31.19
 - Sept 22, 2023
 - add `rescaleWeightsAfter` argument to `voomWithDreamWeights()`

# variancePartition 1.31.18
 - Sept 5, 2023
 - improved error handling for `fitVarPartModel()`, `fitExtractVarPartModel()`, and `voomWithDreamWeights()`

# variancePartition 1.31.16
 - August 18, 2023
 - in `dream()`, if `"Kenward-Roger"` is specified but gives covariance matrix that has poor condition number or is not positive definite, then fall back to `"Satterthwaite"` for hypothesis testing in linear mixed models
 - Update documentation, and reformat code

# variancePartition 1.31.15
 - August 10, 2023
 - `fit = dream()` now returns `fit$loglik` (the log-likelihood for each gene), and `fit$edf` (the effective degreees of freedom for each gene)

# variancePartition 1.31.13
 - August 7, 2023
 - fix bug in `calcVarPart()` where weights was ignored in some cases
  - add additional tests to check this

# variancePartition 1.31.12
 - July 3, 2023
 - `makeContrastsDream()` converts `NA` `contrasts` to `NULL`

# variancePartition 1.31.11
 - setting `voomWithDreamWeights(..., span="auto")` now estimates tuning parameter from data using `fANCOVA::loess.as()`

# variancePartition 1.31.10
 - `filterInputData()` now ensures EList contains a matrix

# variancePartition 1.31.9
 - Fix issue in `mvTest()` when specifying features with strings

# variancePartition 1.31.8
 - Fix error message when linear mixed model fails

# variancePartition 1.31.7
 - update `mvTest()` to run in parallel

# variancePartition 1.31.6
 - update `mvTest()` to include Hotelling T2 test and `LS.empirical()`

# variancePartition 1.31.1
 - June 4, 2023
 - Rename for Bioconductor compatability

# variancePartition 2.0.5
 - May 31, 2023
 - fix convergence issues
 - fix initialization of `lmer()` fit
 - use 1 OMP thread internally, then restore to original value

# variancePartition 2.0.4
 - May 30, 2023
 - When running `dream()`, ensure model convergence using second fitting with `Nelder_Mead` to avoid edge cases where the approximate hessian from `lmerTest::as_lmerModLT()` has a negative eigenvalue
 - fix issue in `get_prediction()` returning NA values when variables modeled as categorical and levels are omitted
 - fix issue in `voomWithDreamWeights()` when some genes don't converge
- retry `lmer()` model fit with another optimizer after it fails convergence test.

# variancePartition 2.0.3
 - May 13, 2023
 - fix `vcov()`

# variancePartition 2.0.2
 - May 17, 2023
 - add matrix argument to `mvTest()`

# variancePartition 2.0.1
 - May 12, 2023
 - `mvTest()` now shrinks covariance using the Schafer-Strimmer method
 - `vcovSqrt()` returns the matrix whose cross product gives the `vcov()` result from fits with `dream()`

# variancePartition 2.0.0
 - April 20, 2023
 - Major code refactoring to: 
   - improve code reuse
   - simplify debugging and maintaining code
   - simplify addition of new features
   - improve error handling 
   - some linear mixed model analyses are 50% faster
   - enable additional features for `dreamlet` package that depends heavily on `variancePartition`.  

# variancePartition 1.29.2
 - Rescale weights

# variancePartition 1.28.9
 - March 14, 2023
 - fix rounding error in `makeContrastsDream()`
 - add Pearson residuals to `residuals()`

# variancePartition 1.28.8
 - March 8, 2023
 - add `mvTest()` with features as list

# variancePartition 1.28.7
 - March 7, 2023
 - Fix bug in `makeContrastsDream()` by adding `droplevels()`

# variancePartition 1.28.6
 - March 1, 2023
 - `diffVar()` now fits contrasts estimated in first step

# variancePartition 1.28.5
 - Feb 24, 2023
 - Fix error in `vcov()` when samples are dropped due to covariate having `NA` value

# variancePartition 1.28.4
 - Feb 2, 2023
 - Improve documentation for contrasts

# variancePartition 1.28.3
 - Jan 19, 2023
 - Improve checking and documentation for contrasts

# variancePartition 1.28.2
 - Jan 13, 2023
 - `canCorPairs()` now allows random effects in formula
   - but won't change results

# variancePartition 1.28.1
 - Jan 04, 2023
 - bug fixes

# variancePartition 1.27.17
 - in `mvTest()`, more consistent return values when one features is used

# variancePartition 1.27.16
 - bug fix

# variancePartition 1.27.15
 - fix approximation for `vcov()`
 - add checks for `vcov()` in `test_vcov2()`

# variancePartition 1.27.14
 - fix bug in `topTable()`
 - add `deviance()`
 - update docs
 - update `sqrtMatrix()` to have positive diagonal

# variancePartition 1.27.13
 - add `diffVar()` test of differential variance
 - `dream()` now returns `formula`, `data`, and `hatvalues`
 - define `hatvalues()` for result of `dream()`

# variancePartition 1.27.12
 - in `mvTest()`:
   - default method is now `"FE"`

# variancePartition 1.27.11
 - in `mvTest()`:
   - return number of features
   - return stat.FE and stat.het for RE2C
   - return NA for stat if 1 feature

# variancePartition 1.27.10
 - Add check for very large weights in `voomWithDreamWeights()`
   - follows bug report: https://github.com/GabrielHoffman/variancePartition/issues/66
 - in `mvTest()` change option "LS" to "FE"

# variancePartition 1.27.9
 - bug fix in `mvTest()` for contrasts
   - extensive changes to `vcov()` in order to handle contrasts
 - fix by in subsetting `MArrayLM2` object.  `cov.coefficients.list` are now named based on gene names 

# variancePartition 1.27.8
 - bug fix in `mvTest()`

# variancePartition 1.27.7
 - add `vcov()` for results from `dream()`
 - add `mvTest()` for  multivariate tests on results from `dream()` using `vcov()`

# variancePartition 1.27.6
 - remove FMT functions

# variancePartition 1.27.5
 - update dependencies
 - make `topTable()` generic to work with R 4.2.1 and Bioc 3.16

# variancePartition 1.27.4
 - update dependencies

# variancePartition 1.27.3
 - update filtering of covariates, especially for when many samples are dropped 

# variancePartition 1.27.1
 - update `plotPercentBars` arguments

# variancePartition 1.26.1
 - patch to update docs

# variancePartition 1.25.13
 - Fix compatibility issue with lme4 1.1.29
  - reported https://github.com/GabrielHoffman/variancePartition/issues/51

# variancePartition 1.25.12
 - in `makeContrastsDream()`, fix issue where terms with colon cause and error

# variancePartition 1.25.11
 - fix bug in `dream()` for variables with NA values
 - improve handling of invalid contrasts in `makeContrastsDream()`

# variancePartition 1.25.9
 - for `getContrast()` and `makeContrastsDream()` make sure formula argument is a formula and not a string

# variancePartition 1.25.8
 - small bug fixes

# variancePartition 1.25.7
 - `dream()` now drops samples with missing data gracefully

# variancePartition 1.25.6
 - fix small plotting bug in `plotStratify()` and `plotStratifyBy()`

# variancePartition 1.25.5
 - add `getTreat()` to evaluate `treat()`/`topTreat()` seamlessly on results of `dream()`

# variancePartition 1.25.4
 - in `dream()` set default `ddf = "adaptive"`, which uses "KR" for less than 12 samples 
 - all functions default tp `BPPARAM=SerialParam()`
 - add `eBayes()` to vignette for `dream()`

# variancePartition 1.25.3
 - add genes argument to `plotPercentBars()` 

# variancePartition 1.25.2
 - change `plotPercentBars()` to use generic S4

# variancePartition 1.25.1
 - update handling of `weights` in `voomWithDreamWeights()` and add `applyQualityWeights()`

# variancePartition 1.23.6
 - update `calcVarPart()` with argument `scale=TRUE` allowing the user to disable scaling to fractions

# variancePartition 1.23.5
 - use `RhpcBLASctl::omp_set_num_threads(1)` to use only 1 OpenMP thread for linear algebra within each BiocParallel process
 - update `dream()` so `useWeights=FALSE` works with `lmFit()`

# variancePartition 1.23.4
 - convert some warnings to errors
 - add proper handling of weights to `voomWithDreamWeights()`

# variancePartition 1.23.3
 - add calcVarPart() support for a range of GLM/GLMM's

# variancePartition 1.23.2
 - add flag to checkModelStatus() so warnings are thrown immediately
 - fix export of as.data.frame
 - fixed issues where messages were printed even if quiet=TRUE

# variancePartition 1.23.1
 - New freeze for Bioconductor devel branch

# variancePartition 1.21.10
 - add suppressWarnings flag to makeContrastsDream() 
 - ensure that z.std is finite

# variancePartition 1.21.9
 - fix bug with makeContrastsDream() 

# variancePartition 1.21.8
 - update dream vignette
 - update documentation for makeContrastsDream()
 - fix error in rdf_from_matrices() with eigen() failing

# variancePartition 1.21.7
 - Merge changes to contrast code: https://github.com/GabrielHoffman/variancePartition/pull/32
 - Merge improvments to error checking: https://github.com/GabrielHoffman/variancePartition/pull/28
 - New warning/error if variables in formula have missing data
 - add `makeContrastsDream()`

# variancePartition 1.21.6
 - in dream() allow L to be NULL

# variancePartition 1.21.5
 - Fix roxygen2 documentation to pass R CMD check

# variancePartition 1.21.4
 - variance fractions for fixed effects model is now computed using new method
  - fixes subtle issue with previous version where estimates dependend on of terms in the formula
  - https://github.com/GabrielHoffman/variancePartition/issues/30
 - Update documentation

# variancePartition 1.21.3
 - Faster aggregration after running multiple threads
 - Pulled from https://github.com/GabrielHoffman/variancePartition/pull/27
  - by Ryan C. Thompson
 - eBayes() now works with dream for linear mixed models 
  - add rdf.merMod 

# variancePartition 1.21.2
 - Reduce size of data passed to each thread by only including variables used in the formula.  Applies to multiple functions
 - add more unit tests

# variancePartition 1.21.1
 - fix issue with dream and random slopes

# variancePartition 1.19.20
 - fix bug discovered when the number of features is less than the number of chunks in iterBatch()

# variancePartition 1.19.19
 - simple bug fixes to pass R CMD check

# variancePartition 1.19.18
 - simplify calcVarPart for lm and lmer.  Add compatibility for glm
 - Simplify checkModelStatus.merMod to allow formula (A|B) where A is continuous
 - remove unused "adjust" arguments for clarity

# variancePartition 1.19.17
 - add get_prediction() for results of lm()
 - improve documentation of get_prediction()

# variancePartition 1.19.16
 - in canCorPairs() change statistic used to summarize CCA to Cramer's V.  The difference is very subtle, but is now based on first principles.
 - in dream, check that data is a data.frame
 - dream() defaults to computeResiduals=TRUE for compatability with zenith

# variancePartition 1.19.14 
 - fix issue with residuals() where examples fail

# variancePartition 1.19.13 
 - fix issues with residuals() 
  - https://github.com/GabrielHoffman/variancePartition/issues/18
 - fix issue exporting eBayes, topTable, etc

# variancePartition 1.19.12
 - Improve documentation for contrasts in dream.Rmd
 - check that contrasts sum to zero in plotContrasts.

# variancePartition 1.19.11
 - in voomWithDreamWeights() fix issue with not defining design
  - https://github.com/GabrielHoffman/variancePartition/issues/17

# variancePartition 1.19.10
 - in voomWithDreamWeights() fix issue with returning design matrix 
 - better error if counts can't be converted to matrix
  - https://github.com/GabrielHoffman/variancePartition/issues/15

# variancePartition 1.19.7 
 - Round numbers in plotContrasts()
 - fix issues with strings are passed to formula arguments

# variancePartition 1.19.6 
 - New gives meaning full error message for dream(), etc when variable is not found in data.

# variancePartition 1.19.5 
 - Better error catching when running fitVarPartModel() with fxn that fails
 - add get_prediction() function
  - the following code now can be run in parallel
  fitList = fitVarPartModel( Y, ~ (1|Batch), data, fxn = function(fit){
             B = variancePartition::get_prediction(fit, ~(1|Batch))
             fit@resp$y - B
          }, BPPARAM=SnowParam(3))

# variancePartition 1.19.4
 - Update vignette #3, and update documentation of REML argument 

# variancePartition 1.19.3
 - add new FAQ.Rmd

# variancePartition 1.19.2
 - canCorPairs() now returns NA correlation when two variables have 
  no overlapping observed values
 - plotCorrMatrix() now handles NA correlation values 

# variancePartition 1.19.1
 - Bump to next Bioconductor version

# variancePartition 1.18.3
 - Improve documentation
 - move location of eBayesFMT code

# variancePartition 1.18.2
 - Clean up some code and add documentation
 - document ebayesFMT

# variancePartition 1.18.1

 - Clean up some code and add documentation
 - compute effective degrees of freedom for each model

# variancePartition 1.18.0
 - Bioconductor freeze

# variancePartition 1.17.10
 - fix issue returning residuals from limma
 - resolve issue where dream gives error: r[cbind(1L:p, 1L:p)] <- 1 : subscript out of bounds
  - only occured when no fixed effects were used

# variancePartition 1.17.9
 - Fix issues with compatability with R/4.0.0

# variancePartition 1.17.8
 - Better error message when response contains missing data

# variancePartition 1.17.7
 - Better error message when variable in formula does not exist

# variancePartition 1.17.6
 - comply with new Bioconductor check: _R_CHECK_LENGTH_1_LOGIC2_

# variancePartition 1.17.5
 - if contrasts for dream() is data.frame, convert to matrix

# variancePartition 1.17.4
 - Don't print warnings for residuals() when only one argument passed.
 - fix bug with residuals evaluated with only fixed effects

# variancePartition 1.17.3
 - Allow sparseMatrix for gene expression.  Now saves memory by avoiding conversion to matrix.  Processing sparseMatrix will be slower, but memory usage will be low. 
 - dream(..., computeResiduals=TRUE) now computes residuals and allows use of residuals() function

# variancePartition 1.17.2
 - fix error in voomWithDreamWeights() when design matrix is null

# variancePartition 1.17.1
 - topTable(...,sort.by=) now is correct when and F-test is used
 - fixed issue in classifyTestsF.MArrayLM2, now is much faster

# variancePartition 1.15.8
 - Replace cat() with message()
 - add quiet option to a few functions
 - dream() does not call eBayes() when lmFit is used

# variancePartition 1.15.7
 - Official release to development branch

# variancePartition 1.15.6
 - fix convergence errror when recycling parameters values from first gene
 - add column z.std and F.std to topTable

# variancePartition 1.15.4
 - Add error message when scale of fixed effects causes a problem

# variancePartition 1.15.3
 - Try changing order of eBayes 

# variancePartition 1.15.2
 - Update vignette

# variancePartition 1.15.0
 - Push changes to Bioconductor devel

# variancePartition 1.13.11
 - Fixed bug in voomWithDreamWeights() with approxfun() with ties

# variancePartition 1.13.10
 - fix parallel processing issues

# variancePartition 1.13.9
 - update vignette for new parallel processing backend

# variancePartition 1.13.8
 - apply empirical Bayes when doing F-test
 - hypothesis testing for single coefficients is now included by default, so only need to specify contrast matrix if for more complicated contrasts  
 - add voomWithDreamWeights() for computing observation weights using random effects
 - Add BiocParallel capability with BPPARAM argument
  - allows parallel processing with lower memory usage
 - dream() is now compatable with gene set enrichments from pinnacle (software comming soon)
 
# variancePartition 1.13.7
 - in dream(), add support for genes annotation in DGElist()
 - in dream(), automatically evaluate contrasts for all single coefficients
 - add future compatability for gene set enrichments method "pinnacle"

# variancePartition 1.13.5
 - add plotContrasts() to dream vignette

# variancePartition 1.13.4
 - export classes to fix bug with class "varPartResults" not being defined
 - Thanks Megan Behringer

# variancePartition 1.13.3
 - add plotContrasts()

# variancePartition 1.13.2
 - Enable random slope models in dream, but not for estimating variance fractions
 - Thanks Jonas Zierer

# variancePartition 1.13.1
 - Add progress bar at ETA
=======
# variancePartition 1.12.4
 - Update dream vignette

# variancePartition 1.12.3
    o add plotContrasts()
    o Enable random slope models in dream, but not for estimating variance fractions
    o export classes to fix bug with class "varPartResults" not being defined
     - Thanks Megan Behringer

# variancePartition 1.11.13
 - Fix multithreading issue

# variancePartition 1.11.11
 - dream can handle multiple contrasts at the same time

# variancePartition 1.11.10
 - fix typos in dream vignette
 
# variancePartition 1.11.8
 - Check and stop() if response variable has variance of 0
  - in dream(), fitExtractVarPartModel(), and fitVarPartModel() 
 - add standardized_t_stat() implicitly in eBayes() using MArrayLM2 class
  - this transforms moderated t-statistics to have same degrees of freedom

# variancePartition 1.11.7
 - Simplify object return by dream to be more more similar to lmFit
  - now returns MArrayLM instead of MArrayLMM_lmer
 - if a fixed effects formula is specified (i.e. not random terms)
  - dream call lmFit in the backend
  - getContrast() works seamlessly
 - dream() now returns gene annotation if passed to function

# variancePartition 1.11.6
 - add error checing for L in dream
 - fix typoes in dream vignette
 - fix typoes in theory_practice_random_effects.Rnw

# variancePartition 1.11.5 
 - Add dream function for differential expression for repeated measures with a linear mixed model

# variancePartition 1.11.2 
 - Add warnings to canCorPairs for colinear terms

# variancePartition 1.11.1
 - Add vignette: theory_practice_random_effects.Rnw

# variancePartition 1.9.9
 - Fix issue when package is autoloaded when starting R

# variancePartition 1.9.8
     o Fix issue where if info data.frame contained a column name "gene", 
      fitExtractVarPartModel() would not run

# variancePartition 1.9.6
        o Fix tximport issue with eval=FALSE

# variancePartition 1.9.2
        o Fix vignette

# variancePartition 1.9.1
        o Fix formatting of vignette
        o add description of canCorPairs() function

# variancePartition 1.5.7
 - include splines in foreach .packages

# variancePartition 1.5.6
 - compatibility with tximport v1.3.5


# variancePartition 1.5.5
 - Decrease computing time of effective sample size with ESS() by additional ~10x with sparse solver
 -  fix margins for plotPercentBars()
 - Fix bug for getVarianceComponents() when correlated continous variables are included
 - compatibility with ggplot2 2.2.0
  - center plot titles 
  - fix order of bars in plotPercentBars()
  - legend background to transparent
  - set text to be black
 - include lme4 in foreach .packages
 - change residuals color to not be transparent
 - add CITATION information 
 - plotCorrMatrix now shows dendrogram by default 
 - Estimate run time for fitExtractVarPartModel() / fitVarPartModel()
 - improve warnings for plotPercentBar()
 - improve warnings for plotCorrStructure()
 - define ylab for plotVarPart()
 - add as.matrix.varPartResults() (hidden)
 - define isVaryingCoefficientModel() (hidden)      

# variancePartition 1.3.11
 - in canCorPairs() and other functions, convert formula with as.formula()
 - improve error messages for canCorPairs()


# variancePartition 1.3.10
 - Add plotStratify()
 - Update documentation

# variancePartition 1.3.8
 - Add additional examples to vignette
 - show projected memory usage of fitVarPartModel()

# variancePartition 1.3.7
 - fitVarPartModel warns if names in exprObj and data are not identical
 - residuals() and other functions deal with missing values properly

# variancePartition 1.3.6
 - Small changes to vignette

# variancePartition 1.3.5
 - Fix Bioconductor error

# variancePartition 1.3.4
 - Fix typos

# variancePartition 1.3.3
 - Improve documentation

# variancePartition 1.1.9
 - Update sortCols to handle Measurement.error
 - change backend package structure
 - set Residuals to be grey by default in plotVarPart() and plotPercentBars()
 - add control = lme4::lmerControl(calc.derivs=FALSE, check.rankX="stop.deficient" )
 - add plotCorrStructure

# variancePartition 1.1.8
 - Add ESS.R
 - Add fitVarTest.R
 - use lmerTest by default
 - fix bug checkModelStatus() for variables with backticks in name

# variancePartition 1.1.7
 - GPL License

# variancePartition 1.1.6
 - Move packages from Depends to Imports
 - For clarity, replace = with <- in parts of examples and vignette
 - Stop cluster in examples to solve error on Windows machines

# variancePartition 1.1.5
 - Stop cluster in vignette to solve error on Windows machines

# variancePartition 1.1.4
 - Fix Bioconductor check error

# variancePartition 1.1.3
    o Add details to vignette
 - Fix ggplot2 compatibility issues

# variancePartition 1.1.2
 - Add details to vignette

# variancePartition 1.1.1
 - add plotPercentBars() to vizualize variance fractions for a subset of genes
 - add ESS() to compute effective sample size
 - fix x.labels argument in plotStratifyBy().  Previously, this argument was not used correctly

# variancePartition 1.0.0
 - Release to Bioconductor 3.2

# variancePartition 0.99.9
 - add legend argument to plotStratifyBy()
 - improve warnings / errors for varying coefficient models
 - allow user to manually adjust cutoff for determining when design matrix is singular
   - changed default cutoff to 0.999 from 0.99

# variancePartition 0.99.8
 - improve warnings / errors when design matrix is close to or exactly singular

# variancePartition 0.99.7
 - added new class varPartResults to store results of fitExtractVarPartModel() and extractVarPart()
   - the user will not notice any change, only the backend is different
        o Allow computation of adjusted ICC in addition to ICC.
 - add warning when categorical variables are modeled as fixed effects
 - fix computation of variance fractions for varying coefficient models
 - add getVarianceComponents() to return variances from lmer() or lm() model fit
 - showWarnings=FALSE suppresses warning messages
 - add fxn argument to fitVarPartModel to evaluate any function on the model fit

# variancePartition 0.99.6
 - Update DESCRIPTION information

# variancePartition 0.99.5
 - residuals deals with missing data gracefully and returns a matrix

# variancePartition 0.99.4
       o add documentation for example datasets
       o convert calcVarPart() to S4 from S3 function call
       o fix typos in vignette

# variancePartition 0.99.3
       o fitVarPartModel() and fitExractVarPartModel() use S4 instead of S3 calls 

# variancePartition 0.99.2
 - rename sort.varParFrac to sortCols
 - support ExpressionSet
 - change options for plotStratifyBy()

# variancePartition 0.0.12
 - add plotStratifyBy()
 - update documentation

# variancePartition 0.0.11
 - fix sort to work on correct argument

# variancePartition 0.0.10
 - fitExtractVarPartModel() and fitVarPartModel() now take subset argument
 - throw warning when no Intercept is specified
 - if using lmer, warning if categorical variable is modeled as fixed effect
 - fixed calcVarPart bug with reporting too few variances for multicategory fixed effects
 - add colinearityScore

# variancePartition 0.0.9
    o function now use the precision weights when specified
 - remove warning about unspecified weights, when useWeights=TRUE
 - fix issue with sort with only one variable
 - add main argument to plotVarPart

# variancePartition 0.0.8
        o Bug update to vignette and simulated data

# variancePartition 0.0.7
 - remove 'variable' from xlab of plotVarPart

# variancePartition 0.0.6
 - set REML=FALSE to default.  This fixes issues of inaccurate variance estiamtes, and makes lmer() results more concordant with lm() results
 - Fix residuals function when lm or lmer is used
 - fix useWeights argument error for fitExtractVarPartModel()

# variancePartition 0.0.5
        o plotVarPart() as ylim=c(0,100) as default, and can be changed by user
 - fitExtractVarPartModel() labels rows correctly


# variancePartition 0.0.4
 - Add sort() for output of extractVarPart() or fitExtractVarPartModel()

# variancePartition 0.0.3
 - add residuals function
 - add fitExtractVarPartModel()
 - foreach loops us iterators

# variancePartition 0.99.0
 - Initial version
