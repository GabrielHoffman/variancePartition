# Package index

## Core functions

- [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  : Differential expression with linear mixed model

- [`eBayes()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/eBayes-method.md)
  : eBayes generic for for MArrayLM and MArrayLM2

- [`fitExtractVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitExtractVarPartModel-method.md)
  : Fit linear (mixed) model, report variance fractions

- [`makeContrastsDream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/makeContrastsDream.md)
  : Construct Matrix of Custom Contrasts

- [`mvTest()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/mvTest-method.md)
  :

  Multivariate tests on results from
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

- [`topTable()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/toptable-method.md)
  : Table of Top Genes from Linear Model Fit

- [`varpart()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varpart.md)
  : Variance Partitioning Analysis

- [`voomWithDreamWeights()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/voomWithDreamWeights.md)
  :

  Transform RNA-Seq Data Ready for Linear Mixed Modelling with
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)

## Extract results from model fits

- [`BIC(`*`<MArrayLM>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/BIC.md)
  [`BIC(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/BIC.md)
  : BIC from model fit

- [`ESS()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/ESS-method.md)
  : Effective sample size

- [`classifyTestsF(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/classifyTestsF-MArrayLM2-method.md)
  : Multiple Testing Genewise Across Contrasts

- [`classifyTestsF()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/classifyTestsF.md)
  : Multiple Testing Genewise Across Contrasts

- [`hatvalues(`*`<MArrayLM>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/hatvalues-method.md)
  [`hatvalues(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/hatvalues-method.md)
  : Compute hatvalues

- [`logLik(`*`<MArrayLM>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/logLik.md)
  [`logLik(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/logLik.md)
  : Log-likelihood from model fit

- [`rdf()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/rdf.md)
  : Residual degrees of freedom

- [`rdf.merMod()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/rdf.merMod.md)
  : Approximate residual degrees of freedom

- [`rdf_from_matrices()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/rdf_from_matrices.md)
  : Fast approximate residual degrees of freedom

- [`residuals.MArrayLM2()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/residuals.md)
  [`residuals(`*`<MArrayLM>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/residuals.md)
  [`residuals(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/residuals.md)
  : Residuals for result of dream

- [`residuals(`*`<VarParFitList>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/residuals-VarParFitList-method.md)
  : Residuals from model fit

- [`vcov(`*`<MArrayLM>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/vcov-MArrayLM-method.md)
  :

  Co-variance matrix for
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  fit

- [`vcov(`*`<MArrayLM2>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/vcov-MArrayLM2-method.md)
  :

  Co-variance matrix for
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  fit

- [`vcovSqrt()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/vcovSqrt-method.md)
  :

  Sqrt of co-variance matrix for
  [`dream()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dream-method.md)
  fit

## Plots

- [`plotCompareP()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotCompareP-method.md)
  : Compare p-values from two analyses
- [`plotContrasts()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotContrasts.md)
  : Plot representation of contrast matrix
- [`plotCorrMatrix()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotCorrMatrix.md)
  : plotCorrMatrix
- [`plotCorrStructure()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotCorrStructure.md)
  : plotCorrStructure
- [`plotPercentBars()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotPercentBars-method.md)
  : Bar plot of gene fractions
- [`plotStratify()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotStratify.md)
  : plotStratify
- [`plotStratifyBy()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotStratifyBy.md)
  : plotStratifyBy
- [`plotTrendVP()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotTrendVP-methods.md)
  : Plot Trend of Variance Fractions Vs Count Magnitude
- [`plotVarPart()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotVarPart-method.md)
  : Violin plot of variance fractions
- [`plotVarianceEstimates()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotVarianceEstimates.md)
  : Plot Variance Estimates

## Other functions

- [`applyQualityWeights()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/applyQualityWeights.md)
  : Apply pre-specified sample weights
- [`as.data.frame(`*`<varPartResults>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/as.data.frame.varPartResults.md)
  : Convert to data.frame
- [`as.matrix(`*`<varPartResults>`*`)`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/as.matrix-varPartResults-method.md)
  : Convert to matrix
- [`augmentPriorCount()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/augmentPriorCount.md)
  : Augment observed read counts with prior count
- [`calcVarPart()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/calcVarPart.md)
  : Compute variance statistics
- [`canCorPairs()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/canCorPairs.md)
  : canCorPairs
- [`colinearityScore()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/colinearityScore.md)
  : Collinearity score
- [`deviation()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/deviation-method.md)
  : Deviation from expectation for each observation
- [`diffVar()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/diffVar-method.md)
  : Test differential variance
- [`dscchisq()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/dscchisq.md)
  : Scaled chi-square
- [`extractVarPart()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/extractVarPart.md)
  : Extract variance statistics
- [`fitVarPartModel()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/fitVarPartModel-method.md)
  : Fit linear (mixed) model
- [`getContrast()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/getContrast-method.md)
  : Extract contrast matrix for linear mixed model
- [`getTreat()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/getTreat-method.md)
  : Test if coefficient is different from a specified value
- [`get_prediction()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/get_prediction-method.md)
  : Compute predicted value of formula for linear (mixed) model
- [`ggColorHue()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/ggColorHue.md)
  : Default colors for ggplot
- [`isRunableFormula()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/isRunableFormula.md)
  : Test if formula is full rank on this dataset
- [`reOnly()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/reOnly.md)
  : Adapted from lme4:::reOnly
- [`shrinkageMetric()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/shrinkageMetric.md)
  : Shrinkage metric for eBayes
- [`sortCols()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/sortCols-method.md)
  : Sort variance partition statistics
- [`varPartConfInf()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varPartConfInf.md)
  : Linear mixed model confidence intervals
- [`varPartData`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varPartDEdata.md)
  : A simulated dataset of gene counts
- [`varPartData`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varPartData.md)
  : Simulation dataset for examples

## Classes

- [`mvTest_input-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/mvTest_input-class.md)
  : Class mvTest_input
- [`varParFrac-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varParFrac-class.md)
  : Class varParFrac
- [`varPartResults-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/varPartResults-class.md)
  : Class varPartResults
- [`MArrayLM2-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/MArrayLM2-class.md)
  : Class MArrayLM2
- [`VarParCIList-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/VarParCIList-class.md)
  : Class VarParCIList
- [`VarParFitList-class`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/VarParFitList-class.md)
  : Class VarParFitList
