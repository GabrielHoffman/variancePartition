# Additional visualizations of variance structure

The correlation structure between samples in complex study designs can
be decomposed into the contribution of multiple dimensions of variation.
`variancePartition` provides a statistical and visualization framework
to interpret sources of variation. Here I describe a visualization of
the correlation structure between samples for a single gene.

In the example dataset described in the main vignette, samples are
correlated because they can come from the same individual or the same
tissue. The function
[`plotCorrStructure()`](http://DiseaseNeurogenomics.github.io/variancePartition/reference/plotCorrStructure.md)
shows the correlation structure caused by each variable as well and the
joint correlation structure. Figure 1 shows the correlation between
samples from the same individual where (a) shows the samples sorted
based on clustering of the correlation matrix and (b) shows the original
order. Figure 1 c) and d) shows the same type of plot except
demonstrating the effect of tissue. The total correlation structure from
summing individual and tissue correlation matricies is shown in Figure
2. The code to generate these plots is shown below.

## Plot variance structure

``` r

# Fit linear mixed model and examine correlation stucture
# for one gene
data(varPartData)

form <- ~ Age + (1 | Individual) + (1 | Tissue)

fitList <- fitVarPartModel(geneExpr[1:2, ], form, info)

# focus on one gene
fit <- fitList[[1]]
```

### By Individual

#### Reorder samples

``` r

# Figure 1a
# correlation structure based on similarity within Individual
# reorder samples based on clustering
plotCorrStructure(fit, "Individual")
```

![](additional_visualization_files/figure-html/corStructa-1.png)

#### Original order of samples

``` r

# Figure 1b
# use original order of samples
plotCorrStructure(fit, "Individual", reorder = FALSE)
```

![](additional_visualization_files/figure-html/corStructb-1.png)

### By Tissue

#### Reorder samples

``` r

# Figure 1c
# correlation structure based on similarity within Tissue
# reorder samples based on clustering
plotCorrStructure(fit, "Tissue")
```

![](additional_visualization_files/figure-html/corStructc-1.png)

#### Original order of samples

``` r

# Figure 1d
# use original order of samples
plotCorrStructure(fit, "Tissue", reorder = FALSE)
```

![](additional_visualization_files/figure-html/corStructd-1.png)

### By Individual and Tissue

#### Reorder samples

``` r

# Figure 2a
# correlation structure based on similarity within
# Individual *and* Tissue, reorder samples based on clustering
plotCorrStructure(fit)
```

![](additional_visualization_files/figure-html/corStructe-1.png)

#### Original order of samples

``` r

# Figure 2b
# use original order of samples
plotCorrStructure(fit, reorder = FALSE)
```

![](additional_visualization_files/figure-html/corStructf-1.png)

## Session Info

    ## R version 4.5.1 (2025-06-13)
    ## Platform: aarch64-apple-darwin23.6.0
    ## Running under: macOS Sonoma 14.7.1
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /opt/homebrew/Cellar/openblas/0.3.33/lib/libopenblasp-r0.3.33.dylib;  LAPACK version 3.12.0
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
    ## other attached packages:
    ## [1] variancePartition_2.0.0 BiocParallel_1.44.0     limma_3.66.0           
    ## [4] ggplot2_4.0.3           knitr_1.51             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] Rdpack_2.6.6                bitops_1.0-9               
    ##   [3] rlang_1.3.0                 magrittr_2.0.5             
    ##   [5] otel_0.2.0                  matrixStats_1.5.0          
    ##   [7] compiler_4.5.1              reshape2_1.4.5             
    ##   [9] systemfonts_1.3.2           vctrs_0.7.3                
    ##  [11] stringr_1.6.0               pkgconfig_2.0.3            
    ##  [13] fastmap_1.2.0               backports_1.5.1            
    ##  [15] XVector_0.50.0              caTools_1.18.3             
    ##  [17] rmarkdown_2.31              nloptr_2.2.1               
    ##  [19] ragg_1.5.2                  purrr_1.2.2                
    ##  [21] xfun_0.60                   cachem_1.1.0               
    ##  [23] jsonlite_2.0.0              EnvStats_3.1.0             
    ##  [25] remaCor_0.0.20              DelayedArray_0.36.1        
    ##  [27] broom_1.0.13                parallel_4.5.1             
    ##  [29] R6_2.6.1                    stringi_1.8.7              
    ##  [31] bslib_0.11.0                RColorBrewer_1.1-3         
    ##  [33] parallelly_1.48.0           car_3.1-5                  
    ##  [35] boot_1.3-32                 GenomicRanges_1.62.1       
    ##  [37] jquerylib_0.1.4             numDeriv_2016.8-1.1        
    ##  [39] Rcpp_1.1.2                  Seqinfo_1.0.0              
    ##  [41] SummarizedExperiment_1.40.0 iterators_1.0.14           
    ##  [43] IRanges_2.44.0              Matrix_1.7-5               
    ##  [45] splines_4.5.1               tidyselect_1.2.1           
    ##  [47] dichromat_2.0-0.1           abind_1.4-8                
    ##  [49] yaml_2.3.12                 gplots_3.3.0               
    ##  [51] codetools_0.2-20            plyr_1.8.9                 
    ##  [53] lattice_0.22-9              tibble_3.3.1               
    ##  [55] lmerTest_3.2-1              Biobase_2.70.0             
    ##  [57] withr_3.0.3                 S7_0.2.2                   
    ##  [59] evaluate_1.0.5              desc_1.4.3                 
    ##  [61] pillar_1.11.1               MatrixGenerics_1.22.0      
    ##  [63] carData_3.0-6               KernSmooth_2.23-26         
    ##  [65] stats4_4.5.1                reformulas_0.4.4           
    ##  [67] generics_0.1.4              fastglmm_0.4.11            
    ##  [69] S4Vectors_0.48.1            scales_1.4.0               
    ##  [71] aod_1.3.3                   minqa_1.2.8                
    ##  [73] gtools_3.9.5                RhpcBLASctl_0.23-42        
    ##  [75] glue_1.8.1                  tools_4.5.1                
    ##  [77] fANCOVA_0.6-1               lme4_2.1-0                 
    ##  [79] locfit_1.5-9.12             fs_2.1.0                   
    ##  [81] mvtnorm_1.4-2               grid_4.5.1                 
    ##  [83] tidyr_1.3.2                 rbibutils_2.4.1            
    ##  [85] edgeR_4.8.2                 nlme_3.1-170               
    ##  [87] Formula_1.2-5               cli_3.6.6                  
    ##  [89] textshaping_1.0.5           S4Arrays_1.10.1            
    ##  [91] dplyr_1.2.1                 corpcor_1.6.10             
    ##  [93] gtable_0.3.6                DESeq2_1.50.2              
    ##  [95] sass_0.4.10                 digest_0.6.39              
    ##  [97] BiocGenerics_0.56.0         SparseArray_1.10.10        
    ##  [99] pbkrtest_0.5.5              htmlwidgets_1.6.4          
    ## [101] farver_2.1.2                htmltools_0.5.9            
    ## [103] pkgdown_2.2.1               lifecycle_1.0.5            
    ## [105] statmod_1.5.2               MASS_7.3-66
