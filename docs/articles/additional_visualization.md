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
    ## other attached packages:
    ## [1] variancePartition_1.39.5 BiocParallel_1.42.2      limma_3.64.3            
    ## [4] ggplot2_4.0.1            knitr_1.50              
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_1.2.1    dplyr_1.1.4         farver_2.1.2       
    ##  [4] S7_0.2.1            bitops_1.0-9        fastmap_1.2.0      
    ##  [7] digest_0.6.39       lifecycle_1.0.4     statmod_1.5.1      
    ## [10] magrittr_2.0.4      compiler_4.5.1      rlang_1.1.6        
    ## [13] sass_0.4.10         tools_4.5.1         yaml_2.3.10        
    ## [16] htmlwidgets_1.6.4   plyr_1.8.9          RColorBrewer_1.1-3 
    ## [19] KernSmooth_2.23-26  withr_3.0.2         purrr_1.2.0        
    ## [22] numDeriv_2016.8-1.1 BiocGenerics_0.54.1 desc_1.4.3         
    ## [25] grid_4.5.1          aod_1.3.3           caTools_1.18.3     
    ## [28] scales_1.4.0        gtools_3.9.5        iterators_1.0.14   
    ## [31] MASS_7.3-65         dichromat_2.0-0.1   cli_3.6.5          
    ## [34] mvtnorm_1.3-3       rmarkdown_2.30      ragg_1.5.0         
    ## [37] reformulas_0.4.2    generics_0.1.4      reshape2_1.4.5     
    ## [40] minqa_1.2.8         cachem_1.1.0        stringr_1.6.0      
    ## [43] splines_4.5.1       parallel_4.5.1      matrixStats_1.5.0  
    ## [46] vctrs_0.6.5         boot_1.3-32         Matrix_1.7-4       
    ## [49] jsonlite_2.0.0      pbkrtest_0.5.5      systemfonts_1.3.1  
    ## [52] jquerylib_0.1.4     tidyr_1.3.1         glue_1.8.0         
    ## [55] nloptr_2.2.1        pkgdown_2.2.0       codetools_0.2-20   
    ## [58] stringi_1.8.7       gtable_0.3.6        EnvStats_3.1.0     
    ## [61] lme4_1.1-38         lmerTest_3.1-3      tibble_3.3.0       
    ## [64] remaCor_0.0.20      pillar_1.11.1       htmltools_0.5.8.1  
    ## [67] gplots_3.2.0        R6_2.6.1            textshaping_1.0.4  
    ## [70] Rdpack_2.6.4        evaluate_1.0.5      lattice_0.22-7     
    ## [73] Biobase_2.68.0      rbibutils_2.4       backports_1.5.0    
    ## [76] RhpcBLASctl_0.23-42 broom_1.0.10        fANCOVA_0.6-1      
    ## [79] corpcor_1.6.10      bslib_0.9.0         Rcpp_1.1.0         
    ## [82] nlme_3.1-168        xfun_0.54           fs_1.6.6           
    ## [85] pkgconfig_2.0.3
