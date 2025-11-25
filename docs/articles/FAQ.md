# Frequently asked questions

## Interperting the residual variance

In general, I recommend against interpreting the fraction of variance
explained by residuals. This fraction is driven by:

1.  the particulars of the study design
2.  measurement precision (i.e. high read counts give more precise
    measurements)
3.  biological variability
4.  technical variability (i.e. batch effects).

If you have additional variables that explain variation in measured gene
expression, you should include them in order to avoid confounding with
your variable of interest. But a particular residual fraction is not
‘good’ or ‘bad’ and is not a good metric of determining whether more
variables should be included.

## Current GitHub issues

See [GitHub
page](https://github.com/DiseaseNeuroGenomics/variancePartition/issues)
for up-to-date responses to users’ questions.

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
