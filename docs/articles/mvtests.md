# Multivariate tests

Results from the univariate regressions performed using can be combined
in a post-processing step to perform multivariate hypothesis testing. In
this example, we fit on transcript-level counts and then perform
multivariate hypothesis testing by combining transcripts at the
gene-level. This is done with the function.

## Import transcript-level counts

Read in transcript counts from the package.

``` r

library(readr)
library(tximport)
library(tximportData)

# specify directory
path <- system.file("extdata", package = "tximportData")

# read sample meta-data
samples <- read.table(file.path(path, "samples.txt"), header = TRUE)
samples.ext <- read.table(file.path(path, "samples_extended.txt"), header = TRUE, sep = "\t")

# read assignment of transcripts to genes
# remove genes on the PAR, since these are present twice
tx2gene <- read_csv(file.path(path, "tx2gene.gencode.v27.csv"))
tx2gene <- tx2gene[grep("PAR_Y", tx2gene$GENEID, invert = TRUE), ]

# read transcript-level quatifictions
files <- file.path(path, "salmon", samples$run, "quant.sf.gz")
txi <- tximport(files, type = "salmon", txOut = TRUE)

# Create metadata simulating two conditions
sampleTable <- data.frame(condition = factor(rep(c("A", "B"), each = 3)))
rownames(sampleTable) <- paste0("Sample", 1:6)
```

## Standard dream analysis

Perform standard analysis at the transcript-level

``` r

library(variancePartition)
library(edgeR)

# Prepare transcript-level reads
dge <- DGEList(txi$counts)
design <- model.matrix(~condition, data = sampleTable)
isexpr <- filterByExpr(dge, design)
dge <- dge[isexpr, ]
dge <- calcNormFactors(dge)

# Estimate precision weights
vobj <- voomWithDreamWeights(dge, ~condition, sampleTable)

# Fit regression model one transcript at a time
fit <- dream(vobj, ~condition, sampleTable)
fit <- eBayes(fit)
```

## Multivariate analysis

Combine the transcript-level results at the gene-level. The mapping
between transcript and gene is stored in as a list.

``` r

# Prepare transcript to gene mapping
# keep only transcripts present in vobj
# then convert to list with key GENEID and values TXNAMEs
keep <- tx2gene$TXNAME %in% rownames(vobj)
tx2gene.lst <- unstack(tx2gene[keep, ])

# Run multivariate test on entries in each feature set
# Default method is "FE.empirical", but use "FE" here to reduce runtime
res <- mvTest(fit, vobj, tx2gene.lst, coef = "conditionB", method = "FE")

# truncate gene names since they have version numbers
# ENST00000498289.5 -> ENST00000498289
res$ID.short <- gsub("\\..+", "", res$ID)
```

## Gene set analysis

Perform gene set analysis using on the gene-level test statistics.

``` r

# must have zenith > v1.0.2
library(zenith)
library(GSEABase)

gs <- get_MSigDB("C1", to = "ENSEMBL")

df_gsa <- zenithPR_gsa(res$stat, res$ID.short, gs, inter.gene.cor = .05)

head(df_gsa)
```

    ##          NGenes Correlation     delta       se      p.less    p.greater       PValue Direction
    ## chr7p13      28        0.05  7.144240 2.034359 0.999776828 0.0002231723 0.0004463445        Up
    ## chr11p13     32        0.05 -5.752953 1.982804 0.001859931 0.9981400686 0.0037198628      Down
    ## chr4p14      25        0.05 -5.077180 2.084132 0.007428521 0.9925714788 0.0148570424      Down
    ## chr2q37      75        0.05  3.571510 1.758652 0.978855126 0.0211448742 0.0422897483        Up
    ## chr2q36      21        0.05 -4.130558 2.168488 0.028411437 0.9715885626 0.0568228749      Down
    ## chr18q22     18        0.05  4.108195 2.252712 0.965889229 0.0341107714 0.0682215427        Up
    ##                FDR  Geneset     coef
    ## chr7p13  0.1129252  chr7p13 zenithPR
    ## chr11p13 0.4705626 chr11p13 zenithPR
    ## chr4p14  0.9996791  chr4p14 zenithPR
    ## chr2q37  0.9996791  chr2q37 zenithPR
    ## chr2q36  0.9996791  chr2q36 zenithPR
    ## chr18q22 0.9996791 chr18q22 zenithPR

## Session info

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
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] org.Hs.eg.db_3.21.0      msigdbr_25.1.1           GSEABase_1.70.1         
    ##  [4] graph_1.86.0             annotate_1.86.1          XML_3.99-0.20           
    ##  [7] AnnotationDbi_1.70.0     IRanges_2.42.0           S4Vectors_0.48.0        
    ## [10] Biobase_2.68.0           BiocGenerics_0.54.1      generics_0.1.4          
    ## [13] zenith_1.10.0            edgeR_4.6.3              variancePartition_1.39.5
    ## [16] BiocParallel_1.42.2      limma_3.64.3             ggplot2_4.0.1           
    ## [19] tximportData_1.36.0      tximport_1.36.1          readr_2.1.6             
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3          jsonlite_2.0.0              magrittr_2.0.4             
    ##   [4] farver_2.1.2                nloptr_2.2.1                rmarkdown_2.30             
    ##   [7] fs_1.6.6                    ragg_1.5.0                  vctrs_0.6.5                
    ##  [10] memoise_2.0.1               minqa_1.2.8                 RCurl_1.98-1.17            
    ##  [13] progress_1.2.3              S4Arrays_1.8.1              htmltools_0.5.8.1          
    ##  [16] curl_7.0.0                  broom_1.0.10                SparseArray_1.8.1          
    ##  [19] sass_0.4.10                 KernSmooth_2.23-26          bslib_0.9.0                
    ##  [22] htmlwidgets_1.6.4           desc_1.4.3                  pbkrtest_0.5.5             
    ##  [25] plyr_1.8.9                  cachem_1.1.0                lifecycle_1.0.4            
    ##  [28] iterators_1.0.14            pkgconfig_2.0.3             Matrix_1.7-4               
    ##  [31] R6_2.6.1                    fastmap_1.2.0               GenomeInfoDbData_1.2.14    
    ##  [34] rbibutils_2.4               MatrixGenerics_1.20.0       digest_0.6.39              
    ##  [37] numDeriv_2016.8-1.1         textshaping_1.0.4           GenomicRanges_1.60.0       
    ##  [40] RSQLite_2.4.4               abind_1.4-8                 httr_1.4.7                 
    ##  [43] compiler_4.5.1              bit64_4.6.0-1               aod_1.3.3                  
    ##  [46] withr_3.0.2                 S7_0.2.1                    backports_1.5.0            
    ##  [49] DBI_1.2.3                   gplots_3.2.0                MASS_7.3-65                
    ##  [52] DelayedArray_0.34.1         corpcor_1.6.10              gtools_3.9.5               
    ##  [55] caTools_1.18.3              tools_4.5.1                 remaCor_0.0.20             
    ##  [58] glue_1.8.0                  nlme_3.1-168                grid_4.5.1                 
    ##  [61] reshape2_1.4.5              gtable_0.3.6                tzdb_0.5.0                 
    ##  [64] tidyr_1.3.1                 hms_1.1.4                   XVector_0.48.0             
    ##  [67] pillar_1.11.1               stringr_1.6.0               vroom_1.6.6                
    ##  [70] babelgene_22.9              splines_4.5.1               dplyr_1.1.4                
    ##  [73] lattice_0.22-7              bit_4.6.0                   tidyselect_1.2.1           
    ##  [76] locfit_1.5-9.12             Biostrings_2.76.0           knitr_1.50                 
    ##  [79] reformulas_0.4.2            SummarizedExperiment_1.38.1 RhpcBLASctl_0.23-42        
    ##  [82] xfun_0.54                   statmod_1.5.1               matrixStats_1.5.0          
    ##  [85] KEGGgraph_1.68.0            stringi_1.8.7               UCSC.utils_1.4.0           
    ##  [88] yaml_2.3.10                 boot_1.3-32                 evaluate_1.0.5             
    ##  [91] codetools_0.2-20            tibble_3.3.0                Rgraphviz_2.52.0           
    ##  [94] cli_3.6.5                   RcppParallel_5.1.11-1       xtable_1.8-4               
    ##  [97] systemfonts_1.3.1           Rdpack_2.6.4                jquerylib_0.1.4            
    ## [100] dichromat_2.0-0.1           Rcpp_1.1.0                  GenomeInfoDb_1.44.3        
    ## [103] zigg_0.0.2                  EnvStats_3.1.0              png_0.1-8                  
    ## [106] Rfast_2.1.5.2               parallel_4.5.1              assertthat_0.2.1           
    ## [109] pkgdown_2.2.0               blob_1.2.4                  prettyunits_1.2.0          
    ## [112] bitops_1.0-9                lme4_1.1-38                 mvtnorm_1.3-3              
    ## [115] lmerTest_3.1-3              scales_1.4.0                purrr_1.2.0                
    ## [118] crayon_1.5.3                fANCOVA_0.6-1               rlang_1.1.6                
    ## [121] EnrichmentBrowser_2.38.1    KEGGREST_1.48.1

\<\>

## References
