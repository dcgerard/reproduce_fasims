
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2019)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2019).

If you find a bug, please create an
[issue](https://github.com/dcgerard/reproduce_fasims/issues).

## Instructions:

1.  Download and Install PEER: <https://github.com/PMBio/peer/wiki>. I
    had to build the R package from the cloned repository, as installing
    it from the source package didn’t work for me.

2.  Download and install powsimR: <https://github.com/bvieth/powsimR>.

3.  Install the appropriate R packages
    
    ``` r
    install.packages(c("seqgendiff",
                       "tidyverse", 
                       "BiocManager", 
                       "vroom", 
                       "fastICA", 
                       "ssvd",
                       "sparsepca",
                       "elasticnet",
                       "doSNOW", 
                       "ggthemes",
                       "latex2exp"))
    BiocManager::install(c("SummarizedExperiment", 
                           "pcaMethods", 
                           "sva",
                           "edgeR",
                           "DESeq2", 
                           "limma"))
    ```

4.  Download data from the GTEx portal: <https://gtexportal.org/home/>.
    You’ll need to place the following files in the “data/gtex” folder:
    
    1.  GTEx\_Analysis\_2016-01-15\_v7\_RNASeQCv1.1.8\_gene\_reads.gct
    2.  GTEx\_v7\_Annotations\_SubjectPhenotypesDS.txt
    3.  GTEx\_v7\_Annotations\_SampleAttributesDS.txt
    
    I cannot release these data, but it’s free to sign up for access to
    GTEx.

5.  Adjust the Makefile. Change the `nc` and `rexec` variables in the
    Makefile according to your local computing environment. For example,
    you would need to decrease `nc` if you have fewer than 12 CPU cores.
    You can check the number of CPU cores you have by typing the
    following in R:
    
    ``` r
    parallel::detectCores()
    ```

6.  Run `make`. To reproduce all of the results in the paper, just type
    in the terminal
    
    ``` bash
    make
    ```
    
    To just reproduce the results comparing the `powsimR` datasets from
    the `seqgendiff` datasets, run in the terminal:
    
    ``` bash
    make powsimr
    ```
    
    To just reproduce the differential expression simulations, run in
    the terminal:
    
    ``` bash
    make diff_exp
    ```
    
    To just reproduce the factor analysis simulations, run in the
    terminal:
    
    ``` bash
    make FAsims
    ```
    
    To just reproduce the plot of the flexible class of mixtures of
    binomials and negative binomials, run in the terminal:
    
    ``` bash
    make NBplots
    ```
    
    To just reproduce the simulations evaluating the Monte Carlo
    correlation estimator, run in the terminal:
    
    ``` bash
    make corest
    ```

7.  Get coffee. Though most of the results will be generated rather
    quickly, the `make FAsims` call will take a very long time (mostly
    because of PEER). You should get some coffee\! Here is a list of
    some of my favorite places:
    
      - Washington, DC
          - [Colony
            Club](https://www.yelp.com/biz/colony-club-washington)
          - [Grace Street
            Coffee](https://www.yelp.com/biz/grace-street-coffee-georgetown)
          - [Shop Made in
            DC](https://www.yelp.com/biz/shop-made-in-dc-washington)
      - Chicago
          - [Sawada
            Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
          - [Plein Air
            Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
      - Seattle
          - [Bauhaus
            Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
          - [Cafe
            Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
      - Columbus
          - [Yeah, Me
            Too](https://www.yelp.com/biz/yeah-me-too-columbus)
          - [Stauf’s Coffee
            Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)

## Package Versions

If you are having trouble reproducing these results, check your package
versions. These are the ones that I used:

``` r
sessionInfo()
#> R version 3.6.0 (2019-04-26)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.2 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas/libblas.so.3
#> LAPACK: /usr/lib/x86_64-linux-gnu/libopenblasp-r0.2.20.so
#> 
#> locale:
#>  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#>  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#>  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#>  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#>  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
#> [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
#> 
#> attached base packages:
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] powsimR_1.1.4               gamlss.dist_5.1-4          
#>  [3] MASS_7.3-51.4               peer_1.0                   
#>  [5] DESeq2_1.24.0               edgeR_3.26.4               
#>  [7] limma_3.40.2                sva_3.32.1                 
#>  [9] genefilter_1.66.0           mgcv_1.8-28                
#> [11] nlme_3.1-140                pcaMethods_1.76.0          
#> [13] SummarizedExperiment_1.14.0 DelayedArray_0.10.0        
#> [15] BiocParallel_1.18.0         matrixStats_0.54.0         
#> [17] Biobase_2.44.0              GenomicRanges_1.36.0       
#> [19] GenomeInfoDb_1.20.0         IRanges_2.18.1             
#> [21] S4Vectors_0.22.0            BiocGenerics_0.30.0        
#> [23] latex2exp_0.4.0             ggthemes_4.2.0             
#> [25] doSNOW_1.0.16               snow_0.4-3                 
#> [27] iterators_1.0.10            foreach_1.4.4              
#> [29] elasticnet_1.1.1            lars_1.2                   
#> [31] sparsepca_0.1.2             ssvd_1.0                   
#> [33] fastICA_1.2-1               vroom_1.0.1                
#> [35] BiocManager_1.30.4          forcats_0.4.0              
#> [37] stringr_1.4.0               dplyr_0.8.1                
#> [39] purrr_0.3.2                 readr_1.3.1                
#> [41] tidyr_0.8.3                 tibble_2.1.3               
#> [43] ggplot2_3.1.1               tidyverse_1.2.1            
#> [45] seqgendiff_0.3.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0             softImpute_1.4            
#>   [3] minpack.lm_1.2-1           pbapply_1.4-0             
#>   [5] lattice_0.20-38            haven_2.1.0               
#>   [7] penalized_0.9-51           blob_1.1.1                
#>   [9] survival_2.44-1.1          later_0.8.0               
#>  [11] nloptr_1.2.1               DBI_1.0.0                 
#>  [13] R.utils_2.8.0              SingleCellExperiment_1.6.0
#>  [15] Linnorm_2.8.0              dqrng_0.2.1               
#>  [17] zlibbioc_1.30.0            MatrixModels_0.4-1        
#>  [19] pspline_1.0-18             SDMTools_1.1-221.1        
#>  [21] htmlwidgets_1.3            mvtnorm_1.0-10            
#>  [23] future_1.13.0              UpSetR_1.4.0              
#>  [25] scater_1.12.2              irlba_2.3.3               
#>  [27] DEoptimR_1.0-8             Rcpp_1.0.1                
#>  [29] KernSmooth_2.23-15         DT_0.6                    
#>  [31] promises_1.0.1             gdata_2.18.0              
#>  [33] DDRTree_0.1.5              vegan_2.5-5               
#>  [35] Hmisc_4.2-0                ShortRead_1.42.0          
#>  [37] apcluster_1.4.7            RSpectra_0.15-0           
#>  [39] msir_1.3.2                 mnormt_1.5-5              
#>  [41] digest_0.6.19              png_0.1-7                 
#>  [43] qlcMatrix_0.9.7            sctransform_0.2.0         
#>  [45] cowplot_0.9.4              glmnet_2.0-18             
#>  [47] pkgconfig_2.0.2            docopt_0.6.1              
#>  [49] DelayedMatrixStats_1.6.0   ggbeeswarm_0.6.0          
#>  [51] minqa_1.2.4                lavaan_0.6-3              
#>  [53] reticulate_1.12            spam_2.2-2                
#>  [55] beeswarm_0.2.3             modeltools_0.2-22         
#>  [57] xfun_0.7                   zoo_1.8-6                 
#>  [59] tidyselect_0.2.5           ZIM_1.1.0                 
#>  [61] reshape2_1.4.3             kernlab_0.9-27            
#>  [63] ica_1.0-2                  pcaPP_1.9-73              
#>  [65] EDASeq_2.18.0              viridisLite_0.3.0         
#>  [67] rtracklayer_1.44.0         rlang_0.3.4               
#>  [69] hexbin_1.27.3              glue_1.3.1                
#>  [71] metap_1.1                  RColorBrewer_1.1-2        
#>  [73] registry_0.5-1             modelr_0.1.4              
#>  [75] fpc_2.2-2                  pkgmaker_0.27             
#>  [77] fields_9.8-3               SparseM_1.77              
#>  [79] gbRd_0.4-11                httpuv_1.5.1              
#>  [81] class_7.3-15               BPSC_0.99.2               
#>  [83] RMTstat_0.3                BiocNeighbors_1.2.0       
#>  [85] annotate_1.62.0            jsonlite_1.6              
#>  [87] XVector_0.24.0             bit_1.1-14                
#>  [89] mime_0.6                   gridExtra_2.3             
#>  [91] gplots_3.0.1.1             Rsamtools_2.0.0           
#>  [93] zingeR_0.1.0               stringi_1.4.3             
#>  [95] gmodels_2.18.1             gsl_2.1-6                 
#>  [97] bitops_1.0-6               cli_1.1.0                 
#>  [99] Rdpack_0.11-0              maps_3.3.0                
#> [101] RSQLite_2.1.1              pheatmap_1.0.12           
#> [103] data.table_1.12.2          DEDS_1.58.0               
#> [105] energy_1.7-5               rstudioapi_0.10           
#> [107] GenomicAlignments_1.20.0   qvalue_2.16.0             
#> [109] scran_1.12.1               fastcluster_1.1.25        
#> [111] scone_1.8.0                locfit_1.5-9.1            
#> [113] listenv_0.7.0              cobs_1.3-3                
#> [115] R.oo_1.22.0                prabclus_2.3-1            
#> [117] segmented_0.5-4.0          readxl_1.3.1              
#> [119] ROTS_1.12.0                munsell_0.5.0             
#> [121] cellranger_1.1.0           R.methodsS3_1.7.1         
#> [123] moments_0.14               hwriter_1.3.2             
#> [125] caTools_1.17.1.2           codetools_0.2-16          
#> [127] coda_0.19-2                vipor_0.4.5               
#> [129] lmtest_0.9-37              htmlTable_1.13.1          
#> [131] rARPACK_0.11-0             lsei_1.2-0                
#> [133] xtable_1.8-4               SAVER_1.1.1               
#> [135] ROCR_1.0-7                 diptest_0.75-7            
#> [137] lpsymphony_1.12.0          abind_1.4-5               
#> [139] FNN_1.1.3                  RANN_2.6.1                
#> [141] sparsesvd_0.1-4            CompQuadForm_1.4.3        
#> [143] bibtex_0.4.2               ggdendro_0.1-20           
#> [145] cluster_2.0.9              future.apply_1.2.0        
#> [147] Seurat_3.0.1               Matrix_1.2-17             
#> [149] prettyunits_1.0.2          shinyBS_0.61              
#> [151] lubridate_1.7.4.9000       ggridges_0.5.1            
#> [153] NOISeq_2.28.0              shinydashboard_0.7.1      
#> [155] mclust_5.4.3               igraph_1.2.4.1            
#> [157] slam_0.1-45                testthat_2.1.1            
#> [159] htmltools_0.3.6            yaml_2.2.0                
#> [161] GenomicFeatures_1.36.1     plotly_4.9.0              
#> [163] XML_3.98-1.20              DrImpute_1.0              
#> [165] foreign_0.8-71             withr_2.1.2               
#> [167] fitdistrplus_1.0-14        aroma.light_3.14.0        
#> [169] bit64_0.9-7                rngtools_1.3.1.1          
#> [171] doRNG_1.7.1                robustbase_0.93-5         
#> [173] outliers_0.14              scde_2.12.0               
#> [175] Biostrings_2.52.0          combinat_0.0-8            
#> [177] rsvd_1.0.1                 iCOBRA_1.12.1             
#> [179] memoise_1.1.0              evaluate_0.14             
#> [181] VGAM_1.1-1                 nonnest2_0.5-2            
#> [183] geneplotter_1.62.0         permute_0.9-5             
#> [185] fdrtool_1.2.15             acepack_1.4.1             
#> [187] checkmate_1.9.3            npsurv_0.4-0              
#> [189] tensorA_0.36.1             DECENT_1.0.0              
#> [191] ellipse_0.4.1              rjson_0.2.20              
#> [193] ggrepel_0.8.1              distillery_1.0-6          
#> [195] scDD_1.8.0                 stabledist_0.7-1          
#> [197] Lmoments_1.3-1             tools_3.6.0               
#> [199] sandwich_2.5-1             magrittr_1.5              
#> [201] RCurl_1.95-4.12            pbivnorm_0.6.0            
#> [203] ape_5.3                    bayesm_3.1-1              
#> [205] xml2_1.2.0                 EBSeq_1.24.0              
#> [207] httr_1.4.0                 assertthat_0.2.1          
#> [209] rmarkdown_1.13             Rhdf5lib_1.6.0            
#> [211] boot_1.3-22                globals_0.12.4            
#> [213] R6_2.4.0                   nnet_7.3-12               
#> [215] progress_1.2.2             gtools_3.8.1              
#> [217] statmod_1.4.32             Rook_1.1-1                
#> [219] BiocSingular_1.0.0         rhdf5_2.28.0              
#> [221] splines_3.6.0              colorspace_1.4-1          
#> [223] amap_0.8-17                generics_0.0.2            
#> [225] NBPSeq_0.3.0               compositions_1.40-2       
#> [227] base64enc_0.1-3            baySeq_2.18.0             
#> [229] pillar_1.4.1               HSMMSingleCell_1.4.0      
#> [231] GenomeInfoDbData_1.2.1     plyr_1.8.4                
#> [233] extRemes_2.0-10            dotCall64_1.0-0           
#> [235] gtable_0.3.0               rvest_0.3.4               
#> [237] SCnorm_1.6.0               monocle_2.12.0            
#> [239] knitr_1.23                 RcppArmadillo_0.9.500.2.0 
#> [241] latticeExtra_0.6-28        biomaRt_2.40.0            
#> [243] ADGofTest_0.3              copula_0.999-19.1         
#> [245] Cairo_1.5-10               doParallel_1.0.14         
#> [247] pscl_1.5.2                 flexmix_2.3-15            
#> [249] quantreg_5.40              AnnotationDbi_1.46.0      
#> [251] broom_0.5.2                scales_1.0.0              
#> [253] arm_1.10-1                 backports_1.1.4           
#> [255] IHW_1.12.0                 densityClust_0.3          
#> [257] lme4_1.1-21                blme_1.0-4                
#> [259] brew_1.0-6                 hms_0.4.2                 
#> [261] DESeq_1.36.0               Rtsne_0.15                
#> [263] shiny_1.3.2                grid_3.6.0                
#> [265] numDeriv_2016.8-1.1        bbmle_1.0.20              
#> [267] lazyeval_0.2.2             dynamicTreeCut_1.63-1     
#> [269] Formula_1.2-3              tsne_0.1-3                
#> [271] blockmodeling_0.3.4        crayon_1.3.4              
#> [273] MAST_1.10.0                RUVSeq_1.18.0             
#> [275] viridis_0.5.1              rpart_4.1-15              
#> [277] zinbwave_1.6.0             compiler_3.6.0
```

Note that I’ve also only tried this on Ubuntu.

## References

<div id="refs" class="references">

<div id="ref-gerard2019data">

Gerard, David. 2019. “Data-Based RNA-Seq Simulations by Binomial
Thinning.”

</div>

</div>
