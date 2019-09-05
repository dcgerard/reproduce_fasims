
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2019)

## Introduction

This repository contains the code and instructions needed to reproduce
all of the results from Gerard (2019).

If you find a bug, please create an
[issue](https://github.com/dcgerard/reproduce_fasims/issues).

## Instructions

1.  Download and Install PEER: <https://github.com/PMBio/peer/wiki>. I
    had to build the R package from the cloned repository, as installing
    it from the source package didn’t work for me.

2.  Download and install powsimR: <https://github.com/bvieth/powsimR>.
    
    To download the version I used, try:
    
    ``` r
    devtools::install_github("bvieth/powsimR", 
                             ref = "c70c819beeec35a21558cfb717fe91e7704be126")
    ```

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
                       "latex2exp",
                       "devtools"))
    BiocManager::install(c("SummarizedExperiment", 
                           "pcaMethods", 
                           "sva",
                           "edgeR",
                           "DESeq2", 
                           "limma"))
    ```

4.  Download and install flashr:
    <https://github.com/stephenslab/flashr>.
    
    To download the version I used, try:
    
    ``` r
    devtools::install_github("stephenslab/flashr", 
                             ref = "276bca27a634d53697bae550c0566b3d1f12dbfc")
    ```

5.  Download data from the GTEx portal: <https://gtexportal.org/home/>.
    You’ll need to place the following files in the “data/gtex” folder:
    
    1.  GTEx\_Analysis\_2016-01-15\_v7\_RNASeQCv1.1.8\_gene\_reads.gct
    2.  GTEx\_v7\_Annotations\_SubjectPhenotypesDS.txt
    3.  GTEx\_v7\_Annotations\_SampleAttributesDS.txt
    
    I cannot release these data, but it’s free to sign up for access to
    GTEx.

6.  Adjust the Makefile. Change the `nc` and `rexec` variables in the
    Makefile according to your local computing environment. For example,
    you would need to decrease `nc` if you have fewer than 12 CPU cores.
    You can check the number of CPU cores you have by typing the
    following in R:
    
    ``` r
    parallel::detectCores()
    ```

7.  Run `make`. To reproduce all of the results in the paper, just type
    in the terminal
    
    ``` bash
    make
    ```
    
    To just reproduce the results comparing the `powsimR` datasets to
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

8.  Get coffee. Though most of the results will be generated rather
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
#> R version 3.6.1 (2019-07-05)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.3 LTS
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
#>  [1] flashr_0.6-3                devtools_2.1.0             
#>  [3] usethis_1.5.1               powsimR_1.1.4              
#>  [5] gamlss.dist_5.1-4           MASS_7.3-51.4              
#>  [7] peer_1.0                    DESeq2_1.24.0              
#>  [9] edgeR_3.26.8                limma_3.40.6               
#> [11] sva_3.32.1                  genefilter_1.66.0          
#> [13] mgcv_1.8-28                 nlme_3.1-141               
#> [15] pcaMethods_1.76.0           SummarizedExperiment_1.14.1
#> [17] DelayedArray_0.10.0         BiocParallel_1.18.1        
#> [19] matrixStats_0.54.0          Biobase_2.44.0             
#> [21] GenomicRanges_1.36.0        GenomeInfoDb_1.20.0        
#> [23] IRanges_2.18.2              S4Vectors_0.22.0           
#> [25] BiocGenerics_0.30.0         latex2exp_0.4.0            
#> [27] ggthemes_4.2.0              doSNOW_1.0.18              
#> [29] snow_0.4-3                  iterators_1.0.12           
#> [31] foreach_1.4.7               elasticnet_1.1.1           
#> [33] lars_1.2                    sparsepca_0.1.2            
#> [35] ssvd_1.0                    fastICA_1.2-2              
#> [37] vroom_1.0.2                 BiocManager_1.30.4         
#> [39] forcats_0.4.0               stringr_1.4.0              
#> [41] dplyr_0.8.3                 purrr_0.3.2                
#> [43] readr_1.3.1                 tidyr_0.8.3                
#> [45] tibble_2.1.3                ggplot2_3.2.1              
#> [47] tidyverse_1.2.1             seqgendiff_1.1.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] mixtools_1.1.0             softImpute_1.4            
#>   [3] minpack.lm_1.2-1           pbapply_1.4-2             
#>   [5] lattice_0.20-38            haven_2.1.1               
#>   [7] vctrs_0.2.0                penalized_0.9-51          
#>   [9] blob_1.2.0                 survival_2.44-1.1         
#>  [11] later_0.8.0                nloptr_1.2.1              
#>  [13] DBI_1.0.0                  R.utils_2.9.0             
#>  [15] SingleCellExperiment_1.6.0 uwot_0.1.3                
#>  [17] Linnorm_2.8.0              dqrng_0.2.1               
#>  [19] zlibbioc_1.30.0            MatrixModels_0.4-1        
#>  [21] pspline_1.0-18             SDMTools_1.1-221.1        
#>  [23] htmlwidgets_1.3            mvtnorm_1.0-11            
#>  [25] future_1.14.0              UpSetR_1.4.0              
#>  [27] leiden_0.3.1               scater_1.12.2             
#>  [29] irlba_2.3.3                DEoptimR_1.0-8            
#>  [31] Rcpp_1.0.2                 KernSmooth_2.23-15        
#>  [33] DT_0.8                     promises_1.0.1            
#>  [35] gdata_2.18.0               DDRTree_0.1.5             
#>  [37] pkgload_1.0.2              vegan_2.5-6               
#>  [39] RcppParallel_4.4.3         Hmisc_4.2-0               
#>  [41] ShortRead_1.42.0           apcluster_1.4.8           
#>  [43] fs_1.3.1                   RSpectra_0.15-0           
#>  [45] msir_1.3.2                 mnormt_1.5-5              
#>  [47] digest_0.6.20              png_0.1-7                 
#>  [49] qlcMatrix_0.9.7            sctransform_0.2.0         
#>  [51] cowplot_1.0.0              glmnet_2.0-18             
#>  [53] pkgconfig_2.0.2            docopt_0.6.1              
#>  [55] DelayedMatrixStats_1.6.0   ggbeeswarm_0.6.0          
#>  [57] minqa_1.2.4                lavaan_0.6-5              
#>  [59] reticulate_1.13            beeswarm_0.2.3            
#>  [61] spam_2.2-2                 modeltools_0.2-22         
#>  [63] xfun_0.9                   zoo_1.8-6                 
#>  [65] tidyselect_0.2.5           ZIM_1.1.0                 
#>  [67] reshape2_1.4.3             kernlab_0.9-27            
#>  [69] ica_1.0-2                  pcaPP_1.9-73              
#>  [71] EDASeq_2.18.0              viridisLite_0.3.0         
#>  [73] rtracklayer_1.44.3         pkgbuild_1.0.5            
#>  [75] rlang_0.4.0                hexbin_1.27.3             
#>  [77] glue_1.3.1                 metap_1.1                 
#>  [79] RColorBrewer_1.1-2         registry_0.5-1            
#>  [81] modelr_0.1.5               fpc_2.2-3                 
#>  [83] pkgmaker_0.27              fields_9.8-6              
#>  [85] SparseM_1.77               gbRd_0.4-11               
#>  [87] httpuv_1.5.1               class_7.3-15              
#>  [89] BPSC_0.99.2                RMTstat_0.3               
#>  [91] BiocNeighbors_1.2.0        annotate_1.62.0           
#>  [93] jsonlite_1.6               XVector_0.24.0            
#>  [95] bit_1.1-14                 mime_0.7                  
#>  [97] gridExtra_2.3              gplots_3.0.1.1            
#>  [99] Rsamtools_2.0.0            zingeR_0.1.0              
#> [101] stringi_1.4.3              gmodels_2.18.1            
#> [103] processx_3.4.1             gsl_2.1-6                 
#> [105] bitops_1.0-6               cli_1.1.0                 
#> [107] Rdpack_0.11-0              maps_3.3.0                
#> [109] RSQLite_2.1.2              pheatmap_1.0.12           
#> [111] data.table_1.12.2          DEDS_1.58.0               
#> [113] energy_1.7-6               rstudioapi_0.10           
#> [115] GenomicAlignments_1.20.1   qvalue_2.16.0             
#> [117] scran_1.12.1               fastcluster_1.1.25        
#> [119] scone_1.8.0                locfit_1.5-9.1            
#> [121] listenv_0.7.0              cobs_1.3-3                
#> [123] R.oo_1.22.0                prabclus_2.3-1            
#> [125] segmented_1.0-0            sessioninfo_1.1.1         
#> [127] readxl_1.3.1               ROTS_1.12.0               
#> [129] munsell_0.5.0              cellranger_1.1.0          
#> [131] R.methodsS3_1.7.1          moments_0.14              
#> [133] hwriter_1.3.2              caTools_1.17.1.2          
#> [135] codetools_0.2-16           coda_0.19-3               
#> [137] vipor_0.4.5                lmtest_0.9-37             
#> [139] htmlTable_1.13.1           rARPACK_0.11-0            
#> [141] lsei_1.2-0                 xtable_1.8-4              
#> [143] SAVER_1.1.1                ROCR_1.0-7                
#> [145] diptest_0.75-7             lpsymphony_1.12.0         
#> [147] abind_1.4-5                FNN_1.1.3                 
#> [149] RANN_2.6.1                 sparsesvd_0.2             
#> [151] CompQuadForm_1.4.3         bibtex_0.4.2              
#> [153] RcppAnnoy_0.0.12           ggdendro_0.1-20           
#> [155] cluster_2.1.0              future.apply_1.3.0        
#> [157] Seurat_3.1.0               zeallot_0.1.0             
#> [159] Matrix_1.2-17              SQUAREM_2017.10-1         
#> [161] prettyunits_1.0.2          shinyBS_0.61              
#> [163] lubridate_1.7.4.9000       ggridges_0.5.1            
#> [165] NOISeq_2.28.0              shinydashboard_0.7.1      
#> [167] mclust_5.4.5               igraph_1.2.4.1            
#> [169] remotes_2.1.0              slam_0.1-45               
#> [171] ashr_2.2-38                testthat_2.2.1            
#> [173] htmltools_0.3.6            yaml_2.2.0                
#> [175] GenomicFeatures_1.36.4     plotly_4.9.0              
#> [177] XML_3.98-1.20              DrImpute_1.0              
#> [179] foreign_0.8-72             withr_2.1.2               
#> [181] fitdistrplus_1.0-14        aroma.light_3.14.0        
#> [183] bit64_0.9-7                rngtools_1.4              
#> [185] doRNG_1.7.1                robustbase_0.93-5         
#> [187] outliers_0.14              scde_2.12.0               
#> [189] Biostrings_2.52.0          combinat_0.0-8            
#> [191] rsvd_1.0.2                 iCOBRA_1.12.1             
#> [193] memoise_1.1.0              evaluate_0.14             
#> [195] VGAM_1.1-1                 nonnest2_0.5-2            
#> [197] callr_3.3.1                geneplotter_1.62.0        
#> [199] permute_0.9-5              ps_1.3.0                  
#> [201] fdrtool_1.2.15             acepack_1.4.1             
#> [203] checkmate_1.9.4            desc_1.2.0                
#> [205] npsurv_0.4-0               truncnorm_1.0-8           
#> [207] tensorA_0.36.1             DECENT_1.0.0              
#> [209] ellipse_0.4.1              rjson_0.2.20              
#> [211] ggrepel_0.8.1              distillery_1.0-6          
#> [213] rprojroot_1.3-2            scDD_1.8.0                
#> [215] stabledist_0.7-1           Lmoments_1.3-1            
#> [217] tools_3.6.1                sandwich_2.5-1            
#> [219] magrittr_1.5               RCurl_1.95-4.12           
#> [221] pbivnorm_0.6.0             ape_5.3                   
#> [223] bayesm_3.1-3               xml2_1.2.2                
#> [225] EBSeq_1.24.0               httr_1.4.1                
#> [227] assertthat_0.2.1           rmarkdown_1.15            
#> [229] Rhdf5lib_1.6.0             boot_1.3-23               
#> [231] globals_0.12.4             R6_2.4.0                  
#> [233] nnet_7.3-12                progress_1.2.2            
#> [235] gtools_3.8.1               statmod_1.4.32            
#> [237] Rook_1.1-1                 BiocSingular_1.0.0        
#> [239] rhdf5_2.28.0               splines_3.6.1             
#> [241] colorspace_1.4-1           amap_0.8-17               
#> [243] generics_0.0.2             NBPSeq_0.3.0              
#> [245] compositions_1.40-2        base64enc_0.1-3           
#> [247] baySeq_2.18.0              mixsqp_0.1-97             
#> [249] pillar_1.4.2               HSMMSingleCell_1.4.0      
#> [251] GenomeInfoDbData_1.2.1     plyr_1.8.4                
#> [253] extRemes_2.0-10            dotCall64_1.0-0           
#> [255] gtable_0.3.0               rvest_0.3.4               
#> [257] SCnorm_1.6.0               monocle_2.12.0            
#> [259] knitr_1.24                 RcppArmadillo_0.9.700.2.0 
#> [261] latticeExtra_0.6-28        biomaRt_2.40.4            
#> [263] ADGofTest_0.3              copula_0.999-19.1         
#> [265] Cairo_1.5-10               doParallel_1.0.15         
#> [267] pscl_1.5.2                 flexmix_2.3-15            
#> [269] quantreg_5.51              AnnotationDbi_1.46.1      
#> [271] broom_0.5.2                scales_1.0.0              
#> [273] arm_1.10-1                 backports_1.1.4           
#> [275] IHW_1.12.0                 densityClust_0.3          
#> [277] lme4_1.1-21                blme_1.0-4                
#> [279] brew_1.0-6                 hms_0.5.1                 
#> [281] DESeq_1.36.0               Rtsne_0.15                
#> [283] shiny_1.3.2                grid_3.6.1                
#> [285] numDeriv_2016.8-1.1        bbmle_1.0.20              
#> [287] lazyeval_0.2.2             dynamicTreeCut_1.63-1     
#> [289] Formula_1.2-3              tsne_0.1-3                
#> [291] blockmodeling_0.3.4        crayon_1.3.4              
#> [293] MAST_1.10.0                RUVSeq_1.18.0             
#> [295] viridis_0.5.1              rpart_4.1-15              
#> [297] zinbwave_1.6.0             compiler_3.6.1
```

Note that I’ve also only tried this on Ubuntu.

## References

<div id="refs" class="references">

<div id="ref-gerard2019data">

Gerard, David. 2019. “Data-Based RNA-Seq Simulations by Binomial
Thinning.” *bioRxiv*. Cold Spring Harbor Laboratory.
<https://doi.org/10.1101/758524>.

</div>

</div>
