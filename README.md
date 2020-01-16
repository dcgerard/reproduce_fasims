
<!-- README.md is generated from README.Rmd. Please edit that file -->
Reproduce the Results of Gerard (2019)
======================================

Introduction
------------

This repository contains the code and instructions needed to reproduce all of the results from Gerard (2019).

If you find a bug, please create an [issue](https://github.com/dcgerard/reproduce_fasims/issues).

Instructions
------------

1.  Download and Install PEER: <https://github.com/PMBio/peer/wiki>. I had to build the R package from the cloned repository, as installing it from the source package didn't work for me.

    When I installed it (in Ubuntu), I had to use an older version of g++ (see <https://stackoverflow.com/questions/1616983/building-r-packages-using-alternate-gcc>). The steps I used were

    1.  Install g++-5

        ``` bash
        sudo apt-get install g++-5
        ```

    2.  Create a Makevars file

        ``` bash
        touch ~/.R/Makevars
        ```

    3.  Edit the Makevars file to set the default g++ version to use

            CXX=g++-5

    4.  Install peer from the cloned repo using their instructions.

    5.  Delete the Makevars file so that future packages won't use the older version of g++.

        ``` bash
        rm ~/.R/Makevars
        ```

2.  Download and install powsimR: <https://github.com/bvieth/powsimR>.

    To download the version I used, try:

    ``` r
    devtools::install_github("bvieth/powsimR", 
                             ref = "c70c819beeec35a21558cfb717fe91e7704be126")
    ```

    I had to install a few packages using apt in my Ubuntu system before all of the dependencies would install:

    ``` bash
    sudo apt-get install libcairo2-dev
    sudo apt-get install libgsl-dev
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
                       "devtools",
                       "SimSeq"))
    BiocManager::install(c("SummarizedExperiment", 
                           "pcaMethods", 
                           "sva",
                           "edgeR",
                           "DESeq2", 
                           "limma",
                           "Seurat"))
    ```

4.  Download and install flashr: <https://github.com/stephenslab/flashr>.

    To download the version I used, try:

    ``` r
    devtools::install_github("stephenslab/ebnm")
    devtools::install_github("stephenslab/flashr", 
                             ref = "276bca27a634d53697bae550c0566b3d1f12dbfc")
    ```

5.  Download data from the GTEx portal (“GTEx Consortium” 2017): <https://gtexportal.org/home/>. You'll need to place the following files in the "data/gtex" folder:

    1.  GTEx\_Analysis\_2016-01-15\_v7\_RNASeQCv1.1.8\_gene\_reads.gct
    2.  GTEx\_v7\_Annotations\_SubjectPhenotypesDS.txt
    3.  GTEx\_v7\_Annotations\_SampleAttributesDS.txt

    I cannot release these data, but it's free to sign up for access to GTEx.

6.  Download the 10x Genomics PBMC dataset (Zheng et al. 2017).

    We used the same dataset as from Seurat tutorial (<https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html>). The data can be directly downloaded from here:

    <https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz>

    After unzipping, you should then place the following files in the "data/pbmc" folder:

    1.  barcodes.tsv
    2.  genes.tsv
    3.  matrix.mtx

7.  Adjust the Makefile. Change the `nc` and `rexec` variables in the Makefile according to your local computing environment. For example, you would need to decrease `nc` if you have fewer than 12 CPU cores. You can check the number of CPU cores you have by typing the following in R:

    ``` r
    parallel::detectCores()
    ```

8.  Run `make`. To reproduce all of the results in the paper, just type in the terminal

    ``` bash
    make
    ```

    To just reproduce the results comparing the `powsimR` datasets to the `seqgendiff` datasets, run in the terminal:

    ``` bash
    make powsimr
    ```

    To just reproduce the differential expression simulations, run in the terminal:

    ``` bash
    make diff_exp
    ```

    To just reproduce the factor analysis simulations, run in the terminal:

    ``` bash
    make FAsims
    ```

    To just reproduce the plot of the flexible class of mixtures of binomials and negative binomials, run in the terminal:

    ``` bash
    make NBplots
    ```

    To just reproduce the simulations evaluating the Monte Carlo correlation estimator, run in the terminal:

    ``` bash
    make corest
    ```

    To just reproduce the simulations using the PBMC data, run in the terminal

    ``` bash
    make sc_fa
    ```

9.  Get coffee. Though most of the results will be generated rather quickly, the `make FAsims` call will take a very long time (mostly because of PEER). You should get some coffee! Here is a list of some of my favorite places:

    -   Washington, DC
        -   [Colony Club](https://www.yelp.com/biz/colony-club-washington)
        -   [Grace Street Coffee](https://www.yelp.com/biz/grace-street-coffee-georgetown)
        -   [Maketto](https://www.yelp.com/biz/maketto-washington-2)
    -   Chicago
        -   [Sawada Coffee](https://www.yelp.com/biz/sawada-coffee-chicago)
        -   [Plein Air Cafe](https://www.yelp.com/biz/plein-air-cafe-and-eatery-chicago-2)
    -   Seattle
        -   [Bauhaus Ballard](https://www.yelp.com/biz/bauhaus-ballard-seattle)
        -   [Cafe Solstice](https://www.yelp.com/biz/cafe-solstice-seattle)
    -   Columbus
        -   [Yeah, Me Too](https://www.yelp.com/biz/yeah-me-too-columbus)
        -   [Stauf's Coffee Roasters](https://www.yelp.com/biz/staufs-coffee-roasters-columbus-2)

Package Versions
----------------

If you are having trouble reproducing these results, check your package versions. These are the ones that I used:

``` r
sessionInfo()
#> R version 3.6.2 (2019-12-12)
#> Platform: x86_64-pc-linux-gnu (64-bit)
#> Running under: Ubuntu 18.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
#> LAPACK: /home/dgerard/.local/share/r-miniconda/envs/r-reticulate/lib/libmkl_rt.so
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> attached base packages:
#> [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] Seurat_3.1.2                SimSeq_1.4.0               
#>  [3] flashr_0.6-7                devtools_2.2.1             
#>  [5] usethis_1.5.1               powsimR_1.1.4              
#>  [7] gamlss.dist_5.1-5           MASS_7.3-51.5              
#>  [9] peer_1.0                    DESeq2_1.26.0              
#> [11] edgeR_3.28.0                limma_3.42.0               
#> [13] sva_3.34.0                  genefilter_1.68.0          
#> [15] mgcv_1.8-31                 nlme_3.1-143               
#> [17] pcaMethods_1.78.0           SummarizedExperiment_1.16.1
#> [19] DelayedArray_0.12.2         BiocParallel_1.20.1        
#> [21] matrixStats_0.55.0          Biobase_2.46.0             
#> [23] GenomicRanges_1.38.0        GenomeInfoDb_1.22.0        
#> [25] IRanges_2.20.1              S4Vectors_0.24.1           
#> [27] BiocGenerics_0.32.0         latex2exp_0.4.0            
#> [29] ggthemes_4.2.0              doSNOW_1.0.18              
#> [31] snow_0.4-3                  iterators_1.0.12           
#> [33] foreach_1.4.7               elasticnet_1.1.1           
#> [35] lars_1.2                    sparsepca_0.1.2            
#> [37] ssvd_1.0                    fastICA_1.2-2              
#> [39] vroom_1.2.0                 BiocManager_1.30.10        
#> [41] forcats_0.4.0               stringr_1.4.0              
#> [43] dplyr_0.8.3                 purrr_0.3.3                
#> [45] readr_1.3.1                 tidyr_1.0.0                
#> [47] tibble_2.1.3                ggplot2_3.2.1              
#> [49] tidyverse_1.3.0             seqgendiff_1.2.0           
#> 
#> loaded via a namespace (and not attached):
#>   [1] spam_2.5-1                 locfit_1.5-9.1            
#>   [3] remotes_2.1.0              slam_0.1-47               
#>   [5] ZIM_1.1.0                  ashr_2.2-39               
#>   [7] lattice_0.20-38            EDASeq_2.20.0             
#>   [9] vctrs_0.2.1                blob_1.2.0                
#>  [11] R.oo_1.23.0                withr_2.1.2               
#>  [13] foreign_0.8-74             registry_0.5-1            
#>  [15] readxl_1.3.1               lifecycle_0.1.0           
#>  [17] pkgmaker_0.27              cellranger_1.1.0          
#>  [19] munsell_0.5.0              SCnorm_1.8.2              
#>  [21] codetools_0.2-16           SparseM_1.78              
#>  [23] lmtest_0.9-37              gbRd_0.4-11               
#>  [25] RMTstat_0.3                densityClust_0.3          
#>  [27] annotate_1.64.0            apcluster_1.4.8           
#>  [29] fs_1.3.1                   brew_1.0-6                
#>  [31] ellipse_0.4.1              Rtsne_0.15                
#>  [33] DESeq_1.38.0               stringi_1.4.5             
#>  [35] qlcMatrix_0.9.7            distillery_1.0-6          
#>  [37] sctransform_0.2.1          grid_3.6.2                
#>  [39] sandwich_2.5-1             cluster_2.1.0             
#>  [41] blockmodeling_0.3.4        ape_5.3                   
#>  [43] pbivnorm_0.6.0             pkgconfig_2.0.3           
#>  [45] pheatmap_1.0.12            SQUAREM_2020.1            
#>  [47] prettyunits_1.1.0          data.table_1.12.8         
#>  [49] ggridges_0.5.2             lubridate_1.7.4           
#>  [51] httr_1.4.1                 igraph_1.2.4.2            
#>  [53] progress_1.2.2             fastcluster_1.1.25        
#>  [55] scone_1.10.0               modeltools_0.2-22         
#>  [57] haven_2.2.0                amap_0.8-18               
#>  [59] htmltools_0.4.0            viridisLite_0.3.0         
#>  [61] yaml_2.2.0                 baySeq_2.20.0             
#>  [63] pillar_1.4.3               hexbin_1.28.0             
#>  [65] later_1.0.0                fitdistrplus_1.0-14       
#>  [67] glue_1.3.1                 DBI_1.1.0                 
#>  [69] doRNG_1.7.1                plyr_1.8.5                
#>  [71] robustbase_0.93-5          outliers_0.14             
#>  [73] dotCall64_1.0-0            gtable_0.3.0              
#>  [75] rsvd_1.0.2                 caTools_1.17.1.4          
#>  [77] nonnest2_0.5-2             latticeExtra_0.6-29       
#>  [79] fastmap_1.0.1              AnnotationDbi_1.48.0      
#>  [81] broom_0.5.3                rARPACK_0.11-0            
#>  [83] SAVER_1.1.2                arm_1.10-1                
#>  [85] checkmate_1.9.4            promises_1.1.0            
#>  [87] IHW_1.14.0                 truncnorm_1.0-8           
#>  [89] lpsymphony_1.14.0          FNN_1.1.3                 
#>  [91] mnormt_1.5-5               hms_0.5.3                 
#>  [93] askpass_1.1                png_0.1-7                 
#>  [95] lazyeval_0.2.2             Formula_1.2-3             
#>  [97] crayon_1.3.4               bayesm_3.1-4              
#>  [99] RUVSeq_1.20.0              reprex_0.3.0              
#> [101] boot_1.3-24                softImpute_1.4            
#> [103] tidyselect_0.2.5           xfun_0.12                 
#> [105] BiocSingular_1.2.1         TFisher_0.2.0             
#> [107] kernlab_0.9-29             splines_3.6.2             
#> [109] survival_3.1-8             penalized_0.9-51          
#> [111] rappdirs_0.3.1             bit64_0.9-7               
#> [113] segmented_1.1-0            rngtools_1.4              
#> [115] modelr_0.1.5               jpeg_0.1-8.1              
#> [117] extRemes_2.0-11            scde_2.14.0               
#> [119] ROTS_1.14.0                fields_10.0               
#> [121] MatrixModels_0.4-1         combinat_0.0-8            
#> [123] R.methodsS3_1.7.1          SDMTools_1.1-221.2        
#> [125] VGAM_1.1-2                 UpSetR_1.4.0              
#> [127] permute_0.9-5              pscl_1.5.2                
#> [129] quantreg_5.54              fdrtool_1.2.15            
#> [131] htmlTable_1.13.3           BPSC_0.99.2               
#> [133] xtable_1.8-4               DT_0.11                   
#> [135] DDRTree_0.1.5              gdata_2.18.0              
#> [137] vegan_2.5-6                abind_1.4-5               
#> [139] ShortRead_1.44.1           mime_0.8                  
#> [141] tensorA_0.36.1             rjson_0.2.20              
#> [143] ggrepel_0.8.1              sparsesvd_0.2             
#> [145] processx_3.4.1             numDeriv_2016.8-1.1       
#> [147] bibtex_0.4.2.2             tools_3.6.2               
#> [149] Rdpack_0.11-1              cli_2.0.1                 
#> [151] magrittr_1.5               future.apply_1.4.0        
#> [153] Matrix_1.2-18              shinyBS_0.61              
#> [155] DEDS_1.60.0                DelayedMatrixStats_1.8.0  
#> [157] ggbeeswarm_0.6.0           assertthat_0.2.1          
#> [159] qvalue_2.18.0              mixtools_1.1.0            
#> [161] ica_1.0-2                  pbapply_1.4-2             
#> [163] NBPSeq_0.3.0               BiocFileCache_1.10.2      
#> [165] rtracklayer_1.46.0         plotly_4.9.1              
#> [167] R.utils_2.9.2              metap_1.2                 
#> [169] HSMMSingleCell_1.6.0       multcomp_1.4-12           
#> [171] zlibbioc_1.32.0            hwriter_1.3.2             
#> [173] monocle_2.14.0             biomaRt_2.42.0            
#> [175] geneplotter_1.64.0         ps_1.3.0                  
#> [177] mutoss_0.1-12              fansi_0.4.1               
#> [179] lsei_1.2-0                 TH.data_1.0-10            
#> [181] ROCR_1.0-7                 KernSmooth_2.23-16        
#> [183] backports_1.1.5            XVector_0.26.0            
#> [185] bit_1.1-15                 Rsamtools_2.2.1           
#> [187] gplots_3.0.1.2             RANN_2.6.1                
#> [189] zingeR_0.1.0               shiny_1.4.0               
#> [191] RcppAnnoy_0.0.14           maps_3.3.0                
#> [193] viridis_0.5.1              rstudioapi_0.10           
#> [195] minqa_1.2.4                Rhdf5lib_1.8.0            
#> [197] lavaan_0.6-5               gtools_3.8.1              
#> [199] beeswarm_0.2.3             Rook_1.1-1                
#> [201] listenv_0.8.0              reshape2_1.4.3            
#> [203] rhdf5_2.30.1               generics_0.0.2            
#> [205] colorspace_1.4-1           base64enc_0.1-3           
#> [207] compositions_1.40-3        GenomicFeatures_1.38.0    
#> [209] cobs_1.3-3                 XML_3.98-1.20             
#> [211] pkgbuild_1.0.6             DrImpute_1.0              
#> [213] sn_1.5-4                   aroma.light_3.16.0        
#> [215] dbplyr_1.4.2               RColorBrewer_1.1-2        
#> [217] Linnorm_2.10.0             dqrng_0.2.1               
#> [219] GenomeInfoDbData_1.2.2     Biostrings_2.54.0         
#> [221] moments_0.14               evaluate_0.14             
#> [223] memoise_1.1.0              iCOBRA_1.14.0             
#> [225] coda_0.19-3                knitr_1.26                
#> [227] doParallel_1.0.15          vipor_0.4.5               
#> [229] httpuv_1.5.2               class_7.3-15              
#> [231] irlba_2.3.3                Rcpp_1.0.3                
#> [233] acepack_1.4.1              openssl_1.4.1             
#> [235] diptest_0.75-7             pkgload_1.0.2             
#> [237] jsonlite_1.6               Hmisc_4.3-0               
#> [239] RSpectra_0.16-0            msir_1.3.2                
#> [241] digest_0.6.23              gmodels_2.18.1            
#> [243] CompQuadForm_1.4.3         scDD_1.10.0               
#> [245] rprojroot_1.3-2            cowplot_1.0.0             
#> [247] Lmoments_1.3-1             bitops_1.0-6              
#> [249] RSQLite_2.2.0              tsne_0.1-3                
#> [251] EBSeq_1.26.0               NOISeq_2.30.0             
#> [253] rmarkdown_2.0              globals_0.12.5            
#> [255] compiler_3.6.2             nnet_7.3-12               
#> [257] multtest_2.42.0            reticulate_1.14           
#> [259] statmod_1.4.33             scran_1.14.5              
#> [261] zoo_1.8-7                  minpack.lm_1.2-1          
#> [263] testthat_2.3.1             rlang_0.4.2               
#> [265] mixsqp_0.2-2               nloptr_1.2.1              
#> [267] prabclus_2.3-2             uwot_0.1.5                
#> [269] SingleCellExperiment_1.8.0 sessioninfo_1.1.1         
#> [271] fpc_2.2-3                  rvest_0.3.5               
#> [273] bdsmatrix_1.3-4            future_1.15.1             
#> [275] mvtnorm_1.0-12             htmlwidgets_1.5.1         
#> [277] RcppArmadillo_0.9.800.3.0  callr_3.4.0               
#> [279] leiden_0.3.1               Cairo_1.5-10              
#> [281] flexmix_2.3-15             curl_4.3                  
#> [283] scater_1.14.6              DEoptimR_1.0-8            
#> [285] BiocNeighbors_1.4.1        scales_1.1.0              
#> [287] plotrix_3.7-7              desc_1.2.0                
#> [289] npsurv_0.4-0               RcppParallel_4.4.4        
#> [291] lme4_1.1-21                gridExtra_2.3             
#> [293] DECENT_1.1.0               bbmle_1.0.22              
#> [295] RCurl_1.95-4.12            ggdendro_0.1-20           
#> [297] zeallot_0.1.0              docopt_0.6.1              
#> [299] ellipsis_0.3.0             MAST_1.12.0               
#> [301] xml2_1.2.2                 shinydashboard_0.7.1      
#> [303] rpart_4.1-15               R6_2.4.1                  
#> [305] mclust_5.4.5               zinbwave_1.8.0            
#> [307] GenomicAlignments_1.22.1
```

Note that I've also only tried this on Ubuntu.

References
----------

Gerard, David. 2019. “Data-Based RNA-Seq Simulations by Binomial Thinning.” *bioRxiv*. Cold Spring Harbor Laboratory. doi:[10.1101/758524](https://doi.org/10.1101/758524).

“GTEx Consortium”. 2017. “Genetic Effects on Gene Expression Across Human Tissues.” *Nature* 550 (7675). Nature Publishing Group: 204. doi:[10.1038/nature24277](https://doi.org/10.1038/nature24277).

Zheng, Grace XY, Jessica M Terry, Phillip Belgrader, Paul Ryvkin, Zachary W Bent, Ryan Wilson, Solongo B Ziraldo, et al. 2017. “Massively Parallel Digital Transcriptional Profiling of Single Cells.” *Nature Communications* 8. Nature Publishing Group: 14049. doi:[10.1038/ncomms14049](https://doi.org/10.1038/ncomms14049).
