
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Reproduce the Results of Gerard (2019)

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

5.  Run `make`. To reproduce all of the results in the paper, just type
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
