################################
## Estimate parameters for powsimR simulations
################################

suppressPackageStartupMessages(library(SummarizedExperiment))
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
epout <- powsimR::estimateParam(countData = assay(musc),
                                Distribution = "NB",
                                Normalisation = "TMM",
                                RNAseq = "bulk")
saveRDS(object = epout, file = "./output/compare_powsimR/powsim_params.RDS")
