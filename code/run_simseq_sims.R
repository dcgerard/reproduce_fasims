########################################
## Simulate data according to SimSeq and fit basic DE methods
########################################

## Load packages ---------------------------------------------------------------
library(SimSeq)
source("./code/de_methods.R")
suppressPackageStartupMessages(library(doSNOW))
suppressPackageStartupMessages(library(SummarizedExperiment))

# Number of threads to use for multithreaded computing. This must be
# specified in the command-line shell; e.g., to use 8 threads, run
# command
#
#  R CMD BATCH '--args nc=8' mouthwash_sims.R
#
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}

## Simulation Settings --------------------------------------------------------
prop_null <- 0.9
nsamp     <- 10
ngene     <- 10000
itermax   <- 500
fdr_control <- 0.05

## Load in muscle data and filter ---------------------------------------------
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
sexvec    <- as.factor(colData(musc)$SEX)
fullcounts <- assay(musc)

## Set up parallel computing environment --------------------------------------
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

retmat <- foreach (iterindex = seq_len(itermax),
                   .combine = rbind,
                   .export = c("SimData")) %dopar% {
                     
                     set.seed(iterindex)
                     
                     ## Simulate data -----------------------------------------
                     time_run <- system.time({
                       simout <- SimData(counts = fullcounts, 
                                         treatment = sexvec,
                                         sort.method = "unpaired",
                                         k.ind = nsamp/2, 
                                         n.genes = ngene,
                                         n.diff = round((1 - prop_null) * ngene))                       
                     })
                     
                     countdat <- simout$counts
                     design_mat <- cbind(1, simout$treatment)
                     which_null <- !simout$DE.ind
                     
                     ## Fit methods -------------------------------------------
                     fitlist <- list(
                       vout = get_voom(countdat = countdat, design_mat = design_mat),
                       dout = get_DESeq2(countdat = countdat, design_mat = design_mat),
                       eout = get_edgeR(countdat = countdat, design_mat = design_mat)
                     )
                     
                     ## Assess fits -------------------------------------------
                     fitlist <- lapply(fitlist, FUN = function(obj) {
                       obj$qval <- p.adjust(obj$pval, method = "BH")
                       obj$discovery <- obj$qval < fdr_control
                       return(obj)
                     })
                     
                     fprvec <- sapply(fitlist, FUN = function(obj) {
                       if (any(obj$discovery, na.rm = TRUE)) {
                         mean(which_null[obj$discovery], na.rm = TRUE)
                       } else {
                         0
                       }
                     })
                     names(fprvec) <- paste0("fpr_", names(fprvec))
                     
                     powervec <- sapply(fitlist, FUN = function(obj) {
                       mean(obj$discovery[!which_null], na.rm = TRUE)
                     })
                     names(powervec) <- paste0("power_", names(powervec))
                     
                     ## Summary stat of count matrix --------------------------
                     varvec <- apply(log2(countdat + 0.5)[!which_null, , drop = FALSE], 1, var)
                     betavarvec <- apply(tcrossprod(fitlist$vout$bhat[!which_null], design_mat[, 2]), 1, var)
                     mpve <- median(betavarvec / varvec, na.rm = TRUE)
                     
                     
                     retvec <- c(fprvec, powervec, mpve = mpve, time = time_run[["elapsed"]])
                     
                     retvec
                   }
stopCluster(cl)

saveRDS(object = retmat, file = "./output/diff_exp_out/simseq_sims.RDS")




