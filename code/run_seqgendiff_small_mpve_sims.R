##########################
## Copy of run_seqgendiff_sims.R but with a smaller log-fold change to
## decrease the MPVE. Also measure time.
##########################

## load packages --------------------------------------------------------------
library(seqgendiff)
suppressPackageStartupMessages(library(SummarizedExperiment))
source("./code/de_methods.R")
suppressPackageStartupMessages(library(doSNOW))

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
ash_df <- read.csv("./output/simseq_vs_seqgendiff/ash_ests.csv")

prop_null <- 0.9
lfc_sd    <- ash_df$sdhat
nsamp     <- 10
ngene     <- 10000
itermax   <- 500
fdr_control <- 0.05

## Load in muscle data and filter ---------------------------------------------
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
fullcounts <- assay(musc)

## Set up parallel computing environment --------------------------------------
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

retmat <- foreach (iterindex = seq_len(itermax),
         .combine = rbind,
         .export = c("thin_2group")) %dopar% {


           set.seed(iterindex)
           ## Simulate data ---------------------------------------------------
           time_run <- system.time({
             which_gene <- sort(sample(seq_len(BiocGenerics::nrow(musc)), ngene))
             which_samp <- sort(sample(seq_len(BiocGenerics::ncol(musc)), nsamp))
             submusc   <- fullcounts[which_gene, which_samp]
             thout <- thin_2group(mat = submusc,
                                  prop_null = prop_null,
                                  signal_fun = stats::rnorm,
                                  signal_params = list(mean = 0, sd = lfc_sd),
                                  group_prop = 0.5)
           })
           countdat <- thout$mat
           design_mat <- cbind(thout$design_obs, thout$designmat)
           beta <- c(thout$coefmat)
           which_null <- abs(beta) < 10^-6

           # Fit methods ----------------------------------------------------
           fitlist <- list(
             vout = get_voom(countdat = countdat, design_mat = design_mat),
             dout = get_DESeq2(countdat = countdat, design_mat = design_mat),
             eout = get_edgeR(countdat = countdat, design_mat = design_mat)
           )

           ## Assess fits -----------------------------------------------------
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

           msevec <- sapply(fitlist, FUN = function(obj) {
             mean((obj$bhat - beta)^2, na.rm = TRUE)
           })
           names(msevec) <- paste0("mse_", names(msevec))

           ## Summary stat of count matrix ------------------------------------
           varvec <- apply(log2(countdat + 0.5)[!which_null, , drop = FALSE], 1, var)
           betavarvec <- apply(tcrossprod(beta[!which_null], design_mat[, 2]), 1, var)
           mpve <- median(betavarvec / varvec)

           ## Return ----------------------------------------------------------
           retvec <- c(fprvec, powervec, msevec, mpve = mpve, time = time_run[["elapsed"]])

           retvec
}
stopCluster(cl)

saveRDS(object = retmat, file = "./output/diff_exp_out/seqgendiff_small_mpve_sims.RDS")


