#################
## Simulations evaluating various factor analyses
## Created by: David Gerard
## Created on: 06/07/2019
## NOTE: nc (the number of threads to use in the parallelization)
##       should be specified in Makefile. If you are running this in
##       interactive mode, you need to modify the below code so that nc
##       is not specified by commandArgs().
#################
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

one_sim_iter <- function(obj, bigmat) {
  stopifnot(nrow(obj) == 1)
  set.seed(obj$seed)

  ## functions to do FA -------------------------------------------------------
  source("./code/fa_methods.R")
  source("./code/signal_funs.R")

  ## subsample bigmat ---------------------------------------------------------
  which_gene <- sample(seq_len(nrow(bigmat)), size = obj$ngene)
  which_samp <- sample(seq_len(ncol(bigmat)), size = obj$nsamp)
  mat <- bigmat[which_gene, which_samp]
  lmat <- log2(mat + 0.5)

  ## Estimate number of factors -----------------------------------------------
  nsv <- sva::num.sv(dat = lmat, mod = matrix(1, nrow = obj$nsamp))

  ## Get factors and loadings -------------------------------------------------
  fl <- gen_fl(nsamp = obj$nsamp,
               ngene = obj$ngene,
               pz    = obj$pz,
               lsd   = obj$lsd)

  ## Median proportion of variance explained ----------------------------------
  flprod <- tcrossprod(fl$loadings, fl$factors)
  mpve <- median(apply(flprod, 1, stats::var) / apply(lmat, 1, stats::var))

  ## Add factors and loadings to mat ------------------------------------------
  thout <- seqgendiff::thin_diff(mat         = mat,
                                 design_perm = fl$factors,
                                 coef_perm   = fl$loadings,
                                 target_cor  = matrix(obj$corval[[1]], nrow = 1))

  lmat2 <- log2(thout$mat + 0.5)

  ## Do the factor analyses after adding signal -------------------------------
  faafter <- list(
    pca   = get_pca(mat = lmat2, k = nsv + 1),
    ica   = get_ica(mat = lmat2, k = nsv + 1),
    ssvd  = get_ssvd(mat = lmat2, k = nsv + 1),
    peer  = get_peer(mat = lmat2, k = nsv + 1),
    flash = get_flashr(mat = lmat2, k = nsv + 1))

  ## Max correlation metric ---------------------------------------------------
  cordiff <- sapply(faafter, FUN = function(obj) {
    corabs <- sort(abs(cor(thout$designmat, obj$factors)), decreasing = TRUE)
    corabs[1] - corabs[2]
  })

  ## Angle between X3 and columnspace of estimate -----------------------------
  anglevec <- sapply(faafter, FUN = function(obj) {
    svout <- svd(obj$factors)
    smat <- svout$u[, svout$d > 10^-12, drop = FALSE] ## make sure projection matrix is not computationally singular.
    X3hat <- smat %*% solve(crossprod(smat)) %*% crossprod(smat, thout$designmat)
    X3 <- thout$designmat
    X3 <- X3 / sqrt(sum(X3 ^ 2))
    X3hat <- X3hat / sqrt(sum(X3hat ^ 2))
    acos(crossprod(X3, X3hat))
  })

  ## MSE between scaled added factor and scaled estimated factors -------------
  msevec <- sapply(faafter, FUN = function(obj) {
    min(colMeans((scale(obj$factors, center = FALSE) - scale(thout$designmat, center = FALSE)[, 1]) ^ 2))
  })

  ## Change names and combine for output --------------------------------------
  names(cordiff)  <- paste0("cordiff_", names(cordiff))
  names(anglevec) <- paste0("angle_", names(anglevec))
  names(msevec)   <- paste0("mse_", names(msevec))
  retvec <- c(cordiff, anglevec, msevec, nsv = nsv + 1, unlist(obj),
              mpve = mpve)

  return(retvec)
}

## Read in data and filter out low-expressed genes ----------------------------
suppressPackageStartupMessages(library(SummarizedExperiment))
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
bigmat    <- assay(musc)

## Set simulation parameters --------------------------------------------------
ngene     <- 1000             ## Number of genes
nsamp_vec <- c(6, 10, 20)     ## Number of samples
cor_list  <- list(c(0, 0), c(0.5, 0)) ## Target correlation between new factor and two old factors
prop_zero <- c(0, 0.5, 0.9)   ## Proportion of loadings that are 0.
load_sd   <- c(0.4, 0.7, 1)   ## Standard deviation of loadings.
itermax   <- 200              ## Number of iterations per condition
seedindex <- seq_len(itermax)

pardf <- expand.grid(seed = seedindex,
                     nsamp  = nsamp_vec,
                     ngene  = ngene,
                     pz     = prop_zero,
                     lsd    = load_sd,
                     corval = cor_list)

## Randomly permute rows of pardf to evenly
## distribute computation between cores ---------------------------------------
pardf <- pardf[sample(seq_len(nrow(pardf))), ]

## Run simulations ------------------------------------------------------------
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1) ## make sure cluster is set up.
simout <- foreach(i = seq_len(nrow(pardf)), .combine = rbind) %dopar% {
  one_sim_iter(obj = pardf[i, ], bigmat = bigmat)
}
stopCluster(cl)

## Save output ----------------------------------------------------------------
write.csv(x = simout,
          file = "./output/fa_sims/fa_results.csv",
          row.names = FALSE)













