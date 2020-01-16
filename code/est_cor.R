#############################
## Simulations to evaluate effective_cor()
## Created by: David Gerard
## Created on: 06/10/2019
#############################
library(seqgendiff)
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


## Set up simulation settings ------------------------------------------------
itermax <- 100
nvec <- c(6, 10, 20, 500)
cor_list <- list(c(0, 0), c(0.5, 0), c(0.9, 0), c(0.5, 0.5))
pardf <- expand.grid(seed = seq_len(itermax),
                     nsamp = nvec,
                     target_cor = cor_list)

## Set up parallelization -----------------------------------------------------
cl <- parallel::makeCluster(nc)
doParallel::registerDoParallel(cl = cl)
stopifnot(foreach::getDoParWorkers() > 1)

## Run simulations ------------------------------------------------------------
retmat <- foreach(index = seq_len(nrow(pardf)),
        .combine = rbind) %dopar% {
          ## Get current parameters -------------------------------------------
          obj        <- pardf[index, ]
          target_cor <- matrix(obj$target_cor[[1]], ncol = 1)
          nsamp      <- obj$nsamp
          set.seed(obj$seed)

          ## Generate a two-group and a normal design -------------------------
          design_perm1 <- cbind(rep(c(0, 1), each = nsamp / 2),
                                rep(c(0, 1), length.out = nsamp))
          design_perm2 <- matrix(rnorm(2 * nsamp), ncol = 2)
          sv <- matrix(rnorm(nsamp), ncol = 1)

          ## Estimate correlations --------------------------------------------
          corest1 <- seqgendiff::effective_cor(design_perm = design_perm1,
                                               sv = sv,
                                               target_cor = target_cor,
                                               iternum = 100)
          corest2 <- seqgendiff::effective_cor(design_perm = design_perm2,
                                               sv = sv,
                                               target_cor = target_cor,
                                               iternum = 100)

          ## Return estimates and empiricals ----------------------------------
          retvec <- c(corest1, corest2)
          names(retvec) <- c("ce11", "ce12", "ce21", "ce22")
          retvec <- c(retvec, unlist(pardf[index, ]))
          retvec
        }

## Stop parallelization -------------------------------------------------------
stopCluster(cl)

## Save output ----------------------------------------------------------------
saveRDS(object = retmat, file = "./output/est_cor/corsimout.RDS")






