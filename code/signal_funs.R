################
## Functions for signal generation
################

## Function to generate loadings and factors ----------------------------------
gen_fl <- function(nsamp, ngene, pz, lsd) {
  lvec <- stats::rnorm(n = ngene)
  lvec <- scale(lvec) * lsd
  which_zero <- sample(seq_len(ngene), round(pz * ngene))
  lvec[which_zero] <- 0
  lvec <- matrix(lvec)
  fvec <- rnorm(n = nsamp)
  fvec <- matrix(scale(fvec))
  return(list(factors = fvec, loadings = lvec))
}
