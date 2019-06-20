######
## factor analysis methods used, with standardized output.
######

# library(pcaMethods)
# library(fastICA)
# library(elasticnet)
# library(peer)
# library(sparsepca)
# library(SummarizedExperiment)
# musc <- readRDS(file = "../output/tissue_data/muscle_skeletal.RDS")
# which_bad <- rowMeans(assay(musc)) < 10
# musc <- musc[!which_bad, ]
# mat <- log2(assay(musc)[1:1000, 1:10] + 0.5)


get_pca <- function(mat, k) {
  pcout <- pcaMethods::pca(object = t(mat), method = "svd", nPcs = k)
  retlist <- list(factors = pcaMethods::scores(pcout),
                  loadings = pcaMethods::loadings(pcout))
  return(retlist)
}

get_ica <- function(mat, k) {
  trash <- capture.output({
    icout <- fastICA::fastICA(X = t(mat), n.comp = k, method = "C", row.norm = TRUE)
  })
  retlist <- list(factors = icout$S,
                  loadings = t(icout$A))
  return(retlist)
}

get_spca <- function(mat, k) {
  mat <- mat - rowMeans(mat)
  vvec <- eigen(crossprod(mat))$values
  vvec <- vvec / sum(vvec)
  parval <- 10^-6
  while (TRUE) {
    parval <- parval + 1
    spout <- elasticnet::arrayspc(x    = t(mat),
                                  K    = k,
                                  para = rep(parval, length.out = k))
    if (sum(spout$pev) / sum(vvec[seq_len(k)]) < 0.9) {
      break
    }
  }

  retlist <- list(factors = crossprod(mat, spout$loadings),
                  loadings = spout$loadings)
  return(retlist)
}

get_peer <- function(mat, k) {
  model = peer::PEER()
  peer::PEER_setPhenoMean(model, t(mat))
  peer::PEER_setNk(model, k)
  peer::PEER_setAdd_mean(model, TRUE)
  peer::PEER_update(model)
  ## Get rid of intercept term
  retlist <- list(factors = peer::PEER_getX(model)[, -1, drop = FALSE],
                  loadings = peer::PEER_getW(model)[, -1, drop = FALSE])
  return(retlist)
}

get_ssvd <- function(mat, k) {
  sout <- ssvd::ssvd(x = mat - rowMeans(mat), method = "theory", r = k)
  retlist <- list(factors = sout$v,
                  loadings = sout$u)
  return(retlist)
}

get_spca <- function(mat, k, alpha) {
  spout <- sparsepca::spca(X = t(mat - rowMeans(mat)), k = k, alpha = alpha)
  retlist <- list(factors = spout$scores,
                  loadings = spout$loadings)
  return(alpha)
}




get_flashr <- function(mat, k) {
  fdat <- flashr::flash_set_data(mat - rowMeans(mat))

  fout <- flashr::flash(data     = fdat,
                        Kmax     = k,
                        var_type = "by_column",
                        verbose  = FALSE,
                        backfit  = TRUE)

  fitted_val <- flashr:::flash_get_ldf(f = fout, drop_zero_factors = FALSE)

  retlist <- list(factors = fitted_val$f,
                  loadings = fitted_val$l)

  return(retlist)
}




















