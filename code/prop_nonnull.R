suppressPackageStartupMessages(library(SummarizedExperiment))
library(limma)
library(ashr)

## Load in muscle data and filter ---------------------------------------------
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]
countdat  <- assay(musc)
design_mat <- model.matrix(~as.factor(colData(musc)$SEX))

vout <- limma::voom(counts = countdat, design = design_mat)
lout <- limma::lmFit(vout)
eout <- limma::eBayes(lout)

betahat <- eout$coefficients[, 2]
sevec <- eout$sigma * eout$stdev.unscaled[, 2]

aout <- ashr::ash(betahat = betahat, sebetahat = sevec)

pi0hat <- get_pi0(aout)

ghat <- get_fitted_g(aout)
varhat <- sum(ghat$b^2 * ghat$pi) * 2/3 ## variance of mixture of symmetric uniforms
sdhat <- sqrt(varhat)

# library(qvalue)
# qout <- qvalue::qvalue(p = eout$p.value[, 2])
# weights <- (1 - qout$lfdr)
# weights <- weights / sum(weights)
# meanabs <- sum(abs(betahat) * weights)

df <- data.frame(pi0hat = pi0hat, sdhat = sdhat)
write.csv(x = df, file = "./output/simseq_vs_seqgendiff/ash_ests.csv", row.names = FALSE)
