#################################
## Look at data features using powsimR and seqgendiff
#################################

## Load in muscle data and filter based on edgeR's criterions
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(SummarizedExperiment))
library(ggthemes)
library(latex2exp)
musc      <- readRDS("./output/tissue_data/muscle_skeletal.RDS")
which_bad <- rowMeans(assay(musc)) < 10
musc      <- musc[!which_bad, ]

## Set simulation parameters
set.seed(1)
prop_null <- 0.75
lfc_sd    <- 0.8
nsamp     <- ncol(musc)
ngene     <- nrow(musc)

## Simulate powsimR
epout <- readRDS("./output/compare_powsimR/powsim_params.RDS")
psout <- powsimR::simulateCounts(n = c(nsamp / 2, nsamp / 2),
                                 ngenes = ngene,
                                 p.DE = 1 - prop_null,
                                 params = epout,
                                 pLFC   = function(n) rnorm(n, mean = 0, sd = lfc_sd))
countdat_ps <- psout$GeneCounts
designmat_ps <- cbind(1, c(rep(0, nsamp / 2), rep(1, nsamp / 2)))
beta_ps <- c(psout$pLFC)
whichnull_ps <- abs(beta_ps) < 10^-6

## Simulate seqgendiff
which_gene <- sort(sample(seq_len(nrow(musc)), ngene))
which_samp <- sort(sample(seq_len(ncol(musc)), nsamp))
submusc   <- musc[which_gene, which_samp]
thout <- seqgendiff::thin_2group(mat = assay(submusc),
                                 prop_null = prop_null,
                                 signal_fun = stats::rnorm,
                                 signal_params = list(mean = 0, sd = lfc_sd),
                                 group_prop = 0.5)
countdat_sgd <- thout$mat
designmat_sgd <- cbind(thout$design_obs, thout$designmat)
beta_sgd <- c(thout$coefmat)
whichnull_sgd <- abs(beta_sgd) < 10^-6


## Explore features of data
lmusc      <- log2(assay(musc) + 0.5)
lcount_ps  <- log2(countdat_ps + 0.5)
lcount_sgd <- log2(countdat_sgd + 0.5)


## Mean gene counts are the same
gene_count_df <- tibble(GTEx       = rowMeans(lmusc),
                        powersimR  = rowMeans(lcount_ps),
                        seqgendiff = rowMeans(lcount_sgd))
gene_count_df %>%
  gather(key = "Dataset", value = "Depth") %>%
  ggplot(aes(x = Dataset, y = Depth)) +
  geom_boxplot() +
  theme_bw() +
  ylab("Mean log Gene Depth") ->
  pl

ggsave(filename = "./output/figures/powsimr_vs_seqgendiff/mean_gene_depth.pdf",
       plot = pl,
       family = "Times",
       width = 4,
       height = 2)

## Mean gene counts are the same
gene_count_df <- tibble(GTEx       = c(lmusc),
                        powersimR  = c(lcount_ps),
                        seqgendiff = c(lcount_sgd))
gene_count_df %>%
  gather(key = "Dataset", value = "Depth") %>%
  ggplot(aes(color = Dataset, x = Depth, lty = Dataset)) +
  geom_freqpoly(lwd = 1) +
  theme_bw() +
  ylab("log(count + 0.5)") +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/figures/powsimr_vs_seqgendiff/count_dist.pdf",
       plot = pl,
       height = 3,
       width = 6,
       family = "Times")

## PCA
sv_ps <- svd(lcount_ps - rowMeans(lcount_ps))
sv_m  <- svd(lmusc - rowMeans(lmusc))
sv_sg <- svd(lcount_sgd - rowMeans(lcount_sgd))


tibble(powersimR = sv_ps$d,
                  GTEx = sv_m$d,
                  seqgendiff = sv_sg$d) %>%
  mutate(svnum = row_number()) ->
  sval_df

sval_df %>%
  gather(-svnum, key = "Dataset", value = "SV") %>%
  ggplot(aes(x = svnum, y = SV, color = Dataset, lty = Dataset)) +
  geom_line(lwd = 1) +
  theme_bw() +
  scale_color_colorblind() +
  xlab("Index") +
  ylab("Singular Value") ->
  pl

ggsave(filename = "./output/figures/powsimr_vs_seqgendiff/scree.pdf",
       plot = pl,
       family = "Times",
       height = 2.6,
       width = 6)

pcdf <- bind_rows(
  tibble(dataset = "powsimR",
         PC1 = sv_ps$v[, 1],
         PC2 = sv_ps$v[, 2],
         group = factor(designmat_ps[, 2])),
  tibble(dataset = "GTEx",
         PC1 = sv_m$v[, 1],
         PC2 = sv_m$v[, 2],
         group = NA),
  tibble(dataset = "seqgendiff",
         PC1 = sv_sg$v[, 1],
         PC2 = sv_sg$v[, 2],
         group = factor(designmat_sgd[, 2])))

ggplot(pcdf, aes(x = PC1, y = PC2, color = group)) +
  facet_wrap(.~dataset) +
  geom_point() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_colorblind(name = "Group", na.translate = TRUE, na.value = ggthemes::colorblind_pal()(3)[3]) ->
  pl

ggsave(filename = "./output/figures/powsimr_vs_seqgendiff/pc_plot.pdf",
       plot = pl,
       family = "Times",
       width = 6,
       height = 2.6)

as_tibble(sv_sg$v)[, 1:10] %>%
  mutate(group = factor(designmat_sgd[, 2])) %>%
  ggplot(aes(x = V1, y = V4, color = group)) +
  geom_point() +
  theme_bw() +
  scale_color_colorblind(name = "Group") +
  xlab("PC1") +
  ylab("PC4") ->
  pl

ggsave(filename = "./output/figures/powsimr_vs_seqgendiff/pc1_pc4.pdf",
       plot = pl,
       family = "Times",
       width = 4,
       height = 2.6)

## Same coefficient distribution
qplot(sort(thout$coefmat), sort(psout$pLFC)) +
  geom_abline() -> pl


## Do voom limma ebayes on each
dge_sgd  <- edgeR::DGEList(counts = countdat_sgd)
dge_sgd  <- edgeR::calcNormFactors(dge_sgd)
vout_sgd <- limma::voom(counts = dge_sgd, design = designmat_sgd, save.plot = TRUE)
lout_sgd <- limma::lmFit(vout_sgd)
eout_sgd <- limma::eBayes(lout_sgd)


dge_ps  <- edgeR::DGEList(counts = countdat_ps)
dge_ps  <- edgeR::calcNormFactors(dge_ps)
vout_ps <- limma::voom(counts = countdat_ps, design = designmat_ps, save.plot = TRUE)
lout_ps <- limma::lmFit(vout_ps)
eout_ps <- limma::eBayes(lout_ps)

plot(beta_ps, eout_ps$coefficients[, 2])
abline(0, 1)
qplot(beta_sgd, eout_sgd$coefficients[, 2]) +
  geom_abline(color = 2, lty = 2, lwd = 1) +
  xlab(TeX("$b_1$")) +
  ylab(TeX("$\\hat{b}_1$")) +
  theme_bw() ->
  pl

ggsave(file = "./output/figures/powsimr_vs_seqgendiff/seqgendiff_truevsfits.pdf",
       plot = pl,
       family = "Times",
       height = 3,
       width = 6)

## Look at mean variance trend of muscle data
dge_m <- edgeR::DGEList(counts = assay(musc))
dge_m <- edgeR::calcNormFactors(dge_m)
vout_m <- limma::voom(counts = dge_m, design = matrix(1, nrow = ncol(musc)), save.plot = TRUE)

## Voom plots
voomdf <- bind_rows(
  tibble(dataset = "GTEx",
         logcount = vout_m$voom.xy$x,
         sqsd = vout_m$voom.xy$y),
  tibble(dataset = "powsimR",
         logcount = vout_ps$voom.xy$x,
         sqsd = vout_ps$voom.xy$y),
  tibble(dataset = "seqgendiff",
         logcount = vout_sgd$voom.xy$x,
         sqsd = vout_sgd$voom.xy$y))
ggplot(voomdf, aes(x = logcount, y = sqsd)) +
  geom_point(size = 0.1, alpha = 1/5) +
  geom_smooth(se = FALSE) +
  facet_grid(dataset ~ .) +
  theme_bw() +
  xlab(TeX("$\\log_2$(count size + 0.5)")) +
  ylab("Square Root Standard Deviation") +
  theme(strip.background = element_rect(fill = "white")) ->
  pl

ggsave("./output/figures/powsimr_vs_seqgendiff/voom_plots.pdf",
       plot = pl,
       family = "Times",
       height = 8,
       width = 6)


