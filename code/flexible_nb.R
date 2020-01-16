#######################
## Demonstrate flexibility of mixture of binomials and negative binomials
#######################

suppressPackageStartupMessages(library(tidyverse))
mu <- 10
maxlook <- 20
dispvec <- seq(0.001, 0.1, length = 5)
nvec <- round(seq(mu + 1, maxlook, length = 5))
xvec <- seq(0, maxlook, by = 1)
dmat <- matrix(NA, nrow = length(nvec) + length(dispvec), ncol = length(xvec))
for (index in seq_along(nvec)) {
  dmat[index, ] <- stats::dbinom(x = xvec, size = nvec[index], prob = mu / nvec[index])
}
for (index2 in seq_along(dispvec)) {
  dmat[index + index2, ] <- stats::dnbinom(x = xvec, mu = mu, size = 1 / dispvec[index2])
}

ndist <- nrow(dmat)

mixmat <- matrix(NA, nrow = 9, ncol = ndist)
mixmat[1, ] <- dbinom(x = seq_len(ndist) - 1, size = ndist - 1, prob = 0)
mixmat[2, ] <- c(1, rep(0, length = ndist - 2), 1)
mixmat[3, ] <- c(1, 0, 0, 0, 1, 1, 0, 0, 0, 1)
mixmat[4, ] <- dbinom(x = seq_len(ndist) - 1, size = ndist - 1, prob = 1/3)
mixmat[5, ] <- rep(1, length = ndist)
mixmat[6, ] <- dbinom(x = seq_len(ndist) - 1, size = ndist - 1, prob = 1/2)
mixmat[7, ] <- c(0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
mixmat[8, ] <- dbinom(x = seq_len(ndist) - 1, size = ndist - 1, prob = 2/3)
mixmat[9, ] <- dbinom(x = seq_len(ndist) - 1, size = ndist - 1, prob = 1)
mixmat <- mixmat / rowSums(mixmat)


mixdistmat <- mixmat %*% dmat
rownames(mixdistmat) <- seq_len(nrow(mixdistmat))

t(mixdistmat) %>%
  as_tibble() %>%
  mutate(x = xvec) %>%
  gather(-x, key = "distribution", value = "probability") %>%
  ggplot(aes(x = x, xend = x, y = 0, yend = probability)) +
  geom_segment() +
  facet_wrap(.~distribution) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("x") +
  ylab("Pr(x)") ->
  pl

ggsave(filename = "./output/figures/mix_dists.pdf", plot = pl, height = 4, width = 6)

