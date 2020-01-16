##############################
## Plot output of run_sims.R
##############################

suppressPackageStartupMessages(library(tidyverse))
library(latex2exp)
library(ggthemes)
fadf <- read_csv("./output/fa_sims/fa_results.csv")

## Angle metric ---------------------------------------------------------------
fadf %>%
  select(starts_with("angle"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(starts_with("angle"), key = "method", value = "angle") %>%
  mutate(method = str_replace(method, "angle_", "")) %>%
  nest(-corval1, -corval2) ->
  sepangledf

for (index in seq_len(nrow(sepangledf))) {
  sepangledf$data[[index]] %>%
    mutate(pz = str_c("Proportion Zeros = ", pz),
           lsd = str_c("Loadings SD = ", lsd),
           nsamp = factor(nsamp),
           method = factor(method, levels = c("ssvd", "pca", "flash", "ica", "peer"))) %>%
    ggplot(aes(x = nsamp, y = angle, color = method)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_grid(pz ~ lsd) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Sample Size") +
    ylab("Angle (radians)") +
    scale_color_colorblind(name = "Method") +
    scale_y_log10() ->
    pl

  ggsave(filename = str_c("./output/figures/fasim_plots/angle_",
                          sepangledf$corval1[index] * 10,
                          "_",
                          sepangledf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}

## MSE metric -----------------------------------------------------------------
fadf %>%
  select(starts_with("mse"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(starts_with("mse"), key = "method", value = "mse") %>%
  mutate(method = str_replace(method, "mse_", "")) %>%
  nest(-corval1, -corval2) ->
  sepmsedf

for (index in seq_len(nrow(sepmsedf))) {
  sepmsedf$data[[index]] %>%
    mutate(pz = str_c("Proportion Zeros = ", pz),
           lsd = str_c("Loadings SD = ", lsd),
           nsamp = factor(nsamp),
           method = factor(method, levels = c("ssvd", "pca", "flash", "ica", "peer"))) %>%
    ggplot(aes(x = nsamp, y = mse, color = method)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_grid(pz ~ lsd) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Sample Size") +
    ylab("Minimum MSE") +
    scale_color_colorblind(name = "Method") ->
    pl

  ggsave(filename = str_c("./output/figures/fasim_plots/minmse_",
                          sepmsedf$corval1[index] * 10,
                          "_",
                          sepmsedf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}


## Loadings MSE metric --------------------------------------------------------
fadf %>%
  select(starts_with("loadmse"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(starts_with("loadmse"), key = "method", value = "loadmse") %>%
  mutate(method = str_replace(method, "loadmse_", "")) %>%
  nest(-corval1, -corval2) ->
  seploadmsedf

for (index in seq_len(nrow(seploadmsedf))) {
  seploadmsedf$data[[index]] %>%
    mutate(pz = str_c("Proportion Zeros = ", pz),
           lsd = str_c("Loadings SD = ", lsd),
           nsamp = factor(nsamp),
           method = factor(method, levels = c("ssvd", "pca", "flash", "ica", "peer"))) %>%
    ggplot(aes(x = nsamp, y = loadmse, color = method)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_grid(pz ~ lsd) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Sample Size") +
    ylab("Minimum MSE of Loadings") +
    scale_color_colorblind(name = "Method") ->
    pl

  ggsave(filename = str_c("./output/figures/fasim_plots/minloadmse_",
                          seploadmsedf$corval1[index] * 10,
                          "_",
                          seploadmsedf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}
