##############################
## Plot output of run_sims.R
##############################

suppressPackageStartupMessages(library(tidyverse))
library(latex2exp)
library(ggthemes)
fadf <- read_csv("./output/fa_sims/fa_results.csv")

## Correlation difference metric ----------------------------------------------
fadf %>%
  select(contains("cordiff"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(contains("cordiff"), key = "method", value = "cordiff") %>%
  mutate(method = str_replace(method, "cordiff_", "")) %>%
  nest(-corval1, -corval2) ->
  sepcordf

for (index in seq_len(nrow(sepcordf))) {
  sepcordf$data[[index]] %>%
    mutate(pz = str_c("Proportion Zeros = ", pz),
           lsd = str_c("Loadings SD = ", lsd),
           nsamp = factor(nsamp),
           method = factor(method, levels = c("ssvd", "pca", "flash", "ica", "peer"))) %>%
    ggplot(aes(x = nsamp, y = cordiff, color = method)) +
    geom_boxplot(outlier.size = 0.1) +
    facet_grid(pz ~ lsd) +
    theme_bw() +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Sample Size") +
    ylab("Correlation Difference") +
    scale_color_colorblind(name = "Method") ->
    pl

  ggsave(filename = str_c("./output/figures/fasim_plots/cordiff_",
                          sepcordf$corval1[index] * 10,
                          "_",
                          sepcordf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}

## Angle metric ---------------------------------------------------------------
fadf %>%
  select(contains("angle"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(contains("angle"), key = "method", value = "angle") %>%
  mutate(method = str_replace(method, "angle_", "")) %>%
  nest(-corval1, -corval2) ->
  sepangledf

for (index in seq_len(nrow(sepcordf))) {
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
                          sepcordf$corval1[index] * 10,
                          "_",
                          sepcordf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}

## MSE metric -----------------------------------------------------------------
fadf %>%
  select(contains("mse"), nsamp, pz, lsd, corval1, corval2) %>%
  gather(contains("mse"), key = "method", value = "mse") %>%
  mutate(method = str_replace(method, "mse_", "")) %>%
  nest(-corval1, -corval2) ->
  sepmsedf

for (index in seq_len(nrow(sepcordf))) {
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
                          sepcordf$corval1[index] * 10,
                          "_",
                          sepcordf$corval2[index] * 10,
                          ".pdf"),
         plot = pl,
         family = "Times",
         height = 6,
         width = 6)
}






