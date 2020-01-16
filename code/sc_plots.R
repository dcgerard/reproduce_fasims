##########################
## Plot the results of the single cell factor analysis sims
##########################

library(tidyverse)
library(ggthemes)
fadf <- read_csv("./output/sc/sc_fa_sims.csv")

## Angle -----------------------------------------------------------------------
fadf %>%
  select(nsamp, ngene, pz, lsd, corval1, starts_with("angle_")) %>%
  gather(starts_with("angle_"), key = "Method", value = "Angle (radians)") %>%
  mutate(Method = str_replace(Method, "angle_", ""),
         pz = recode(pz,
                     "0.0" = "Proportion Zeros = 0",
                     "0.9" = "Proportion Zeros = 0.9"),
         lsd = recode(lsd,
                      "0.8" = "Loadings SD = 0.8",
                      "0.4" = "Loadings SD = 0.4"),
         corval1 = as.factor(corval1)) ->
  angledf

ggplot(angledf, aes(x = corval1, y = `Angle (radians)`, color = Method)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(lsd ~ pz) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Correlation") +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/figures/sc_plots/sc_angle.pdf", 
       plot = pl,
       family = "Times",
       height = 6,
       width = 6)

## MSE of factors --------------------------------------------------------------
fadf %>%
  select(nsamp, ngene, pz, lsd, corval1, starts_with("mse_")) %>%
  gather(starts_with("mse_"), key = "Method", value = "Minimum MSE") %>%
  mutate(Method = str_replace(Method, "mse_", ""),
         pz = recode(pz,
                     "0.0" = "Proportion Zeros = 0",
                     "0.9" = "Proportion Zeros = 0.9"),
         lsd = recode(lsd,
                      "0.8" = "Loadings SD = 0.8",
                      "0.4" = "Loadings SD = 0.4"),
         corval1 = as.factor(corval1)) ->
  msedf

ggplot(msedf, aes(x = corval1, y = `Minimum MSE`, color = Method)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(lsd ~ pz) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Correlation") +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/figures/sc_plots/sc_mse.pdf", 
       plot = pl,
       family = "Times",
       height = 6,
       width = 6)

## MSE of loadings -------------------------------------------------------------
fadf %>%
  select(nsamp, ngene, pz, lsd, corval1, starts_with("loadmse_")) %>%
  gather(starts_with("loadmse_"), key = "Method", value = "Minimum MSE of Loadings") %>%
  mutate(Method = str_replace(Method, "loadmse_", ""),
         pz = recode(pz,
                     "0.0" = "Proportion Zeros = 0",
                     "0.9" = "Proportion Zeros = 0.9"),
         lsd = recode(lsd,
                      "0.8" = "Loadings SD = 0.8",
                      "0.4" = "Loadings SD = 0.4"),
         corval1 = as.factor(corval1)) ->
  loadmsedf

ggplot(loadmsedf, aes(x = corval1, y = `Minimum MSE of Loadings`, color = Method)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_grid(lsd ~ pz) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Correlation") +
  scale_color_colorblind() ->
  pl

ggsave(filename = "./output/figures/sc_plots/sc_loadmse.pdf", 
       plot = pl,
       family = "Times",
       height = 6,
       width = 6)
