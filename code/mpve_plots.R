#####################
## Plot median pve from output of FA simulations in run_sims.R
#####################

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
fadf <- read_csv("./output/fa_sims/fa_results.csv")

fadf %>%
  select(nsamp, ngene, pz, lsd, mpve) %>%
  filter(pz == 0) %>%
  mutate(nsamp = factor(nsamp),
         lsd = as.character(lsd),
         lsd = recode(lsd,
                      "0.4" = "Loading SD = 0.4",
                      "0.7" = "Loading SD = 0.7",
                      "1"   = "Loading SD = 1")) %>%
  ggplot(aes(x = nsamp, y = mpve)) +
  facet_wrap(~lsd) +
  geom_boxplot(outlier.size = 0.5) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  xlab("Sample Size") +
  ylab("Median PVE") ->
  pl

ggsave("./output/figures/fasim_plots/mpve_fasims.pdf",
       plot   = pl,
       height = 2,
       width  = 6,
       family = "Times")
