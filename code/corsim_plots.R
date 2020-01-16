#######################
## Plot output of est_cor.R
#######################

suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
corsimdf <- as_tibble(readRDS("./output/est_cor/corsimout.RDS"))

## Get true correlations
corsimdf %>%
  select(seed, nsamp, target_cor1, target_cor2, starts_with("ce")) %>%
  mutate(corcombo = str_c("(", target_cor1, ", ", target_cor2, ")")) %>%
  gather(contains("ce"), key = "design_est", value = "corest") %>%
  mutate(design_est = str_sub(design_est, 3, 4)) %>%
  separate(col = design_est, into = c("designmat", "column"), sep = 1) %>%
  mutate(column = factor(column),
         designmat = recode(designmat,
                            `1` = "Indicator",
                            `2` = "Normal"),
         target_cor = ifelse(column == 1, target_cor1, target_cor2),
         nsamp = factor(nsamp),
         corcombo = str_c("Target Correlation = ", corcombo)) ->
  dflong

dflong %>%
  distinct(corcombo, designmat, target_cor1, target_cor2) ->
  dfsmall

dflong %>%
  ggplot(aes(x = nsamp, y = corest, color = column)) +
  facet_grid(corcombo ~ designmat) +
  geom_boxplot() +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_colorblind(name = "Column") +
  xlab("Sample Size") +
  ylab("Correlation Estimate") +
  geom_hline(data = dfsmall, aes(yintercept = target_cor1), lty = 2, alpha = 0.5) +
  geom_hline(data = dfsmall, aes(yintercept = target_cor2), lty = 2, alpha = 0.5)->
  pl

ggsave(filename = "./output/figures/cor_est.pdf",
       plot = pl,
       family = "Times",
       height = 7.3,
       width = 6)
