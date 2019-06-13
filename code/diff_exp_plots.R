##################
## Plots from differential expression analysis simulations
##################
suppressPackageStartupMessages(library(tidyverse))
library(ggthemes)
powsimr_df <- as_tibble(readRDS(file = "./output/diff_exp_out/powsimr_sims.RDS"))
seqgendiff_df <- as_tibble(readRDS(file = "./output/diff_exp_out/seqgendiff_sims.RDS"))

## Compare median PVE on non-null genes ---------------------------------------
powsimr_df %>%
  select(lfc_sd, mpve) %>%
  mutate(method = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(mpve) %>%
  mutate(method = "seqgendiff", lfc_sd = 0.8) %>%
  select(lfc_sd, mpve, method) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(method_lfc = paste0(method, "\nsd(lfc) = ", lfc_sd)) %>%
  ggplot(aes(x = method_lfc, y = mpve)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Simulation Method") +
  ylab("Median PVE") ->
  pl

ggsave(filename = "./output/figures/diff_exp/mpve.pdf",
       family = "Times",
       height = 2.6,
       width = 3)

## Compare fdp ----------------------------------------------------------------
powsimr_df %>%
  select(contains("fpr"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "fdp") %>%
  mutate(method = str_replace(method, "fpr_", ""),
         sim = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(contains("fpr")) %>%
  gather(key = "method", value = "fdp") %>%
  mutate(method = str_replace(method, "fpr_", ""),
         sim = "seqgendiff",
         lfc_sd = 0.8) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = paste0(sim, "\nsd(lfc) = ", lfc_sd)) %>%
  mutate(method = recode(method,
                         "dout" = "DESeq2",
                         "eout" = "edgeR",
                         "vout" = "voom+limma")) %>%
  ggplot(aes(x = method, y = fdp, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  geom_hline(yintercept = 0.05, lty = 2) +
  xlab("Method") +
  ylab("FDP") ->
  pl

ggsave(filename = "./output/figures/diff_exp/fdp.pdf",
       family = "Times",
       height = 2.6,
       width  = 4)

## Compare power ----------------------------------------------------------------
powsimr_df %>%
  select(contains("power"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "power") %>%
  mutate(method = str_replace(method, "power_", ""),
         sim = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(contains("power")) %>%
  gather(key = "method", value = "power") %>%
  mutate(method = str_replace(method, "power_", ""),
         sim = "seqgendiff",
         lfc_sd = 0.8) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = paste0(sim, "\nsd(lfc) = ", lfc_sd)) %>%
  mutate(method = recode(method,
                         "dout" = "DESeq2",
                         "eout" = "edgeR",
                         "vout" = "voom+limma")) %>%
  ggplot(aes(x = method, y = power, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  xlab("Method") +
  ylab("Power") ->
  pl

ggsave(filename = "./output/figures/diff_exp/power.pdf",
       family = "Times",
       height = 2.6,
       width  = 4)

## Compare mse ----------------------------------------------------------------
powsimr_df %>%
  select(contains("mse"), lfc_sd) %>%
  gather(-lfc_sd, key = "method", value = "mse") %>%
  mutate(method = str_replace(method, "mse_", ""),
         sim = "powsimR") ->
  pstemp

seqgendiff_df %>%
  select(contains("mse")) %>%
  gather(key = "method", value = "mse") %>%
  mutate(method = str_replace(method, "mse_", ""),
         sim = "seqgendiff",
         lfc_sd = 0.8) ->
  sgdtemp

bind_rows(pstemp, sgdtemp) %>%
  mutate(sim_lfc = paste0(sim, "\nsd(lfc) = ", lfc_sd)) %>%
  mutate(method = recode(method,
                         "dout" = "DESeq2",
                         "eout" = "edgeR",
                         "vout" = "voom+limma")) %>%
  ggplot(aes(x = method, y = mse, color = sim_lfc)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_color_colorblind(name = "Simulation\nMethod") +
  theme_bw() +
  xlab("Method") +
  ylab("MSE") ->
  pl

ggsave(filename = "./output/figures/diff_exp/mse.pdf",
       family = "Times",
       height = 2.6,
       width  = 4)













