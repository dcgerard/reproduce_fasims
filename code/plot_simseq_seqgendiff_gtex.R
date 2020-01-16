############################
## Plot the result of running SimSeq on the GTEx data
############################

library(tidyverse)
library(ggthemes)
library(gridExtra)

seqgendiff_small <- as_tibble(readRDS("./output/diff_exp_out/seqgendiff_small_mpve_sims.RDS"))
seqgendiff_small$type <- "seqgendiff"
simseq <- as_tibble(readRDS("./output/diff_exp_out/simseq_sims.RDS"))
simseq$type <- "SimSeq"

df <- bind_rows(seqgendiff_small, simseq)
outsize <- 0.1

df %>%
  select(type, starts_with("fpr_")) %>%
  gather(-type, key = "method", value = "fpr") %>%
  mutate(method = recode(method,
                         fpr_dout = "DESeq2",
                         fpr_eout = "edgeR",
                         fpr_vout = "voom-limma")) %>%
  ggplot(aes(x = method, y = fpr, color = type)) +
  geom_boxplot(outlier.size = outsize) +
  theme_bw() +
  scale_color_colorblind(name = "Simulation\nType") +
  theme(axis.text = element_text(angle = 90, hjust = 1)) +
  xlab("Method") +
  ylab("FPR") + 
  ggtitle("(A)") +
  guides(color = FALSE) ->
  pl1

df %>%
  select(type, starts_with("power_")) %>%
  gather(-type, key = "method", value = "power") %>%
  mutate(method = recode(method,
                         power_dout = "DESeq2",
                         power_eout = "edgeR",
                         power_vout = "voom-limma")) %>%
  ggplot(aes(x = method, y = power, color = type)) +
  geom_boxplot(outlier.size = outsize) +
  theme_bw() +
  scale_color_colorblind(name = "Simulation\nType") +
  theme(axis.text = element_text(angle = 90, hjust = 1)) +
  xlab("Method") +
  ylab("Power") +
  ggtitle("(B)") ->
  pl2

pdf(file = "./output/figures/simseq_plots/simseq_gtex.pdf",
    width = 6, 
    height = 3, 
    family = "Times")
gridExtra::grid.arrange(pl1, pl2, nrow = 1, widths = c(2.4, 3.6))
dev.off()


df %>%
  select(type, time) %>%
  ggplot(aes(x = type, y = time)) +
  geom_boxplot(outlier.size = 0.3) +
  theme_bw() +
  xlab("Simulation Method") +
  ylab("Time (seconds)") ->
  pl

ggsave(filename = "./output/figures/simseq_plots/simseq_time.pdf", plot = pl, family = "Times", height = 3, width = 3)

df %>%
  group_by(type) %>%
  summarize(avetime = mean(time))
