library(tidyverse)
library(ggplot2)
library(latex2exp)
library(cowplot)
library(Metrics)
library(extrafont)
font_import()

## MEPS data ----
meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]

subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                     'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
subgroup1 <- which(meps$age < 34 & meps$race == 'White')
subgroup2 <- which(meps$age < 34 & meps$race != 'White')
subgroup3 <- which(meps$age >= 67)
subgroup4 <- which(meps$age >= 34 & meps$age < 67 & meps$race == 'White')
subgroup5 <- which(meps$age >= 34 & meps$age < 67 & meps$race != 'White')

## Simulation Results ----
simulation_bart <- readRDS('data/simulation_bart.rds')
simulation_lsem <- readRDS('data/simulation_lsem.rds')

# True Values
out_lsem <- readRDS('data/out_lsem.rds')
zeta_hat_lsem <- out_lsem$zeta
delta_hat_lsem <- out_lsem$delta

out_bart <- readRDS('data/out_bart.rds')
zeta_hat_bart    <- colMeans(out_bart$zeta_samples)
delta_hat_bart   <- colMeans(out_bart$d_samples) * colMeans(out_bart$tau_samples)
rm(out_bart)
gc()

## Make Truth/Fit data frames ----
bb_avg            <- simulation_bart$average_bart_bart %>%
  mutate(group = 'average', .after = seed) %>%
  mutate(truth = 'BART', fit = 'BART',
         zeta_true = rep(mean(zeta_hat_bart), 200),
         delta_true = rep(mean(delta_hat_bart), 200))
bb_indv           <- simulation_bart$individual_bart_bart %>%
  mutate(truth = 'BART', fit = 'BART',
         zeta_true = rep(zeta_hat_bart[i_test], 200),
         delta_true = rep(delta_hat_bart[i_test], 200))
bb_subgroups      <- simulation_bart$subgroups_bart_bart %>%
  mutate(truth = 'BART', fit = 'BART',
         zeta_true = rep(c(mean(zeta_hat_bart[subgroup1]),
                           mean(zeta_hat_bart[subgroup2]),
                           mean(zeta_hat_bart[subgroup3]),
                           mean(zeta_hat_bart[subgroup4]),
                           mean(zeta_hat_bart[subgroup5])), 200),
         delta_true = rep(c(mean(delta_hat_bart[subgroup1]),
                            mean(delta_hat_bart[subgroup2]),
                            mean(delta_hat_bart[subgroup3]),
                            mean(delta_hat_bart[subgroup4]),
                            mean(delta_hat_bart[subgroup5])), 200))
bb_tree_subgroups <- simulation_bart$tree_subgroups_bart_bart %>%
  mutate(truth = 'BART', fit = 'BART',
         zeta_true = NA,
         delta_true = NA)
lb_avg            <- simulation_bart$average_lsem_bart %>%
  mutate(group = 'average', .after = seed) %>%
  mutate(truth = 'LSEM', fit = 'BART',
         zeta_true = rep(mean(zeta_hat_lsem), 200),
         delta_true = rep(mean(delta_hat_lsem), 200))
lb_indv           <- simulation_bart$individual_lsem_bart %>%
  mutate(truth = 'LSEM', fit = 'BART',
         zeta_true = rep(zeta_hat_lsem[i_test], 200),
         delta_true = rep(delta_hat_lsem[i_test], 200))
ll_avg            <- simulation_lsem$average_lsem_lsem %>%
  mutate(group = 'average', .after = seed) %>%
  mutate(truth = 'LSEM', fit = 'LSEM',
         zeta_true = rep(mean(zeta_hat_lsem), 200),
         delta_true = rep(mean(delta_hat_lsem), 200))
ll_indv           <- simulation_lsem$individual_lsem_lsem %>%
  mutate(truth = 'LSEM', fit = 'LSEM',
         zeta_true = rep(zeta_hat_lsem[i_test], 200),
         delta_true = rep(delta_hat_lsem[i_test], 200))
bl_avg            <- simulation_lsem$average_bart_lsem %>%
  mutate(group = 'average', .after = seed) %>%
  mutate(truth = 'BART', fit = 'LSEM',
         zeta_true = rep(mean(zeta_hat_bart), 200),
         delta_true = rep(mean(delta_hat_bart), 200))
bl_indv           <- simulation_lsem$individual_bart_lsem %>%
  mutate(truth = 'BART', fit = 'LSEM',
         zeta_true = rep(zeta_hat_bart[i_test], 200),
         delta_true = rep(delta_hat_bart[i_test], 200))

indv_df           <- rbind(bb_indv, lb_indv, ll_indv, bl_indv)
avg_subgroups_df  <- rbind(bb_avg, bb_subgroups, bb_tree_subgroups,
                           lb_avg, ll_avg, bl_avg)

# Coverage, RMSE, bias, interval length for average zeta & delta
avg_df_zeta <- avg_subgroups_df %>%
  filter(group == 'average') %>%
  group_by(truth, fit) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len))

avg_df_delta <- avg_subgroups_df %>%
  filter(group == 'average') %>%
  group_by(truth, fit) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            len_delta = mean(delta_len))


# Coverage, RMSE, bias, interval length for individual zeta & delta
indv_df_zeta <- indv_df %>%
  group_by(truth, fit, subj_id) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len))

indv_df_delta <- indv_df %>%
  group_by(truth, fit, subj_id) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            len_delta = mean(delta_len))

# Coverage, RMSE, bias, interval length for subgroup zeta & delta
avg_subgroups_df %>%
  filter(group %in% subgroup_labels) %>%
  group_by(truth, fit, group) %>%
  summarize(cov_zeta = mean(zeta_catch),
            rmse_zeta = rmse(zeta_true, zeta_mean),
            bias_zeta = abs(bias(zeta_true, zeta_mean)),
            len_zeta = mean(zeta_len))

avg_subgroups_df %>%
  filter(group %in% subgroup_labels) %>%
  group_by(truth, fit, group) %>%
  summarize(cov_delta = mean(delta_catch),
            rmse_delta = rmse(delta_true, delta_mean),
            bias_delta = abs(bias(delta_true, delta_mean)),
            len_delta = mean(delta_len))

avg_subgroups_df %>%
  group_by(truth, fit) %>%
  filter(!(group %in% subgroup_labels) & group != 'average') %>%
  summarize(cov_zeta = mean(zeta_catch),
            len_zeta = mean(zeta_len))

avg_subgroups_df %>%
  group_by(truth, fit) %>%
  filter(!(group %in% subgroup_labels) & group != 'average') %>%
  summarize(cov_delta = mean(delta_catch),
            len_delta = mean(delta_len))


## Plots ----
# Zeta
plot_cov_zeta <- ggplot(indv_df_zeta, aes(x = truth, y = cov_zeta, fill = fit)) +
  geom_boxplot(alpha = 0.5) +
  ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman")) +
  geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')
legend <- get_legend(plot_cov_zeta)
plot_cov_zeta <- plot_cov_zeta + theme(legend.position = 'none')

plot_rmse_zeta <- ggplot(indv_df_zeta, aes(x = truth, y = rmse_zeta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('RMSE') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

plot_bias_zeta <- ggplot(indv_df_zeta, aes(x = truth, y = bias_zeta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Bias') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

plot_len_zeta <- ggplot(indv_df_zeta, aes(x = truth, y = len_zeta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Interval Length') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

grid_zeta <- plot_grid(plot_cov_zeta, plot_rmse_zeta, plot_bias_zeta, plot_len_zeta)
plots_zeta <- plot_grid(grid_zeta, legend, ncol = 2, rel_widths = c(1, .1))
ggdraw(add_sub(plots_zeta, TeX('$\\zeta_a(x)$'), fontfamily = "Times New Roman", size = 12,
               vpadding = grid::unit(0, "lines"), y = 5, x = 0.48, vjust = 4.5))

# Delta
plot_cov_delta <- ggplot(indv_df_delta, aes(x = truth, y = cov_delta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Coverage') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman")) +
  geom_hline(yintercept = 0.95, lty = 2, color = 'darkgray')

plot_rmse_delta <- ggplot(indv_df_delta, aes(x = truth, y = rmse_delta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('RMSE') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

plot_bias_delta <- ggplot(indv_df_delta, aes(x = truth, y = bias_delta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Bias') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

plot_len_delta <- ggplot(indv_df_delta, aes(x = truth, y = len_delta, fill = fit)) +
  geom_boxplot(alpha = 0.5, show.legend = FALSE) +
  ylab('Interval Length') + xlab('Truth') + labs(fill = 'Fit') +
  theme_bw() + theme(text = element_text(family = "Times New Roman"))

grid_delta <- plot_grid(plot_cov_delta, plot_rmse_delta, plot_bias_delta, plot_len_delta)
plots_delta <- plot_grid(grid_delta, legend, ncol = 2, rel_widths = c(1, .1))
ggdraw(add_sub(plots_delta, TeX('$\\delta_a(x)$'), fontfamily = "Times New Roman", size = 12,
               vpadding = grid::unit(0, "lines"), y = 5, x = 0.48, vjust = 4.5))
