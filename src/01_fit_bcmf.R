library(SoftBart)
library(tidyverse)
library(ggplot2)
library(cowplot)

source('lib/bart_mediate2.R')
source('lib/get_clever_cov.R')
source('lib/get_ps.R')

# MEPS data
set.seed(123)
meps <- readRDS('data/meps.rds')

# Define formulas
formula_m <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

# Get clever covariates and propensity score
clever_cov <- get_clever_cov(meps, meps, formula_m,
                             'phealth', 'logY', 'smoke')
m0_hat <- clever_cov$m0_hat
m1_hat <- clever_cov$m1_hat

pi_hat <- get_ps(meps, meps, formula_ps)

# Fit BCMF
out_bart <- bart_mediate2(meps, meps, formula_m, formula_y, pi_hat, pi_hat,
                         m0_hat, m0_hat, m1_hat, m1_hat,
                         'phealth', 'logY', 'smoke',
                         5000, 4000, 4)
saveRDS(out_bart, 'data/out_bart.rds')

# Trace Plots
chain <- as.factor(c(rep(1, 1000), rep(2, 1000), rep(3, 1000), rep(4, 1000)))
indv_zeta <- out_bart$zeta_samples
indv_delta <- out_bart$tau_samples * out_bart$d_samples
avg_zeta <- data.frame(i = rep(1:1000, times = 4), avg_zeta = rowMeans(indv_zeta), chain = chain)
avg_delta <- data.frame(i = rep(1:1000, times = 4), avg_delta = rowMeans(indv_delta), chain = chain)
sigma_m <- data.frame(i = rep(1:1000, times = 4), sigma_m = out_bart$sigma_m_samples, chain = chain)
sigma_y <- data.frame(i = rep(1:1000, times = 4), sigma_y = out_bart$sigma_y_samples, chain = chain)

traceplot_avg_zeta <- ggplot(avg_zeta, aes(x = i, y = avg_zeta, group = chain)) +
  geom_line(aes(color = chain)) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'zeta_avg') +
  theme(text = element_text(family = "Times New Roman"))
legend <- get_legend(traceplot_avg_zeta)
traceplot_avg_zeta <- traceplot_avg_zeta + theme(legend.position = 'none')

traceplot_avg_delta <- ggplot(avg_delta, aes(x = i, y = avg_delta, group = chain)) +
  geom_line(aes(color = chain), show.legend = FALSE) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'delta_avg') +
  theme(text = element_text(family = "Times New Roman"))

traceplot_sigma_m <- ggplot(sigma_m, aes(x = i, y = sigma_m, group = chain)) +
  geom_line(aes(color = chain), show.legend = FALSE) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'sigma_m') +
  theme(text = element_text(family = "Times New Roman"))

traceplot_sigma_y <- ggplot(sigma_y, aes(x = i, y = sigma_y, group = chain)) +
  geom_line(aes(color = chain), show.legend = FALSE) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'sigma_y') +
  theme(text = element_text(family = "Times New Roman"))

grid <- plot_grid(traceplot_avg_zeta, traceplot_avg_delta, traceplot_sigma_m, traceplot_sigma_y, align = 'vh')
plot_grid(grid, legend, ncol = 2, rel_widths = c(1, .1))


i <- sample.int(nrow(meps), 100)
zeta_subset <- data.frame(subj_id = as.factor(rep(i, each = 4000)),
                          i = rep(1:1000, times = 4 * length(i)),
                          zeta = as.vector(indv_zeta[,i]),
                          chain = rep(chain, length(i)))
delta_subset <- data.frame(subj_id = as.factor(rep(i, each = 4000)),
                           i = rep(1:1000, times = 4 * length(i)),
                           delta = as.vector(indv_delta[,i]),
                           chain = rep(chain, length(i)))

traceplot_zeta_subset <- ggplot(zeta_subset, aes(x = i, y = zeta, group = chain)) +
  geom_line(aes(color = chain)) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'zeta_indv') +
  theme(text = element_text(family = "Times New Roman"))
legend_subset <- get_legend(traceplot_zeta_subset)
traceplot_zeta_subset <- traceplot_zeta_subset + theme(legend.position = 'none')

traceplot_delta_subset <- ggplot(delta_subset, aes(x = i, y = delta, group = chain)) +
  geom_line(aes(color = chain), show.legend = FALSE) +
  theme_bw() + xlab('') + ylab('') + facet_wrap(~'delta_indv') +
  theme(text = element_text(family = "Times New Roman"))

grid_subset <- plot_grid(traceplot_zeta_subset, traceplot_delta_subset, align = 'vh')
plot_grid(grid_subset, legend_subset, ncol = 2, rel_widths = c(1, .1))
