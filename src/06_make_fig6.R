library(ggplot2)
library(rpart)
library(mgcv)
library(cowplot)
library(latex2exp)

source('lib/projection_tree.R')
source('lib/projection_gam.R')
source('lib/r_sq.R')

out_bart <- readRDS('data/out_bart.rds')
indv_indirect <- out_bart$d_samples * out_bart$tau_samples

meps <- readRDS('data/meps.rds')
formula_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth

# Tree
proj_tree_delta <- projection_tree(meps, formula_y, indv_indirect)
r_sq_tree <- r_sq(indv_indirect, proj_tree_delta)

plot_r2_tree <- ggplot(data.frame(r_sq = r_sq_tree$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') + 
  xlab("") + ylab("Frequency") + facet_wrap(~'Tree') +
  geom_vline(xintercept = r_sq_tree$r2_mean) + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# GAM
proj_gam_delta <- projection_gam(meps, indv_indirect)
r_sq_gam <- r_sq(indv_indirect, proj_gam_delta)

plot_r2_gam <- ggplot(data.frame(r_sq = r_sq_gam$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') +
  xlab("") + ylab("") + facet_wrap(~'GAM') +
  geom_vline(xintercept = r_sq_gam$r2_mean) + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_r2 <- plot_grid(plot_r2_tree, plot_r2_gam, align = 'v')
ggdraw(add_sub(plot_r2, TeX("$R^2$"), fontfamily = "Times New Roman", size = 11,
               vpadding = grid::unit(0, "lines"), y = 6, x = 0.54, vjust = 4.5))
