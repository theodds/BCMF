library(ggplot2)
library(latex2exp)
library(cowplot)

# Average direct and indirect effects
out_bart <- readRDS('data/out_bart.rds')
avg_zeta <- rowMeans(out_bart$zeta_samples)
avg_delta <- rowMeans(out_bart$tau_samples * out_bart$d_samples)

plot_avg_zeta <- ggplot(data.frame(avg_zeta = avg_zeta), aes(x = avg_zeta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'cadetblue3') +
  xlab(TeX('$\\zeta_a$')) + ylab("Frequency") + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot_avg_delta <- ggplot(data.frame(avg_delta = avg_delta), aes(x = avg_delta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'coral1') +
  xlab(TeX('$\\delta_a$')) + ylab("") + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot_grid(plot_avg_zeta, plot_avg_delta, align = 'v')
