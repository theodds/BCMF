## Packages ----

library(SoftBart)
library(glmnet)
library(truncnorm)
library(rpart)
library(gam)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(distributional)
library(ggdist)
library(mgcv)
library(latex2exp)
library(cowplot)
library(grid)
library(gridExtra)

## A function ----

bart_mediate <- function(data, model_m, model_y, pi_hat, m0_hat, m1_hat,
                         mediator_name, outcome_name, treat_name,
                         n_iter, burnin) {
  X_m <- quantile_normalize_bart(preprocess_df(model.frame(model_m, data = data)
                                               %>% select(-mediator_name))[[1]])
  X_y <- quantile_normalize_bart(preprocess_df(
    cbind(model.frame(model_y, data = data) %>% select(-outcome_name),
          m0_hat, m1_hat))[[1]])
  
  ## Extracting outcomes and treatment and scaling
  
  m <- data[[mediator_name]]
  m_scale <- (m - mean(m)) / sd(m)
  
  Y <- data[[outcome_name]]
  Y_scale <- (Y - mean(Y)) / sd(Y)
  
  A <- data[[treat_name]]
  
  # Hypers for mu_m and tau
  
  hypers_tau             <- Hypers(X_m, m_scale, normalize_Y = FALSE)
  hypers_mu_m            <- hypers_tau
  opts_mu_m              <- Opts(update_s = FALSE)
  opts_tau               <- opts_mu_m
  opts_tau$update_sigma  <- FALSE
  
  # Hypers for mu_y, zeta, d
  hypers_mu_y            <- Hypers(X_y, Y_scale, normalize_Y = FALSE)
  hypers_zeta            <- hypers_mu_y
  hypers_d               <- hypers_mu_y
  opts_mu_y              <- Opts(update_s = FALSE)
  opts_zeta              <- opts_d <- opts_mu_y
  opts_zeta$update_sigma <- FALSE
  opts_d$update_sigma    <- FALSE
  
  # Forests
  forest_mu_m <- MakeForest(hypers_mu_m, opts_mu_m)
  forest_tau  <- MakeForest(hypers_tau, opts_tau)
  forest_mu_y <- MakeForest(hypers_mu_y, opts_mu_y)
  forest_zeta <- MakeForest(hypers_zeta, opts_zeta)
  forest_d    <- MakeForest(hypers_d, opts_d)
  
  # Initialize mu_m and tau
  X_mu_m <- cbind(X_m, pi_hat)
  mu_m <- forest_mu_m$do_predict(X_mu_m)
  tau <- forest_tau$do_predict(X_m)
  
  # Initialize mu_y, zeta, d
  X_mu_y <- cbind(X_y, pi_hat)
  mu_y <- forest_mu_y$do_predict(X_mu_y)
  zeta <- forest_zeta$do_predict(X_y)
  d <- forest_d$do_predict(X_y)
  
  # Matrices for individual/average direct & indirect effects
  zeta_indv <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y))
  delta_indv <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y))
  
  zeta_avg <- rep(NA, n_iter - burnin)
  delta_avg <- rep(NA, n_iter - burnin)
  
  for(i in 1:n_iter) {
    print(i)
    
    mu_m  <- forest_mu_m$do_gibbs(X_mu_m, m_scale - A * tau, X_mu_m, 1)
    tau <- forest_tau$do_gibbs(X_m[A == 1,], m_scale[A == 1] - mu_m[A == 1], X_m, 1)
    
    # Yi(a)
    mu_y  <- forest_mu_y$do_gibbs(X_mu_y, Y_scale - A * zeta - m_scale * d, X_mu_y, 1)
    
    sigma_y <- forest_mu_y$get_sigma()
    forest_zeta$set_sigma(sigma_y)
    forest_d$set_sigma(sigma_y)
    
    zeta  <-
      forest_zeta$do_gibbs(X_y[A == 1,], Y_scale[A == 1] - mu_y[A == 1] -
                             (m_scale * d)[A == 1], X_y, 1)
    R <- (Y_scale - mu_y - A * zeta) / m_scale
    d <- forest_d$do_gibbs_weighted(X_y, R, m_scale^2, X_y, 1)
    
    # Direct and Indirect Effects
    zeta_unscale <- zeta * sd(Y) # direct
    d_unscale <- d * sd(Y) / sd(m)
    tau_unscale <- sd(m) * tau
    # delta <- (pnorm(mu_m + tau) - pnorm(mu_m)) * d_unscale # indirect
    delta <- tau_unscale * d_unscale # indirect
    
    if (i > burnin){
      zeta_indv[i - burnin,] <- zeta_unscale
      delta_indv[i - burnin,] <- delta
      
      omega <- MCMCpack::rdirichlet(1, rep(1, nrow(X_m)))
      zeta_avg[i - burnin] <- sum(zeta_unscale * omega)
      delta_avg[i - burnin] <- sum(delta * omega)
    }
    
  }
  
  return(list(indv_direct = zeta_indv, indv_indirect = delta_indv,
              avg_direct = zeta_avg, avg_indirect = delta_avg))
}


# Model for m and y
model_m <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
model_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth

# Preprocess data
set.seed(123)
meps <- read.csv("data/meps2011.csv") %>%
  filter(totexp > 0 & bmi > 0) %>%
  mutate(logY = log(totexp)) %>%
  mutate(smoke = ifelse(smoke == "No", 0, 1)) %>%
  select(-totexp) %>% select(logY, everything())
phealth <- meps$phealth
phealth <- case_when(phealth == "Poor" ~ 1, phealth == "Fair" ~ 2, 
                     phealth == "Good" ~ 3, phealth == "Very Good" ~ 4, 
                     phealth == "Excellent" ~ 5)
meps$phealth <- phealth
rm(phealth)

# Estimate clever covariates with BART
X_m <- model.frame(model_m, data = meps) %>%
  select(-phealth) %>% preprocess_df() %>% pluck("X")

X_m0 <- X_m[meps$smoke == 0,]
X_m1 <- X_m[meps$smoke == 1,]
m0   <- meps$phealth[meps$smoke == 0]
m1   <- meps$phealth[meps$smoke == 1]

bart_m0 <- softbart(X = X_m0, Y = m0, X_test = X_m)
bart_m1 <- softbart(X = X_m1, Y = m1, X_test = X_m)

m0_hat <- bart_m0$y_hat_test_mean
m1_hat <- bart_m1$y_hat_test_mean

# Estimate of propensity score
model_ps <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
glm_logit <- glm(model_ps, data = meps, family = binomial)
pi_hat <- predict(glm_logit, type = 'response')

out_var_coef <- bart_mediate(meps, model_m, model_y, pi_hat, m0_hat,
                             m1_hat, 'phealth', 'logY', 'smoke', 8000, 4000)

# Traceplots
avg_direct <- rowMeans(out_var_coef$zeta_samples)
avg_indirect <- rowMeans(out_var_coef$tau_samples * out_var_coef$d_samples)
indv_direct <- out_var_coef$zeta_samples
indv_indirect <- out_var_coef$d_samples * out_var_coef$tau_samples
plot(avg_direct)
plot(avg_indirect)

# Histograms
plot_avg_zeta <- ggplot(data.frame(avg_zeta = avg_direct), aes(x = avg_zeta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'cadetblue3') +
  xlab(TeX('$\\bar{\\zeta}$')) + ylab("Frequency") + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot_avg_delta <- ggplot(data.frame(avg_delta = avg_indirect), aes(x = avg_delta)) +
  geom_histogram(bins = 40, color = 'white', fill = 'coral1') +
  xlab(TeX('$\\bar{\\delta}$')) + ylab("") + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))

plot_grid(plot_avg_zeta, plot_avg_delta, align = 'v')

## Waterfall Plots ----
get_quantiles <- function(samples, probs) {
  n <- ncol(samples)
  quantiles <- matrix(unlist(lapply(1:n,
                                    function(i) quantile(samples[,i], probs = probs))), nrow = n, byrow = TRUE)
  return (quantiles)
}

zeta_quantiles <- get_quantiles(indv_direct, probs = c(0.025, 0.5, 0.975))
delta_quantiles <- get_quantiles(indv_indirect, probs = c(0.025, 0.5, 0.975))

df_zeta_q <- data.frame(
  index = 1:nrow(zeta_quantiles),
  start = zeta_quantiles[,1],
  med = zeta_quantiles[,2],
  end = zeta_quantiles[,3]
) %>% arrange(med)

df_delta_q <- data.frame(
  index = 1:nrow(delta_quantiles),
  start = delta_quantiles[,1],
  med = delta_quantiles[,2],
  end = delta_quantiles[,3]
) %>% arrange(med)

ggplot(df_zeta_q) +
  geom_segment(aes(x = start, y = 1:nrow(zeta_quantiles), xend = end, yend = 1:nrow(zeta_quantiles)), color = "dark gray") +
  geom_point(data = df_zeta_q, aes(x = med, y = 1:nrow(zeta_quantiles)), size = 0.2, color = 'red') +
  xlab('Direct Effect (zeta)') +
  ylab('Index') +
  geom_vline(xintercept = 0, size = 1)

ggplot(df_delta_q) +
  geom_segment(aes(x = start, y = 1:nrow(zeta_quantiles), xend = end, yend = 1:nrow(zeta_quantiles)), color = "dark gray") +
  geom_point(data = df_delta_q, aes(x = med, y = 1:nrow(zeta_quantiles)), size = 0.2, color = 'red') +
  xlab('Indirect Effect (delta)') +
  ylab('Index') +
  geom_vline(xintercept = 0, size = 1)

## Posterior projection, summarizing with trees ----

projection_tree <- function(data, model_y, samples) {
  projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
  for (i in 1:nrow(samples)) {
    X <- model.matrix(model_y, data = data)
    formula <- samples[i,] ~ .
    tree_fit <- rpart(formula, data = as.data.frame(X))
    projections[i,] <- predict(tree_fit)
  }
  return (projections)
}

meps_post <- meps %>%
  mutate(delta_hat = colMeans(indv_indirect),
         zeta_hat = colMeans(indv_direct)) %>%
  select(age, bmi, edu, income, povlev, region, sex, marital, race, seatbelt, delta_hat, zeta_hat)

rpart.plot(rpart(delta_hat ~ . - zeta_hat, data = meps_post, 
                 control = rpart.control(cp = 0.03)),
           family = 'Times New Roman', nn.family = 'Times New Roman')

# subgroup1 <- meps_post %>% filter(race != 'White' & age < 34)
# subgroup2 <- meps_post %>% filter(race != 'White' & age >= 34)
# subgroup3 <- meps_post %>% filter(race == 'White' & age < 32)
# subgroup4 <- meps_post %>% filter(race == 'White' & age >= 32 & sex == 'Female')
# subgroup5 <- meps_post %>% filter(race == 'White' & age >= 32 & sex == 'Male')

subgroup1 <- meps_post %>% filter(age < 34 & race == 'White')
subgroup2 <- meps_post %>% filter(age < 34 & race != 'White')
subgroup3 <- meps_post %>% filter(age >= 67)
subgroup4 <- meps_post %>% filter(age >= 34 & age < 67 & race == 'White')
subgroup5 <- meps_post %>% filter(age >= 34 & age < 67 & race != 'White')

plot(NULL, xlim = c(-0.02, 0.13), ylim = c(0, 60), ylab = 'Density', xlab = 'delta')
lines(density(subgroup1$delta_hat), col = 'blue')
lines(density(subgroup2$delta_hat), col = 'blue', lty = 2)
lines(density(subgroup3$delta_hat), col = 'green')
lines(density(subgroup4$delta_hat), col = 'red')
lines(density(subgroup5$delta_hat), col = 'red', lty = 2)
legend('topright',
       legend = c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                  'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67'),
       col = c('blue', 'blue', 'green', 'red', 'red'),
       lty = c(1,2,1,1,2), cex = 0.7)

# Half-eye plot
# subgroup_labels <- c('non-white, age < 34', 'non-white, age ≥ 34', 'white, age < 32',
#                      'white, age ≥ 32, female', 'white, age ≥ 32, male')
subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                     'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
df_subgroups = data.frame(
  subgroup = c(rep(subgroup_labels[1], nrow(subgroup1)),
               rep(subgroup_labels[2], nrow(subgroup2)),
               rep(subgroup_labels[3], nrow(subgroup3)),
               rep(subgroup_labels[4], nrow(subgroup4)),
               rep(subgroup_labels[5], nrow(subgroup5))),
  value = c(subgroup1$delta_hat, subgroup2$delta_hat, subgroup3$delta_hat,
            subgroup4$delta_hat, subgroup5$delta_hat)
)

df_subgroups %>%
  ggplot(aes(x = value, y = subgroup, fill = subgroup)) +
  stat_halfeye(alpha = 1, show.legend = FALSE) + xlab(TeX('$\\delta_a(A_i)$')) +
  ylab("") + scale_fill_manual(values = RColorBrewer::brewer.pal(6, "Blues")[2:6]) +
  theme_bw() + theme(text = element_text(family = "Times New Roman", size = 16)) +
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))


## Posterior projection, summarizing with a GAM ----

projection_gam <- function(data, samples) {
  projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
  for (i in 1:nrow(samples)) {
    formula <- samples[i,] ~ s(age, k = 10) + s(bmi, k = 10) + edu +
               s(income, k = 10) + s(povlev, k = 10) + region + sex + 
               marital + race + seatbelt
    gam_fit <- gam(formula, data = data)
    projections[i,] <- predict(gam_fit, type = 'response')
  }
  return(projections)
}

pdf(NULL)
project_gam <- gam(delta_hat ~ s(age, k = 10) + s(bmi, k = 10) + s(edu, k = 10) +
                     s(income, k = 10) + s(povlev, k = 10) + region + sex + 
                     marital + race + seatbelt, data = meps_post)

smooths <- plot(project_gam)

smooth_age <- list()
smooth_bmi <- list()
smooth_edu <- list()
smooth_income <- list()
smooth_povlev <- list()
region_list <- list()
sex_list <- list()
marital_list <- list()
race_list <- list()
seatbelt_list <- list()
# idx <- floor(seq(from = 1, to = 4000, length = 200))
for(i in 1:4000) {
  print(i)
  project_gam <- gam(delta_hat ~ s(age, k = 10) + s(bmi, k = 10) + s(edu, k = 10) +
                     s(income, k = 10) + s(povlev, k = 10) + region + sex + 
                     marital + race + seatbelt,
                     data = meps_post %>% mutate(delta_hat = indv_indirect[i,]))
  
  # Constraint for binary/categorical variables to sum to 0
  coefs <- coef(project_gam)
  coefs_region0 <- c(regionMidwest = 0, coefs[grep('region', names(coefs))])
  coefs_region <- coefs_region0 - mean(coefs_region0)
  region_list[[i]] <- data.frame(coefs = coefs_region, 
                                 category = names(coefs_region), iteration = i)
  
  coefs_sex0 <- c(sexFemale = 0, coefs[grep('sex', names(coefs))])
  coefs_sex <- coefs_sex0 - mean(coefs_sex0)
  sex_list[[i]] <- data.frame(coefs = coefs_sex,
                              category = names(coefs_sex), iteration = i)
  
  coefs_marital0 <- c(maritalDivorced = 0, coefs[grep('marital', names(coefs))])
  coefs_marital <- coefs_marital0 - mean(coefs_marital0)
  marital_list[[i]] <- data.frame(coefs = coefs_marital,
                                  category = names(coefs_marital), iteration = i)
  
  coefs_race0 <- c(raceWhite = 0, coefs[grep('race', names(coefs))])
  coefs_race <- coefs_race0 - mean(coefs_race0)
  race_list[[i]] <- data.frame(coefs = coefs_race, 
                               category = names(coefs_race), iteration = i)
  
  coefs_seatbelt0 <- c(seatbeltAlmostAlways = 0, coefs[grep('seatbelt', names(coefs))])
  coefs_seatbelt <- coefs_seatbelt0 - mean(coefs_seatbelt0)
  seatbelt_list[[i]] <- data.frame(coefs = coefs_seatbelt, 
                                   category = names(coefs_seatbelt), iteration = i)
  
  p <- plot(project_gam)
  smooth_age[[i]] <- data.frame(x = p[[1]]$x, y = p[[1]]$fit, iteration = i)
  smooth_bmi[[i]] <- data.frame(x = p[[2]]$x, y = p[[2]]$fit, iteration = i)
  smooth_edu[[i]] <- data.frame(x = p[[3]]$x, y = p[[3]]$fit, iteration = i)
  smooth_income[[i]] <- data.frame(x = p[[4]]$x, y = p[[4]]$fit, iteration = i)
  smooth_povlev[[i]] <- data.frame(x = p[[5]]$x, y = p[[5]]$fit, iteration = i)
  
}
dev.off()

smooth_age_df <- do.call(rbind, smooth_age)
smooth_bmi_df <- do.call(rbind, smooth_bmi)
smooth_edu_df <- do.call(rbind, smooth_edu)
smooth_income_df <- do.call(rbind, smooth_income)
smooth_povlev_df <- do.call(rbind, smooth_povlev)
region_df <- do.call(rbind, region_list)
sex_df <- do.call(rbind, sex_list)
marital_df <- do.call(rbind, marital_list)
race_df <- do.call(rbind, race_list)
seatbelt_df <- do.call(rbind, seatbelt_list)

saveRDS(smooth_age_df, 'smooth_age_df.rds')
saveRDS(smooth_bmi_df, 'smooth_bmi_df.rds')
saveRDS(smooth_edu_df, 'smooth_edu_df.rds')
saveRDS(smooth_income_df, 'smooth_income_df.rds')
saveRDS(smooth_povlev_df, 'smooth_povlev_df.rds')
saveRDS(region_df, 'region_df.rds')
saveRDS(sex_df, 'sex_df.rds')
saveRDS(marital_df, 'marital_df.rds')
saveRDS(race_df, 'race_df.rds')
saveRDS(seatbelt_df, 'seatbelt_df.rds')

readRDS('smooth_age_df.rds')
readRDS('smooth_bmi_df.rds')
readRDS('smooth_edu_df.rds')
readRDS('smooth_income_df.rds')
readRDS('smooth_povlev_df.rds')
readRDS('region_df.rds')
readRDS('sex_df.rds')
readRDS('marital_df.rds')
readRDS('race_df.rds')
readRDS('seatbelt_df.rds')

smooth_age_summary <- smooth_age_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_bmi_summary <- smooth_bmi_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_edu_summary <- smooth_edu_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_income_summary <- smooth_income_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))
smooth_povlev_summary <- smooth_povlev_df %>%
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025),
                            UCL = quantile(y, 0.975))

# ggplot(smooth_age_df, aes(x = x, y = y)) +
#   geom_line(aes(group = factor(iteration)), alpha = .1) +
#   stat_summary(geom = "line", color = "chartreuse3",
#                alpha = 3, size = 2)

plot_age <- ggplot(smooth_age_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('age') + ylab('s(age)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_bmi <- ggplot(smooth_bmi_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('bmi') + ylab('s(bmi)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_edu <- ggplot(smooth_edu_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('edu') + ylab('s(edu)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_income <- ggplot(smooth_income_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('income') + ylab('s(income)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_povlev <- ggplot(smooth_povlev_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) +
  geom_line(size = 1, lty = 2) + geom_ribbon(alpha = 0.5, fill = 'cornflowerblue') +
  geom_hline(yintercept = 0) + xlab('poverty level') + ylab('s(poverty level)') + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_grid(plot_age, plot_bmi, plot_edu, plot_income, plot_povlev, align = 'vh')

# Boxplots for binary/categorical coefficients

level_order <- c('seatbeltNoCar', 'seatbeltNever', 'seatbeltAlways', 'seatbeltAlmostAlways', 'seatbeltSometimes', 'seatbeltSeldom',
                 'raceWhite', 'raceBlack', 'raceIndig', 'racePacificIslander', 'raceMulti',
                 'maritalMarried', 'maritalSeparated', 'maritalDivorced', 'maritalWidowed',
                 'sexMale', 'sexFemale',
                 'regionNortheast', 'regionMidwest', 'regionSouth', 'regionWest')
colors <- c('blue', rep('white', 119))
boxplot_coefs <- rbind(cbind(region_df, id = 1), cbind(sex_df, id = 2),
                       cbind(marital_df, id = 3), cbind(race_df, id = 4),
                       cbind(seatbelt_df, id = 5))
ggplot(boxplot_coefs, aes(y = factor(category, level = level_order), x = coefs, fill = as.factor(id))) +
  geom_boxplot(alpha = 0.6, show.legend = FALSE) + scale_fill_brewer(palette =  'Set3') +
  xlab('Coefficients') + ylab('Category') + theme_classic() + 
  theme(text = element_text(family = "Times New Roman")) + scale_colour_identity() +
  geom_hline(yintercept = c(6.5, 11.5, 15.5, 17.5), lty = 2, col = 'darkgray', size = 0.4)


# q <- region_df %>% arrange(factor(category, levels = c('regionWest', 'regionSouth', 'regionMidwest', 'regionNortheast')))
# ggplot(q, aes(y=category, x=coefs)) + geom_boxplot()


## R^2 ----
r_sq <- function(samples, samples_proj) {
  r2 <- rep(NA, nrow(samples))
  for (i in 1:nrow(samples)) {
    num <- sum((samples[i,] - samples_proj[i,])^2)
    denom <- sum((samples[i,] - mean(samples[i,]))^2)
    r2[i] <- 1 - (num / denom)
  }
  
  samples_mean <- colMeans(samples)
  samples_proj_mean <- colMeans(samples_proj)
  num <- sum((samples_mean - samples_proj_mean)^2)
  denom <- sum((samples_mean - mean(samples))^2)
  r2_mean <- 1 - (num / denom)
  
  return(list(r2 = r2, r2_mean = r2_mean))
}

# Tree
proj_tree_delta <- projection_tree(meps, model_y, indv_indirect)
r_sq_tree <- r_sq(indv_indirect, proj_tree_delta)

hist(r_sq_tree$r2)
abline(v = r_sq_tree$r2_mean, col = 'red')

plot_r2_tree <- ggplot(data.frame(r_sq = r_sq_tree$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') + 
  xlab("") + ylab("Frequency") + facet_wrap(~'Tree') +
  geom_vline(xintercept = r_sq_tree$r2_mean) + theme_bw() +
  theme(text = element_text(family = "Times New Roman")) +
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))

# GAM
proj_gam_delta <- projection_gam(meps, indv_indirect)
r_sq_gam <- r_sq(indv_indirect, proj_gam_delta)

hist(r_sq_gam$r2)
abline(v = r_sq_gam$r2_mean, col = 'red')

plot_r2_gam <- ggplot(data.frame(r_sq = r_sq_gam$r2), aes(x = r_sq)) +
  geom_histogram(bins = 30, color = 'white', fill = 'aquamarine3') +
  xlab("") + ylab("") + facet_wrap(~'GAM') +
  geom_vline(xintercept = r_sq_gam$r2_mean) + theme_bw() +
  theme(text = element_text(family = "Times New Roman"))

plot_r2 <- plot_grid(plot_r2_tree, plot_r2_gam, align = 'v')
ggdraw(add_sub(plot_r2, TeX("$R^2$"), fontfamily = "Times New Roman", size = 11,
               vpadding = grid::unit(0, "lines"), y = 6, x = 0.5, vjust = 4.5))
