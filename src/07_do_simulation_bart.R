library(SoftBart)
library(tidyverse)
library(matrixStats)
library(rpart)

source('lib/get_clever_cov.R')
source('lib/get_ps.R')
source('lib/bart_mediate.R')
source('lib/lsem_mediate.R')
source('lib/do_simulation_bart.R')

# MEPS train and test set
meps <- readRDS('data/meps.rds')
n <- nrow(meps)
set.seed(12345)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]

# True effects with BART
formula_m_bart <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_y_bart <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth
formula_ps_bart <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

out_bart <- readRDS('data/out_bart.rds')
mu_y_hat_bart    <- colMeans(out_bart$mu_y_samples)
zeta_hat_bart    <- colMeans(out_bart$zeta_samples)
d_hat_bart       <- colMeans(out_bart$d_samples)
mu_m_hat_bart    <- colMeans(out_bart$mu_m_samples)
tau_hat_bart     <- colMeans(out_bart$tau_samples)
sigma_y_hat_bart <- mean(out_bart$sigma_y_samples)
sigma_m_hat_bart <- mean(out_bart$sigma_m_samples)
rm(out_bart)
gc()

# True effects with LSEM
formula_m_lsem <- phealth ~ smoke * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_y_lsem <- logY ~ (smoke + phealth) * (age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt)
formula_X_lsem <- phealth ~ age + bmi + edu + log(income + 1000) + povlev + region + sex + marital + race + seatbelt

fit_m_lsem <- lm(formula_m_lsem, data = meps)
fit_y_lsem <- lm(formula_y_lsem, data = meps)

out_lsem <- lsem_mediate(fit_m_lsem, fit_y_lsem, formula_X_lsem, meps, 'phealth', 'smoke')
saveRDS(out_lsem, 'data/out_lsem.rds')
zeta_hat_lsem <- out_lsem$zeta
delta_hat_lsem <- out_lsem$delta
sigma_m_hat_lsem <- out_lsem$sigma_m
sigma_y_hat_lsem <- out_lsem$sigma_y
rm(out_lsem)
gc()

set.seed(1)
seeds <- sample.int(10e6, 200)
simulation_bart <- do_simulation_bart(meps, i_train, i_test, formula_m_bart,
                                      formula_y_bart, formula_ps_bart,
                                      fit_m_lsem, fit_y_lsem, formula_X_lsem,
                                      'phealth', 'logY', 'smoke',
                                      mu_y_hat_bart, zeta_hat_bart, d_hat_bart,
                                      mu_m_hat_bart, tau_hat_bart,
                                      sigma_y_hat_bart, sigma_m_hat_bart,
                                      zeta_hat_lsem, delta_hat_lsem,
                                      sigma_y_hat_lsem, sigma_m_hat_lsem,
                                      8000, 4000, 200, seeds)
