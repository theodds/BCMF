library(SoftBart)
library(tidyverse)

source('lib/bart_mediate.R')
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
out_bart <- bart_mediate(meps, meps, formula_m, formula_y, pi_hat, pi_hat,
                         m0_hat, m0_hat, m1_hat, m1_hat,
                         'phealth', 'logY', 'smoke',
                         8000, 4000)
saveRDS(out_bart, 'data/out_bart.rds')
