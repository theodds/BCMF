## Load (120) ----

library(SoftBart)
library(tidyverse)
library(progress)
library(caret)
library(rpart)
library(rpart.plot)
source("R/softbart_bcmf.R")

meps <- read_csv("Data/meps2011.csv") %>%
  filter(totexp > 0 & bmi > 0) %>%
  mutate(logY = log(totexp)) %>%
  mutate(smoke = ifelse(smoke == "No", 0, 1)) %>%
  select(-totexp) %>% select(logY, phealth, smoke, everything())

phealth <- meps$phealth
phealth <- case_when(phealth == "Poor" ~ 1, phealth == "Fair" ~ 2,
                     phealth == "Good" ~ 3, phealth == "Very Good" ~ 4,
                     phealth == "Excellent" ~ 5)
meps$phealth <- phealth
rm(phealth)

## Function for simulating data ----

sim_meps <- function(fit, meps) {
  N <- nrow(meps)
  
  mu_y    <- colMeans(fit$mu_y_train)
  zeta    <- colMeans(fit$zeta_train)
  d       <- colMeans(fit$d_train)
  sigma_y <- mean(fit$sigma_y)
  
  mu_m    <- colMeans(fit$mu_m_train)
  tau_m   <- colMeans(fit$tau_m_train)
  sigma_m <- mean(fit$sigma_m)
  
  phealth <- mu_m + tau_m * meps$smoke + sigma_m * rnorm(N)
  logY    <- mu_y + zeta * meps$smoke + phealth * d + sigma_y * rnorm(N)
  
  new_meps <- meps
  new_meps$phealth <- phealth
  new_meps$logY <- logY
  new_meps$true_delta <- d * tau_m
  new_meps$true_zeta  <- zeta
  
  return(new_meps)
}

## Generate data ----

my_fit <- readRDS("Data/meps_fit.rds")
set.seed(digest::digest2int("generating meps 3"))
fake_meps <- sim_meps(my_fit, meps)
rm(my_fit)
gc()

## Fit model ----

formula_y <- logY ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
formula_m <- phealth ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt

opts <- Opts(num_burn = 2000, num_thin = 1, num_save = 2000, update_s = FALSE, update_sigma_mu = FALSE)

set.seed(digest::digest2int("trying meps 3"))
sim_fit <- softbart_bcmf(formula_y = formula_y, formula_m = formula_m, 
                         trt = fake_meps$smoke, data = fake_meps, 
                         test_data = fake_meps, 
                         opts = opts)

## Checking model fits ----

LCL <- function(x) quantile(x, 0.025)
UCL <- function(x) quantile(x, 0.975)

lcl_delta <- apply(sim_fit$delta_train, 2, LCL)
ucl_delta <- apply(sim_fit$delta_train, 2, UCL)
lcl_zeta  <- apply(sim_fit$zeta_train, 2, LCL)
ucl_zeta  <- apply(sim_fit$zeta_train, 2, UCL)

mean((fake_meps$true_delta >= lcl_delta) & (fake_meps$true_delta <= ucl_delta))
mean((fake_meps$true_zeta >= lcl_zeta) & (fake_meps$true_zeta <= ucl_zeta))

## Looking at subgroups ----

meps_for_subgroup <- fake_meps %>% mutate(delta_hat = colMeans(sim_fit$delta_train))

sim_rpart <- rpart(delta ~ age + bmi + edu + income + povlev + region + sex + 
                     marital + race + seatbelt, data = meps_for_subgroup)

meps_for_subgroup$subgroup <- as.numeric(as.factor(sim_rpart$where))

true_subgroup_deltas <- meps_for_subgroup %>% group_by(subgroup) %>%
  summarise(true_delta = mean(true_delta)) %>% arrange(subgroup)

get_subgroup <- function(i) {
  idx <- which(meps_for_subgroup$subgroup == i)
  delta_samples <- rowMeans(sim_fit$delta_test[,idx])
  lower <- LCL(delta_samples)
  upper <- UCL(delta_samples)
  return(list(delta_hat = mean(delta_samples), LCL = lower, UCL = upper, 
              subgroup = i))
}

get_subgroup_samples <- function(i) {
  idx <- which(meps_for_subgroup$subgroup == i)
  delta_samples <- rowMeans(sim_fit$delta_test[,idx]) - rowMeans(sim_fit$delta_test)
  return(list(delta_samples = delta_samples, subgroup = i, 
              iteration = 1:length(delta_samples)))
}

subgroup_inference <- map_df(1:max(meps_for_subgroup$subgroup), get_subgroup) %>%
  left_join(true_subgroup_deltas) %>%
  mutate(catch = (true_delta >= LCL) & (true_delta <= UCL))

subgroup_samples <- 
  map_df(1:max(meps_for_subgroup$subgroup), get_subgroup_samples)

print(subgroup_inference)
