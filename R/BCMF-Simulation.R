## Packages ----

library(SoftBart)
library(glmnet)
library(truncnorm)
library(rpart)
library(gam)
library(tidyverse)
library(rpart)
library(rpart.plot)
library(mgcv)
library(matrixStats)

bart_mediate_sim <- function(data_train, data_test, model_m, model_y,
                             pi_hat_train, pi_hat_test,
                             m0_hat_train, m0_hat_test,
                             m1_hat_train, m1_hat_test,
                             mediator_name, outcome_name, treat_name,
                             n_iter, burnin) {

  X_m_train <- quantile_normalize_bart(preprocess_df(
    model.frame(model_m, data = data_train) %>% select(-mediator_name))[[1]])
  X_y_train <- quantile_normalize_bart(preprocess_df(
    cbind(model.frame(model_y, data = data_train) %>% select(-outcome_name),
          m0_hat_train, m1_hat_train))[[1]])
  X_mu_m_train <- cbind(X_m_train, pi_hat_train)
  X_mu_y_train <- cbind(X_y_train, pi_hat_train)

  X_m_test <- quantile_normalize_bart(preprocess_df(
    model.frame(model_m, data = data_test) %>% select(-mediator_name))[[1]])
  X_y_test <- quantile_normalize_bart(preprocess_df(
    cbind(model.frame(model_y, data = data_test) %>% select(-outcome_name),
          m0_hat_test, m1_hat_test))[[1]])
  X_mu_m_test <- cbind(X_m_test, pi_hat_test)
  X_mu_y_test <- cbind(X_y_test, pi_hat_test)

  # Extracting mediator, outcomes, treatment and scaling

  m <- data_train[[mediator_name]]
  m_scale <- (m - mean(m)) / sd(m)

  Y <- data_train[[outcome_name]]
  Y_scale <- (Y - mean(Y)) / sd(Y)

  A <- data_train[[treat_name]]

  # Hypers for mu_m and tau

  hypers_tau             <- Hypers(X_m_train, m_scale, normalize_Y = FALSE)
  hypers_mu_m            <- hypers_tau
  opts_mu_m              <- Opts(update_s = FALSE)
  opts_tau               <- opts_mu_m
  opts_tau$update_sigma  <- FALSE

  # Hypers for mu_y, zeta, d
  hypers_mu_y            <- Hypers(X_y_train, Y_scale, normalize_Y = FALSE)
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
  mu_m <- forest_mu_m$do_predict(X_mu_m_train)
  tau <- forest_tau$do_predict(X_m_train)

  # Initialize mu_y, zeta, d
  mu_y <- forest_mu_y$do_predict(X_mu_y_train)
  zeta <- forest_zeta$do_predict(X_y_train)
  d <- forest_d$do_predict(X_y_train)

  # Store samples
  mu_y_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
  zeta_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
  d_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
  mu_m_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
  tau_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
  sigma_y_samples <- rep(NA, n_iter - burnin)
  sigma_m_samples <- rep(NA, n_iter - burnin)

  for(i in 1:n_iter) {
    print(i)

    # Mi(a)
    mu_m  <- forest_mu_m$do_gibbs(X_mu_m_train, m_scale - A * tau, X_mu_m_train, 1)
    mu_m_test <- forest_mu_m$do_predict(X_mu_m_test)
    tau <- forest_tau$do_gibbs(X_m_train[A == 1,], m_scale[A == 1] - mu_m[A == 1], X_m_train, 1)
    tau_test <- forest_tau$do_predict(X_m_test)

    sigma_m <- forest_mu_m$get_sigma()
    forest_tau$set_sigma(sigma_m)

    # Yi(a)
    mu_y  <- forest_mu_y$do_gibbs(X_mu_y_train, Y_scale - A * zeta - m_scale * d, X_mu_y_train, 1)
    mu_y_test <- forest_mu_y$do_predict(X_mu_y_test)

    sigma_y <- forest_mu_y$get_sigma()
    forest_zeta$set_sigma(sigma_y)
    forest_d$set_sigma(sigma_y)

    zeta  <-
      forest_zeta$do_gibbs(X_y_train[A == 1,], Y_scale[A == 1] - mu_y[A == 1] -
                                          (m_scale * d)[A == 1], X_y_train, 1)
    zeta_test <- forest_zeta$do_predict(X_y_test)
    R <- (Y_scale - mu_y - A * zeta) / m_scale
    d <- forest_d$do_gibbs_weighted(X_y_train, R, m_scale^2, X_y_train, 1)
    d_test <- forest_d$do_predict(X_y_test)

    # Unscale
    mu_y_unscale <- mean(Y) + (mu_y_test * sd(Y)) - mean(m) / sd(m) * d_test * sd(Y)
    zeta_unscale <- zeta_test * sd(Y)
    d_unscale <- d_test * sd(Y) / sd(m)
    mu_m_unscale <- mean(m) + (mu_m_test * sd(m))
    tau_unscale <- tau_test * sd(m)
    sigma_y_unscale <- sigma_y * sd(Y)
    sigma_m_unscale <- sigma_m * sd(m)

    if (i > burnin){
      mu_y_samples[i - burnin,] <- mu_y_unscale
      zeta_samples[i - burnin,] <- zeta_unscale
      d_samples[i - burnin,] <- d_unscale
      mu_m_samples[i - burnin,] <- mu_m_unscale
      tau_samples[i - burnin,] <- tau_unscale
      sigma_y_samples[i - burnin] <- sigma_y_unscale
      sigma_m_samples[i - burnin] <- sigma_m_unscale
    }

  }

  return(list(mu_y_samples = mu_y_samples, zeta_samples = zeta_samples,
              d_samples = d_samples, mu_m_samples = mu_m_samples,
              tau_samples = tau_samples, sigma_y_samples = sigma_y_samples,
              sigma_m_samples = sigma_m_samples))
}


get_clever_cov <- function(data_train, data_test, model_m,
                           mediator_name, outcome_name, treat_name) {
  X_m_test  <- model.frame(model_m, data = data_test) %>%
               select(-mediator_name) %>% preprocess_df() %>% pluck("X")
  X_m_train  <- model.frame(model_m, data = data_train) %>%
                select(-mediator_name) %>% preprocess_df() %>% pluck("X")
  X_m0_train <- X_m_train[data_train[[treat_name]] == 0,]
  X_m1_train <- X_m_train[data_train[[treat_name]] == 1,]
  m0_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 0]
  m1_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 1]

  bart_m0 <- softbart(X_m0_train, m0_train, X_m_test)
  bart_m1 <- softbart(X_m1_train, m1_train, X_m_test)

  m0_hat <- bart_m0$y_hat_test_mean
  m1_hat <- bart_m1$y_hat_test_mean

  return(list(m0_hat = m0_hat, m1_hat = m1_hat))
}

get_ps <- function(data_train, data_test, model_ps) {
  glm_logit <- glm(model_ps, data = data_train, family = binomial)
  pi_hat <- predict(glm_logit, newdata = data_test,type = 'response')
  return (pi_hat)
}


# MEPS data
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

# Model for m and y
model_m <- phealth ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
model_y <- logY ~ -1 + age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt + phealth

# Clever covariates for entire set
clever_cov <- get_clever_cov(meps, meps, model_m, 'phealth', 'logY', 'smoke')
m0_hat <- clever_cov$m0_hat
m1_hat <- clever_cov$m1_hat

# Propensity score for entire set
model_ps <- smoke ~ age + bmi + edu + income + povlev + region + sex + marital + race + seatbelt
pi_hat <- get_ps(meps, meps, model_ps)

# Get true values
out <- bart_mediate_sim(meps, meps, model_m, model_y,
                        pi_hat, pi_hat,
                        m0_hat, m0_hat,
                        m1_hat, m1_hat,
                        'phealth', 'logY', 'smoke',
                        8000, 4000)

mu_y_hat    <- colMeans(out$mu_y_samples)
zeta_hat    <- colMeans(out$zeta_samples)
d_hat       <- colMeans(out$d_samples)
mu_m_hat    <- colMeans(out$mu_m_samples)
tau_hat     <- colMeans(out$tau_samples)
sigma_y_hat <- mean(out$sigma_y_samples)
sigma_m_hat <- mean(out$sigma_m_samples)

# Get training & testing set indices
n <- nrow(meps)
i_train <- sample(1:nrow(meps), floor(n/2))
i_test <- c(1:nrow(meps))[-i_train]


do_simulation <- function(data, i_train, i_test, model_m, model_y, model_ps,
                          mediator_name, outcome_name, treat_name,
                          mu_y_hat, zeta_hat, d_hat, 
                          mu_m_hat, tau_hat, sigma_y_hat, sigma_m_hat,
                          n_iter, burnin, n_reps, seeds) {
  
  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]
  
  # Store results
  colnames_indv      <- c('seed', 'subj_id', 'zeta_mean', 'zeta_lower',
                          'zeta_upper', 'zeta_catch', 'delta_mean', 'delta_lower',
                          'delta_upper', 'delta_catch')
  colnames_subgroup <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
                          'zeta_upper', 'zeta_catch', 'delta_mean', 'delta_lower',
                          'delta_upper', 'delta_catch')
  colnames_avg       <- c('seed', 'zeta_mean', 'zeta_lower',
                          'zeta_upper', 'zeta_catch', 'delta_mean', 'delta_lower',
                          'delta_upper', 'delta_catch')
  colnames_tree_subgroups <- c('seed', 'group', 'zeta_mean', 'zeta_lower',
                               'zeta_upper', 'delta_mean', 'delta_lower',
                               'delta_upper', 'zeta_catch', 'delta_catch')
  
  indv_mat      <- matrix(NA, nrow = n_reps * nrow(data_test),
                          ncol = length(colnames_indv),
                          dimnames = list(c(), colnames_indv))
  subgroup_mat <- matrix(NA, nrow = n_reps * 5,
                         ncol = length(colnames_subgroup),
                         dimnames = list(c(), colnames_subgroup))
  avg_mat       <- matrix(NA, nrow = n_reps,
                          ncol = length(colnames_avg),
                          dimnames = list(c(), colnames_avg))
  tree_subgroups_list <- list()
  
  for (i in 1:n_reps) {
    
    file_name_indv <- paste0('Simulation/indv_seed', seeds[i], '.rds')
    file_name_subgroup <- paste0('Simulation/subgroup_seed', seeds[i], '.rds')
    file_name_avg <- paste0('Simulation/avg_seed', seeds[i], '.rds')
    file_name_tree_subgroup <- paste0('Simulation/tree_subgroup_seed', seeds[i], '.rds')
    
    if (file.exists(file_name_indv) & file.exists(file_name_subgroup) &
        file.exists(file_name_avg) & file.exists(file_name_tree_subgroup)) {
      indv_mat_i <- readRDS(file_name_indv)
      subgroup_mat_i <- readRDS(file_name_subgroup)
      avg_mat_i <- readRDS(file_name_avg)
      tree_subgroups_list_i <- readRDS(file_name_tree_subgroup)
    } 
    else {
      set.seed(seeds[i])
      
      # Get true direct/indirect effects
      avg_direct_true    <- mean(zeta_hat)
      avg_indirect_true  <- mean(d_hat * tau_hat)
      indv_direct_true   <- zeta_hat
      indv_indirect_true <- d_hat * tau_hat
      
      subgroup_labels <- c('white, age < 34', 'non-white, age < 34', 'age ≥ 67',
                           'white, 34 ≤ age < 67', 'non-white, 34 ≤ age < 67')
      i1 <- which(data$age < 34 & data$race == 'White')
      i2 <- which(data$age < 34 & data$race != 'White')
      i3 <- which(data$age >= 67)
      i4 <- which(data$age >= 34 & data$age < 67 & data$race == 'White')
      i5 <- which(data$age >= 34 & data$age < 67 & data$race != 'White')
      
      group1_direct_true <- mean(indv_direct_true[i1])
      group2_direct_true <- mean(indv_direct_true[i2])
      group3_direct_true <- mean(indv_direct_true[i3])
      group4_direct_true <- mean(indv_direct_true[i4])
      group5_direct_true <- mean(indv_direct_true[i5])
      subgroup_direct_true <- data.frame(group = subgroup_labels,
                             direct = c(group1_direct_true, group2_direct_true,
                                        group3_direct_true, group4_direct_true,
                                        group5_direct_true))

      group1_indirect_true <- mean(indv_indirect_true[i1])
      group2_indirect_true <- mean(indv_indirect_true[i2])
      group3_indirect_true <- mean(indv_indirect_true[i3])
      group4_indirect_true <- mean(indv_indirect_true[i4])
      group5_indirect_true <- mean(indv_indirect_true[i5])
      subgroup_indirect_true <- data.frame(group = subgroup_labels,
                       indirect = c(group1_indirect_true, group2_indirect_true,
                                    group3_indirect_true, group4_indirect_true,
                                    group5_indirect_true))
      
      # Get effects for training set (for simulating)
      mu_y_hat_train <- mu_y_hat[i_train]
      zeta_hat_train <- zeta_hat[i_train]
      d_hat_train    <- d_hat[i_train]
      mu_m_hat_train <- mu_m_hat[i_train]
      tau_hat_train  <- tau_hat[i_train]
      
      # Make training simulation set
      epsilon_y <- rnorm(nrow(data_train))
      epsilon_m <- rnorm(nrow(data_train))
      A_train <- data_train[[mediator_name]]
      
      m_sim <- mu_m_hat_train + A_train * tau_hat_train + sigma_m_hat * epsilon_m
      Y_sim <- mu_y_hat_train + A_train * zeta_hat_train + m_sim * d_hat_train + sigma_y_hat * epsilon_y
      
      data_train[[outcome_name]]  <- Y_sim
      data_train[[mediator_name]] <- m_sim
      
      # Clever covariates for training/testing set
      clever_cov <- get_clever_cov(data_train, data, model_m,
                                   mediator_name, outcome_name, treat_name)
      m0_hat <- clever_cov$m0_hat
      m1_hat <- clever_cov$m1_hat
      m0_hat_train <- m0_hat[i_train]
      m1_hat_train <- m1_hat[i_train]
      
      # Propensity score for training/testing set
      pi_hat <- get_ps(data, data, model_ps)
      pi_hat_train <- pi_hat[i_train]
      
      # Get output for simulated dataset
      out_sim <- bart_mediate_sim(data_train, data, model_m, model_y,
                                  pi_hat_train, pi_hat,
                                  m0_hat_train, m0_hat,
                                  m1_hat_train, m1_hat,
                                  'phealth', 'logY', 'smoke',
                                  n_iter, burnin)
      
      # Get simulated direct/indirect distributions
      zeta <- out_sim$zeta_samples
      zeta_train <- zeta[,i_train]
      zeta_test <- zeta[,i_test]
      
      d <- out_sim$d_samples
      d_train <- d[,i_train]
      d_test <- d[,i_test]
      
      tau <- out_sim$tau_samples
      tau_train <- tau[,i_train]
      tau_test <- tau[,i_test]
      
      delta <- d * tau
      delta_train <- delta[,i_train]
      delta_test <- delta[,i_test]
      
      avg_direct    <- rowMeans(zeta_train)
      avg_indirect  <- rowMeans(delta_train)
      indv_direct   <- zeta_test
      indv_indirect <- delta_test
      
      i1_train <- which(data_train$age < 34 & data_train$race == 'White')
      i2_train <- which(data_train$age < 34 & data_train$race != 'White')
      i3_train <- which(data_train$age >= 67)
      i4_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race == 'White')
      i5_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race != 'White')
      
      group1_direct <- rowMeans(zeta_train[,i1_train])
      group2_direct <- rowMeans(zeta_train[,i2_train])
      group3_direct <- rowMeans(zeta_train[,i3_train])
      group4_direct <- rowMeans(zeta_train[,i4_train])
      group5_direct <- rowMeans(zeta_train[,i5_train])
      
      group1_indirect <- rowMeans((delta_train)[,i1_train])
      group2_indirect <- rowMeans((delta_train)[,i2_train])
      group3_indirect <- rowMeans((delta_train)[,i3_train])
      group4_indirect <- rowMeans((delta_train)[,i4_train])
      group5_indirect <- rowMeans((delta_train)[,i5_train])
      
      # Get direct/indirect credible intervals
      avg_direct_interval    <- quantile(avg_direct, probs = c(0.025, 0.975))
      avg_indirect_interval  <- quantile(avg_indirect, probs = c(0.025, 0.975))
      
      indv_direct_interval   <- colQuantiles(indv_direct, probs = c(0.025, 0.975))
      indv_indirect_interval <- colQuantiles(indv_indirect, probs = c(0.025, 0.975))
      
      subgroup_direct <- cbind(group1_direct, group2_direct, group3_direct,
                               group4_direct, group5_direct)
      subgroup_indirect <- cbind(group1_indirect, group2_indirect, group3_indirect, 
                                 group4_indirect, group5_indirect)
      subgroup_direct_interval <- colQuantiles(subgroup_direct, probs = c(0.025, 0.975))
      subgroup_indirect_interval <- colQuantiles(subgroup_indirect, probs = c(0.025, 0.975))
      
      # Mean point estimates
      avg_direct_mean    <- mean(avg_direct)
      avg_indirect_mean  <- mean(avg_indirect)
      
      indv_direct_mean   <- colMeans(indv_direct)
      indv_indirect_mean <- colMeans(indv_indirect)
      
      subgroup_direct_mean <- colMeans(subgroup_direct)
      subgroup_indirect_mean <- colMeans(subgroup_indirect)
      
      # Catch T/F
      avg_direct_catch <- (avg_direct_true >= avg_direct_interval[1]) &
                                     (avg_direct_true <= avg_direct_interval[2])
      avg_indirect_catch <- (avg_indirect_true >= avg_indirect_interval[1]) &
                                       (avg_indirect_true <= avg_indirect_interval[2])
      
      indv_direct_catch <- unlist(lapply(1:length(indv_direct_true[i_test]), 
               function(i) (indv_direct_true[i_test][i] >= indv_direct_interval[i,1]) &
                 (indv_direct_true[i_test][i] <= indv_direct_interval[i,2])))
      indv_indirect_catch <- unlist(lapply(1:length(indv_indirect_true[i_test]), 
           function(i) (indv_indirect_true[i_test][i] >= indv_indirect_interval[i,1]) &
             (indv_indirect_true[i_test][i] <= indv_indirect_interval[i,2])))
      
      subgroup_direct_catch <- unlist(lapply(1:nrow(subgroup_direct_true), 
       function(i) (subgroup_direct_true[i,2] >= subgroup_direct_interval[i,1]) &
         (subgroup_direct_true[i,2] <= subgroup_direct_interval[i,2])))
      subgroup_indirect_catch <- unlist(lapply(1:nrow(subgroup_indirect_true), 
       function(i) (subgroup_indirect_true[i,2] >= subgroup_indirect_interval[i,1]) &
         (subgroup_indirect_true[i,2] <= subgroup_indirect_interval[i,2])))
      
      
      # Saving
      # Individual
      indv_mat_i <- cbind(rep(seeds[i], nrow(data_test)), 1:nrow(data_test),
                          indv_direct_mean, indv_direct_interval,
                          as.numeric(indv_direct_catch), indv_indirect_mean,
                          indv_indirect_interval, as.numeric(indv_indirect_catch))
      colnames(indv_mat_i) <- colnames_indv
      saveRDS(indv_mat_i, file_name_indv)
      
      # Subgroup
      subgroup_mat_i <- cbind(rep(seeds[i], length(subgroup_labels)), subgroup_labels,
                              subgroup_direct_mean, subgroup_direct_interval,
                              as.numeric(subgroup_direct_catch), subgroup_indirect_mean,
                              subgroup_indirect_interval, as.numeric(subgroup_indirect_catch))
      colnames(subgroup_mat_i) <- colnames_subgroup
      saveRDS(subgroup_mat_i, file_name_subgroup)
      
      # Average
      avg_mat_i <- c(seeds[i], avg_direct_mean, avg_direct_interval,
                     as.numeric(avg_direct_catch), avg_indirect_mean, avg_indirect_interval,
                     as.numeric(avg_indirect_catch))
      names(avg_mat_i) <- colnames_avg
      saveRDS(avg_mat_i, file_name_avg)
      
      # Tree Projection Subgroups
      data_post <- model.frame(model_m, data = data) %>% 
        select(-mediator_name) %>%
        mutate(zeta_hat = colMeans(zeta), delta_hat = colMeans(delta))
      
      tree <- rpart(delta_hat ~ . - zeta_hat, data = data_post)
      
      # True subroup values
      true_leaf_df <- data.frame(zeta_true = indv_direct_true,
                                 delta_true = indv_indirect_true,
                                 leaf = tree$where)
      
      true_leaf_means <- true_leaf_df %>% group_by(leaf) %>%
                         summarise(zeta_true = mean(zeta_true),
                                   delta_true = mean(delta_true))
      
      # Simulated subgroup values
      sim_leaf_df <- data.frame(zeta = c(t(zeta)),
                                delta = c(t(delta)),
                                iteration = rep(1:n_iter, each = nrow(data)),
                                leaf = rep(tree$where, times = n_iter))
      
      sim_leaf_means_df <- sim_leaf_df %>% group_by(iteration, leaf) %>%
                           summarize(zeta_hat = mean(zeta),
                                     delta_hat = mean(delta))
      
      sim_leaf_quantiles <- sim_leaf_means_df %>% group_by(leaf) %>%
                            summarize(zeta_mean = mean(zeta_hat),
                                      zeta_lower = quantile(zeta_hat, 0.025),
                                      zeta_upper = quantile(zeta_hat, 0.975),
                                      delta_mean = mean(delta_hat),
                                      delta_lower = quantile(delta_hat, 0.025),
                                      delta_upper = quantile(delta_hat, 0.975)) 
      
      zeta_catch <- (true_leaf_means$zeta_true >= sim_leaf_quantiles$zeta_lower) &
        (true_leaf_means$zeta_true <= sim_leaf_quantiles$zeta_upper)
      delta_catch <- (true_leaf_means$delta_true >= sim_leaf_quantiles$delta_lower) &
        (true_leaf_means$delta_true <= sim_leaf_quantiles$delta_upper)
      
      tree_subgroups_list_i <- cbind(seeds[i], sim_leaf_quantiles, zeta_catch, delta_catch)
      colnames(tree_subgroups_list_i) <- colnames_tree_subgroups
      saveRDS(tree_subgroups_list_i, file_name_tree_subgroup)
    }
    
    # Save replications into a matrix
    indv_mat[1:nrow(data_test) + nrow(data_test) * (i-1),] <- indv_mat_i
    subgroup_mat[1:5 + 5 * (i-1),] <- subgroup_mat_i
    avg_mat[i,] <- avg_mat_i
    tree_subgroups_list[[i]] <- tree_subgroups_list_i
  }
  
  tree_subgroups_mat <- do.call(rbind, tree_subgroups_list)

  return(list(
    individual = indv_mat,
    subgroups = subgroup_mat,
    average = avg_mat,
    tree_subgroups = tree_subgroups_mat
  ))
}

set.seed(1)
seeds <- sample.int(10e6, 200)
simulation <- do_simulation(meps, i_train, i_test, model_m, model_y, model_ps,
                            'phealth', 'logY', 'smoke',
                            mu_y_hat, zeta_hat, d_hat,
                            mu_m_hat, tau_hat, sigma_y_hat, sigma_m_hat, 
                            8000, 4000, 10, seeds[1:10])






