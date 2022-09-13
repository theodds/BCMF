# ## Packages ----
# 
# library(SoftBart)
# library(glmnet)
# library(truncnorm)
# library(rpart)
# library(gam)
# library(tidyverse)
# library(rpart)
# library(rpart.plot)
# library(mgcv)
# library(matrixStats)
# 
# bart_mediate_sim <- function(data_train, data_test, model_m, model_y,
#                              pi_hat_train, pi_hat_test,
#                              m0_hat_train, m0_hat_test,
#                              m1_hat_train, m1_hat_test,
#                              mediator_name, outcome_name, treat_name,
#                              n_iter, burnin) {
#   X_m_train <- quantile_normalize_bart(preprocess_df(
#     model.frame(model_m, data = data_train) %>% select(-mediator_name))[[1]])
#   X_y_train <- quantile_normalize_bart(preprocess_df(
#     cbind(model.frame(model_y, data = data_train) %>% select(-outcome_name),
#           m0_hat_train, m1_hat_train))[[1]])
#   X_mu_m_train <- cbind(X_m_train, pi_hat_train)
# 
#   X_m_test <- quantile_normalize_bart(preprocess_df(
#     model.frame(model_m, data = data_test) %>% select(-mediator_name))[[1]])
#   X_y_test <- quantile_normalize_bart(preprocess_df(
#     cbind(model.frame(model_y, data = data_test) %>% select(-outcome_name),
#           m0_hat_test, m1_hat_test))[[1]])
#   X_mu_m_test <- cbind(X_m_test, pi_hat_test)
# 
#   # Extracting mediator, outcomes, treatment and scaling
# 
#   m <- data_train[[mediator_name]]
#   m_scale <- (m - mean(m)) / sd(m)
# 
#   Y <- data_train[[outcome_name]]
#   Y_scale <- (Y - mean(Y)) / sd(Y)
# 
#   A <- data_train[[treat_name]]
# 
#   # Hypers for mu_m and tau
# 
#   hypers_tau             <- Hypers(X_m_train, m_scale, normalize_Y = FALSE)
#   hypers_mu_m            <- hypers_tau
#   opts_mu_m              <- Opts(update_s = FALSE)
#   opts_tau               <- opts_mu_m
#   opts_tau$update_sigma  <- FALSE
# 
#   # Hypers for mu_y, zeta, d
#   hypers_mu_y            <- Hypers(X_y_train, Y_scale, normalize_Y = FALSE)
#   hypers_zeta            <- hypers_mu_y
#   hypers_d               <- hypers_mu_y
#   opts_mu_y              <- Opts(update_s = FALSE)
#   opts_zeta              <- opts_d <- opts_mu_y
#   opts_zeta$update_sigma <- FALSE
#   opts_d$update_sigma    <- FALSE
# 
#   # Forests
#   forest_mu_m <- MakeForest(hypers_mu_m, opts_mu_m)
#   forest_tau  <- MakeForest(hypers_tau, opts_tau)
#   forest_mu_y <- MakeForest(hypers_mu_y, opts_mu_y)
#   forest_zeta <- MakeForest(hypers_zeta, opts_zeta)
#   forest_d    <- MakeForest(hypers_d, opts_d)
# 
#   # Initialize mu_m and tau
#   mu_m <- forest_mu_m$do_predict(X_mu_m_train)
#   tau <- forest_tau$do_predict(X_m_train)
# 
#   # Initialize mu_y, zeta, d
#   mu_y <- forest_mu_y$do_predict(X_y_train)
#   zeta <- forest_zeta$do_predict(X_y_train)
#   d <- forest_d$do_predict(X_y_train)
# 
#   # Store samples
#   mu_y_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
#   zeta_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
#   d_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
#   mu_m_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
#   tau_samples <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y_test))
#   sigma_y_samples <- rep(NA, n_iter - burnin)
#   sigma_m_samples <- rep(NA, n_iter - burnin)
# 
#   for(i in 1:n_iter) {
#     print(i)
# 
#     # Mi(a)
#     mu_m  <- forest_mu_m$do_gibbs(X_mu_m_train, m_scale - A * tau, X_mu_m_train, 1)
#     mu_m_test <- forest_mu_m$do_predict(X_mu_m_test)
#     tau <- forest_tau$do_gibbs(X_m_train[A == 1,], m_scale[A == 1] - mu_m[A == 1], X_m_train, 1)
#     tau_test <- forest_tau$do_predict(X_m_test)
# 
#     sigma_m <- forest_mu_m$get_sigma()
#     forest_tau$set_sigma(sigma_m)
# 
#     # Yi(a)
#     mu_y  <- forest_mu_y$do_gibbs(X_y_train, Y_scale - A * zeta - m_scale * d, X_y_train, 1)
#     mu_y_test <- forest_mu_y$do_predict(X_y_test)
# 
#     sigma_y <- forest_mu_y$get_sigma()
#     forest_zeta$set_sigma(sigma_y)
#     forest_d$set_sigma(sigma_y)
# 
#     zeta  <-
#       forest_zeta$do_gibbs(X_y_train[A == 1,], Y_scale[A == 1] - mu_y[A == 1] -
#                                           (m_scale * d)[A == 1], X_y_train, 1)
#     zeta_test <- forest_zeta$do_predict(X_y_test)
#     R <- (Y_scale - mu_y - A * zeta) / m_scale
#     d <- forest_d$do_gibbs_weighted(X_y_train, R, m_scale^2, X_y_train, 1)
#     d_test <- forest_d$do_predict(X_y_test)
# 
#     # Unscale
#     mu_y_unscale <- mean(Y) + (mu_y_test * sd(Y))
#     zeta_unscale <- zeta_test * sd(Y)
#     d_unscale <- d_test * sd(Y) / sd(m)
#     mu_m_unscale <- mean(m) + (mu_m_test * sd(m))
#     tau_unscale <- tau_test * sd(m)
#     sigma_y_unscale <- sigma_y * sd(Y)
#     sigma_m_unscale <- sigma_m * sd(m)
# 
#     if (i > burnin){
#       mu_y_samples[i - burnin,] <- mu_y_unscale
#       zeta_samples[i - burnin,] <- zeta_unscale
#       d_samples[i - burnin,] <- d_unscale
#       mu_m_samples[i - burnin,] <- mu_m_unscale
#       tau_samples[i - burnin,] <- tau_unscale
#       sigma_y_samples[i - burnin] <- sigma_y_unscale
#       sigma_m_samples[i - burnin] <- sigma_m_unscale
# 
#       # omega <- MCMCpack::rdirichlet(1, rep(1, nrow(X_m)))
#       # zeta_avg[i - burnin] <- sum(zeta_unscale * omega)
#       # delta_avg[i - burnin] <- sum(delta * omega)
#     }
# 
#   }
# 
#   return(list(mu_y_samples = mu_y_samples, zeta_samples = zeta_samples,
#               d_samples = d_samples, mu_m_samples = mu_m_samples,
#               tau_samples = tau_samples, sigma_y_samples = sigma_y_samples,
#               sigma_m_samples = sigma_m_samples))
# }
# 
# 
# get_clever_cov <- function(data_train, data_test, model_m,
#                            mediator_name, outcome_name, treat_name) {
#   X_m_test  <- model.frame(model_m, data = data_test) %>%
#                select(-mediator_name) %>% preprocess_df() %>% pluck("X")
#   X_m_train  <- model.frame(model_m, data = data_train) %>%
#                 select(-mediator_name) %>% preprocess_df() %>% pluck("X")
#   X_m0_train <- X_m_train[data_train[[treat_name]] == 0,]
#   X_m1_train <- X_m_train[data_train[[treat_name]] == 1,]
#   m0_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 0]
#   m1_train   <- data_train[[mediator_name]][data_train[[treat_name]] == 1]
# 
#   bart_m0 <- softbart(X_m0_train, m0_train, X_m_test)
#   bart_m1 <- softbart(X_m1_train, m1_train, X_m_test)
# 
#   m0_hat <- bart_m0$y_hat_test_mean
#   m1_hat <- bart_m1$y_hat_test_mean
# 
#   return(list(m0_hat = m0_hat, m1_hat = m1_hat))
# }
# 
# get_ps <- function(data_train, data_test, model_ps) {
#   glm_logit <- glm(model_ps, data = data_train, family = binomial)
#   pi_hat <- predict(glm_logit, newdata = data_test,type = 'response')
#   return (pi_hat)
# }
# 
# 
# # MEPS data
# meps <- readRDS("Data/meps_logy.rds")
# meps$smoke <- ifelse(meps$smoke == 2, 0, 1)
# 
# # Model for m and y
# model_m <- phealth ~ -1 + age + race_white + inc + bmi + edu + povlev
# model_y <- logY ~ -1 + age + race_white + inc + bmi + edu + povlev + phealth
# 
# # Clever covariates for entire set
# clever_cov <- get_clever_cov(meps, model_m, 'phealth', 'logY', 'smoke')
# m0_hat <- clever_cov$m0_hat
# m1_hat <- clever_cov$m1_hat
# 
# # Propensity score for entire set
# model_ps <- smoke ~ age + race_white + inc + bmi + edu + povlev
# pi_hat <- get_ps(meps, meps, model_ps)
# 
# # Get true values
# out <- bart_mediate_sim(meps, meps, model_m, model_y,
#                         pi_hat, pi_hat,
#                         m0_hat, m0_hat,
#                         m1_hat, m1_hat,
#                         'phealth', 'logY', 'smoke',
#                         8000, 4000)
# 
# mu_y_hat    <- colMeans(out$mu_y_samples)
# zeta_hat    <- colMeans(out$zeta_samples)
# d_hat       <- colMeans(out$d_samples)
# mu_m_hat    <- colMeans(out$mu_m_samples)
# tau_hat     <- colMeans(out$tau_samples)
# sigma_y_hat <- mean(out$sigma_y_samples)
# sigma_m_hat <- mean(out$sigma_m_samples)
# 
# # Ground truth for avg/indv direct and indirect effects
# # omega              <- MCMCpack::rdirichlet(nrow(out$zeta_samples), rep(1, ncol(out$zeta_samples)))
# avg_direct_true    <- mean(zeta_hat)
# avg_indirect_true  <- mean(d_hat * tau_hat)
# indv_direct_true   <- zeta_hat
# indv_indirect_true <- d_hat * tau_hat
# 
# i1 <- which(meps$age < 34 & meps$race_white == 1)
# i2 <- which(meps$age < 34 & meps$race_white == 0)
# i3 <- which(meps$age >= 67)
# i4 <- which(meps$age >= 34 & meps$age < 67 & meps$race_white == 1)
# i5 <- which(meps$age >= 34 & meps$age < 67 & meps$race_white == 0)
# 
# group1_direct_true <- mean(indv_direct_true[i1])
# group2_direct_true <- mean(indv_direct_true[i2])
# group3_direct_true <- mean(indv_direct_true[i3])
# group4_direct_true <- mean(indv_direct_true[i4])
# group5_direct_true <- mean(indv_direct_true[i5])
# 
# group1_indirect_true <- mean(indv_indirect_true[i1])
# group2_indirect_true <- mean(indv_indirect_true[i2])
# group3_indirect_true <- mean(indv_indirect_true[i3])
# group4_indirect_true <- mean(indv_indirect_true[i4])
# group5_indirect_true <- mean(indv_indirect_true[i5])
# 
# 
# # Make training & testing set
# # n <- nrow(meps)
# set.seed(1)
# i_train <- sample(1:nrow(meps), 5000)
# i_test <- c(1:nrow(meps))[-i_train]
# 
# # Effects for training set
# mu_y_hat_train <- mu_y_hat[i_train]
# zeta_hat_train <- zeta_hat[i_train]
# d_hat_train    <- d_hat[i_train]
# mu_m_hat_train <- mu_m_hat[i_train]
# tau_hat_train  <- tau_hat[i_train]


do_simulation <- function(data, i_train, i_test, model_m, model_y, model_ps,
                          mediator_name, outcome_name, treat_name,
                          mu_y_train, zeta_train, d_train, 
                          mu_m_train, tau_train, sigma_y, sigma_m,
                          n_iter, burnin, n_reps) {
  
  # Get training and testing set
  data_train <- data[i_train,]
  data_test <- data[i_test,]
  
  # Store intervals
  avg_direct_intervals    <- matrix(NA, nrow = n_reps, ncol = 2)
  avg_indirect_intervals  <- matrix(NA, nrow = n_reps, ncol = 2)
  indv_direct_intervals   <- array(NA, dim = c(nrow(data_test), 2, n_reps))
  indv_indirect_intervals <- array(NA, dim = c(nrow(data_test), 2, n_reps))
  
  group1_direct_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group2_direct_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group3_direct_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group4_direct_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group5_direct_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  
  group1_indirect_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group2_indirect_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group3_indirect_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group4_indirect_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  group5_indirect_intervals <- matrix(NA, nrow = n_reps, ncol = 2)
  
  # Store Means
  avg_direct_means    <- rep(NA, n_reps)
  avg_indirect_means  <- rep(NA, n_reps)
  indv_direct_means   <- matrix(NA, nrow = n_reps, ncol = nrow(data_test))
  indv_indirect_means <- matrix(NA, nrow = n_reps, ncol = nrow(data_test))
  
  group1_direct_means <- rep(NA, n_reps)
  group2_direct_means <- rep(NA, n_reps)
  group3_direct_means <- rep(NA, n_reps)
  group4_direct_means <- rep(NA, n_reps)
  group5_direct_means <- rep(NA, n_reps)
  
  group1_indirect_means <- rep(NA, n_reps)
  group2_indirect_means <- rep(NA, n_reps)
  group3_indirect_means <- rep(NA, n_reps)
  group4_indirect_means <- rep(NA, n_reps)
  group5_indirect_means <- rep(NA, n_reps)
  
  for (i in 1:n_reps) {
    # Make training simulation set
    epsilon_y <- rnorm(nrow(data_train))
    epsilon_m <- rnorm(nrow(data_train))
    A_train <- data_train[[mediator_name]]
    
    m_sim <- mu_m_train + A_train * tau_train + sigma_m * epsilon_m
    Y_sim <- mu_y_train + A_train * zeta_train + m_sim * d_train + sigma_y * epsilon_y
    
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
    
    avg_direct    <- rowMeans(zeta_train)
    avg_indirect  <- rowMeans(d_train * tau_train)
    indv_direct   <- zeta_test
    indv_indirect <- d_test * tau_test
    
    i1_train <- which(data_train$age < 34 & data_train$race_white == 1)
    i2_train <- which(data_train$age < 34 & data_train$race_white == 0)
    i3_train <- which(data_train$age >= 67)
    i4_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race_white == 1)
    i5_train <- which(data_train$age >= 34 & data_train$age < 67 & data_train$race_white == 0)
    
    group1_direct <- rowMeans(zeta_train[,i1_train])
    group2_direct <- rowMeans(zeta_train[,i2_train])
    group3_direct <- rowMeans(zeta_train[,i3_train])
    group4_direct <- rowMeans(zeta_train[,i4_train])
    group5_direct <- rowMeans(zeta_train[,i5_train])
    
    group1_indirect <- rowMeans((d_train * tau_train)[,i1_train])
    group2_indirect <- rowMeans((d_train * tau_train)[,i2_train])
    group3_indirect <- rowMeans((d_train * tau_train)[,i3_train])
    group4_indirect <- rowMeans((d_train * tau_train)[,i4_train])
    group5_indirect <- rowMeans((d_train * tau_train)[,i5_train])
    
    # Get direct/indirect credible intervals
    avg_direct_interval    <- quantile(avg_direct, probs = c(0.025, 0.975))
    avg_indirect_interval  <- quantile(avg_indirect, probs = c(0.025, 0.975))
    indv_direct_interval   <- colQuantiles(indv_direct, probs = c(0.025, 0.975))
    indv_indirect_interval <- colQuantiles(indv_indirect, probs = c(0.025, 0.975))
    
    group1_direct_interval <- quantile(group1_direct, probs = c(0.025, 0.975))
    group2_direct_interval <- quantile(group2_direct, probs = c(0.025, 0.975))
    group3_direct_interval <- quantile(group3_direct, probs = c(0.025, 0.975))
    group4_direct_interval <- quantile(group4_direct, probs = c(0.025, 0.975))
    group5_direct_interval <- quantile(group5_direct, probs = c(0.025, 0.975))
    
    group1_indirect_interval <- quantile(group1_indirect, probs = c(0.025, 0.975))
    group2_indirect_interval <- quantile(group2_indirect, probs = c(0.025, 0.975))
    group3_indirect_interval <- quantile(group3_indirect, probs = c(0.025, 0.975))
    group4_indirect_interval <- quantile(group4_indirect, probs = c(0.025, 0.975))
    group5_indirect_interval <- quantile(group5_indirect, probs = c(0.025, 0.975))
    
    # Mean point estimates
    avg_direct_mean    <- mean(avg_direct)
    avg_indirect_mean  <- mean(avg_indirect)
    indv_direct_mean   <- colMeans(indv_direct)
    indv_indirect_mean <- colMeans(indv_indirect)
    
    group1_direct_mean <- mean(group1_direct)
    group2_direct_mean <- mean(group2_direct)
    group3_direct_mean <- mean(group3_direct)
    group4_direct_mean <- mean(group4_direct)
    group5_direct_mean <- mean(group5_direct)
    
    group1_indirect_mean <- mean(group1_indirect)
    group2_indirect_mean <- mean(group2_indirect)
    group3_indirect_mean <- mean(group3_indirect)
    group4_indirect_mean <- mean(group4_indirect)
    group5_indirect_mean <- mean(group5_indirect)
    
    # Save each interval
    avg_direct_intervals[i,]     <- avg_direct_interval
    avg_indirect_intervals[i,]   <- avg_indirect_interval
    indv_direct_intervals[,,i]   <- indv_direct_interval
    indv_indirect_intervals[,,i] <- indv_indirect_interval
    
    group1_direct_intervals[i,] <- group1_direct_interval
    group2_direct_intervals[i,] <- group2_direct_interval
    group3_direct_intervals[i,] <- group3_direct_interval
    group4_direct_intervals[i,] <- group4_direct_interval
    group5_direct_intervals[i,] <- group5_direct_interval
    
    group1_indirect_intervals[i,] <- group1_indirect_interval
    group2_indirect_intervals[i,] <- group2_indirect_interval
    group3_indirect_intervals[i,] <- group3_indirect_interval
    group4_indirect_intervals[i,] <- group4_indirect_interval
    group5_indirect_intervals[i,] <- group5_indirect_interval
    
    # Save each mean
    avg_direct_means[i]     <- avg_direct_mean
    avg_indirect_means[i]   <- avg_indirect_mean
    indv_direct_means[i,]   <- indv_direct_mean
    indv_indirect_means[i,] <- indv_indirect_mean
    
    group1_direct_means[i] <- group1_direct_mean
    group2_direct_means[i] <- group2_direct_mean
    group3_direct_means[i] <- group3_direct_mean
    group4_direct_means[i] <- group4_direct_mean
    group5_direct_means[i] <- group5_direct_mean
    
    group1_indirect_means[i] <- group1_indirect_mean
    group2_indirect_means[i] <- group2_indirect_mean
    group3_indirect_means[i] <- group3_indirect_mean
    group4_indirect_means[i] <- group4_indirect_mean
    group5_indirect_means[i] <- group5_indirect_mean
  }
  
  return(list(
    avg_direct_intervals = avg_direct_intervals,
    avg_indirect_intervals = avg_indirect_intervals,
    indv_direct_intervals = indv_direct_intervals,
    indv_indirect_intervals = indv_indirect_intervals,
    group1_direct_intervals = group1_direct_intervals,
    group2_direct_intervals = group2_direct_intervals,
    group3_direct_intervals = group3_direct_intervals,
    group4_direct_intervals = group4_direct_intervals,
    group5_direct_intervals = group5_direct_intervals,
    group1_indirect_intervals = group1_indirect_intervals,
    group2_indirect_intervals = group2_indirect_intervals,
    group3_indirect_intervals = group3_indirect_intervals,
    group4_indirect_intervals = group4_indirect_intervals,
    group5_indirect_intervals = group5_indirect_intervals,
    avg_direct_means = avg_direct_means,
    avg_indirect_means = avg_indirect_means,
    indv_direct_means = indv_direct_means,
    indv_indirect_means = indv_indirect_means,
    group1_direct_means = group1_direct_means,
    group2_direct_means = group2_direct_means,
    group3_direct_means = group3_direct_means,
    group4_direct_means = group4_direct_means,
    group5_direct_means = group5_direct_means,
    group1_indirect_means = group1_indirect_means,
    group2_indirect_means = group2_indirect_means,
    group3_indirect_means = group3_indirect_means,
    group4_indirect_means = group4_indirect_means,
    group5_indirect_means = group5_indirect_means
  ))
}

simulation <- do_simulation(meps, i_train, i_test, model_m, model_y, model_ps,
                            'phealth', 'logY', 'smoke',
                            mu_y_hat_train, zeta_hat_train, d_hat_train,
                            mu_m_hat_train, tau_hat_train, sigma_y_hat, sigma_m_hat,
                            8000, 4000, 5)

coverage_avg_direct <- rep(NA, 5)
coverage_avg_indirect <- rep(NA, 5)
coverage_indirect_subgroup1 <- rep(NA, 5)
coverage_indirect_subgroup2 <- rep(NA, 5)
coverage_indirect_subgroup3 <- rep(NA, 5)
coverage_indirect_subgroup4 <- rep(NA, 5)
coverage_indirect_subgroup5 <- rep(NA, 5)
for (i in 1:5){
  in_avg_direct <- (avg_direct_true >= simulation$avg_direct_intervals[i,1]) &
                   (avg_direct_true <= simulation$avg_direct_intervals[i,2])
  in_avg_indirect <- (avg_indirect_true >= simulation$avg_indirect_intervals[i,1]) &
                    (avg_indirect_true <= simulation$avg_indirect_intervals[i,2])
  in_indirect_subgroup1 <- (group1_indirect_true >= simulation$group1_indirect_intervals[i,1]) &
                         (group1_indirect_true <= simulation$group1_indirect_intervals[i,2])
  in_indirect_subgroup2 <- (group2_indirect_true >= simulation$group2_indirect_intervals[i,1]) &
                         (group2_indirect_true <= simulation$group2_indirect_intervals[i,2]) 
  in_indirect_subgroup3 <- (group3_indirect_true >= simulation$group3_indirect_intervals[i,1]) &
                         (group3_indirect_true <= simulation$group3_indirect_intervals[i,2])
  in_indirect_subgroup4 <- (group4_indirect_true >= simulation$group4_indirect_intervals[i,1]) &
                           (group4_indirect_true <= simulation$group4_indirect_intervals[i,2])
  in_indirect_subgroup5 <- (group5_indirect_true >= simulation$group5_indirect_intervals[i,1]) &
                           (group5_indirect_true <= simulation$group5_indirect_intervals[i,2])
  
  coverage_avg_direct[i] = in_avg_direct
  coverage_avg_indirect[i] = in_avg_indirect
  coverage_indirect_subgroup1[i] = in_indirect_subgroup1
  coverage_indirect_subgroup2[i] = in_indirect_subgroup2
  coverage_indirect_subgroup3[i] = in_indirect_subgroup3
  coverage_indirect_subgroup4[i] = in_indirect_subgroup4
  coverage_indirect_subgroup5[i] = in_indirect_subgroup5
  
}


