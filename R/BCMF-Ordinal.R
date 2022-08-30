library(SoftBart)
library(glmnet)
library(truncnorm)
library(rpart)
library(gam)
library(tidyverse)

make_data <- function(data, model_y, mediator_name, treat_name) {
  X <- model.frame(model_y, data = data)
  X[[treat_name]] <- data[[treat_name]]
  # X[['pi_hat']] <- pi_hat
  # X[['m0_hat']] <- m0_hat
  # X[['m1_hat']] <- m1_hat
  X_rep <- X %>% mutate(count = data[[mediator_name]]) %>% uncount(count)
  order <- 1:5
  m <- unlist(lapply(data[[mediator_name]], function(x) order[1:x]))
  z <- as.numeric(X_rep[[mediator_name]] == m)
  X_final <- X_rep %>% mutate(m = m) %>% select(-mediator_name)
  X_final[[mediator_name]] <- z
  X_final[['m']] <- m
  
  return(X_final)
}

bart_mediate <- function(data, model_m, model_y, pi_hat, m0_hat, m1_hat, mediator_name, outcome_name, treat_name, n_iter, burnin){
  X_m <- quantile_normalize_bart(preprocess_df(model.frame(model_m, data = data) %>% select(-mediator_name))[[1]])
  X_y <- quantile_normalize_bart(preprocess_df(
         cbind(model.frame(model_y, data = data) %>% select(-outcome_name), m0_hat, m1_hat))[[1]])
  m <- data[[mediator_name]]
  Y <- data[[outcome_name]]
  Y_scale <- (Y - mean(Y)) / sd(Y)
  A <- data[[treat_name]]
  
  # Hypers for mu_m and tau
  hypers_mu_m <- hypers_tau <- Hypers(X_m, m)
  opts_mu_m <- Opts(update_s = FALSE)
  opts_tau <- opts_mu_m
  opts_mu_m$update_sigma <- FALSE
  opts_tau$update_sigma <- FALSE
  
  # Hypers for mu_y, zeta, d
  hypers_mu_y <- hypers_zeta <- hypers_d <- Hypers(X_y, Y_scale)
  opts_mu_y <- Opts(update_s = FALSE)
  opts_zeta <- opts_d <- opts_mu_y
  opts_zeta$update_sigma <- FALSE
  opts_d$update_sigma <- FALSE
  
  # Forests
  forest_mu_m  <- MakeForest(hypers_mu_m, opts_mu_m)
  forest_tau <- MakeForest(hypers_tau, opts_tau)
  forest_mu_m$set_sigma(1)
  forest_tau$set_sigma(1)
  forest_mu_y <- MakeForest(hypers_mu_y, opts_mu_y)
  forest_zeta <- MakeForest(hypers_zeta, opts_zeta)
  forest_d <- MakeForest(hypers_d, opts_d)
  
  # Initialize mu_m and tau
  X_mu_m <- cbind(X_m, pi_hat)
  mu_m <- forest_mu_m$do_predict(X_mu_m)
  tau <- forest_tau$do_predict(X_m)
  
  # Initialize mu_y, zeta, d
  mu_y <- forest_mu_y$do_predict(X_y)
  zeta <- forest_zeta$do_predict(X_y)
  d <- forest_d$do_predict(X_y)
  
  # Matrices for individual/average direct & indirect effects
  zeta_indv <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y))
  delta_indv <- matrix(NA, nrow = n_iter - burnin, ncol = nrow(X_y))
  
  zeta_avg <- rep(NA, n_iter - burnin)
  delta_avg <- rep(NA, n_iter - burnin)
  
  for(i in 1:n_iter) {
    print(i)
    
    # Mi(a)
    Z <- rtruncnorm(length(m), a = ifelse(m == 0, -Inf, 0), b = ifelse(m == 0, 0, Inf), mean = mu_m + A * tau, sd = 1)
    mu_m  <- forest_mu_m$do_gibbs(X_mu_m, Z - A * tau, X_mu_m, 1)
    tau <- forest_tau$do_gibbs(X_m[A == 1,], Z[A == 1] - mu_m[A == 1], X_m, 1)
    
    # Yi(a)
    mu_y  <- forest_mu_y$do_gibbs(X_y, Y_scale - A * zeta - m * d, X_y, 1)
    
    sigma_y <- forest_mu_y$get_sigma()
    forest_zeta$set_sigma(sigma_y)
    forest_d$set_sigma(sigma_y)
    
    zeta  <- forest_zeta$do_gibbs(X_y[A == 1,], Y_scale[A == 1] - mu_y[A == 1] - (m * d)[A == 1], X_y, 1)
    d <- forest_d$do_gibbs(X_y[m == 1,], Y_scale[m == 1] - mu_y[m == 1] - (A * zeta)[m == 1], X_y, 1)
    
    # Direct and Indirect Effects
    zeta_unscale <- zeta * sd(Y) # direct
    d_unscale <- d * sd(Y)
    delta <- (pnorm(mu_m + tau) - pnorm(mu_m)) * d_unscale # indirect
    
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

projection_gam <- function(data, samples) {
    projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
    for (i in 1:nrow(samples)) {
      formula <- samples[i,] ~ ns(age, 5) + race_white + ns(inc, 5) + ns(bmi, 5) +
        edu + ns(povlev, 5) + phealth
      gam_fit <- gam(formula, data = data)
      projections[i,] <- predict(gam_fit, type = 'response')
    }
  
  return(projections = projections)
}

projection_tree <- function(data, model, samples) {
    projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
    for (i in 1:nrow(samples)) {
      X <- model.matrix(model, data = data)
      formula <- samples[i,] ~ .
      tree_fit <- rpart(formula, data = as.data.frame(X))
      projections[i,] <- predict(tree_fit)
    }
  
  return(projections = projections)
}

# Model for m and y
model_m <- phealth ~ -1 + age + race_white + inc + bmi + edu + povlev + m
model_y <- logY ~ -1 + age + race_white + inc + bmi + edu + povlev + phealth + m

# Estimate clever covariates with BART
data <- readRDS("Data/meps_logy.rds")
data$smoke <- ifelse(data$smoke == 2, 0, 1)
model <- logY ~ -1 + age + race_white + inc + bmi + edu + povlev + phealth
data_new <- make_data(data, model, 'phealth', 'smoke')

X_m <- quantile_normalize_bart(preprocess_df(model.frame(model_m, data = data_new) %>% select(-phealth))[[1]])
X_m0 <- X_m[data_new$smoke == 0,]
X_m1 <- X_m[data_new$smoke == 1,]
m0 <- data_new$phealth[data_new$smoke == 0]
m1 <- data_new$phealth[data_new$smoke == 1]

bart_m0 <- softbart(X = X_m0, Y = m0, X_test = X_m)
bart_m1 <- softbart(X = X_m1, Y = m1, X_test = X_m)
# bart_m0 <- readRDS("bart_m0.rds")
# bart_m1 <- readRDS("bart_m1.rds")

m0_hat <- bart_m0$y_hat_test_mean
m1_hat <- bart_m1$y_hat_test_mean

# Estimate of propensity score
glm_logit <- glm(smoke ~ age + race_white + inc + bmi + edu + povlev,
                 data = data_new, family = binomial)
pi_hat <- predict(glm_logit, type = 'response')

out_ordinal <- bart_mediate(data_new, model_m, model_y, pi_hat, m0_hat, m1_hat, 'phealth', 'logY', 'smoke', 8000, 4000)
# out_ordinal <- readRDS("Data/out_ordinal.rds")

# Traceplots
plot(out_ordinal$avg_direct)
plot(out_ordinal$avg_indirect)
plot(out_ordinal$indv_direct[,1])
plot(out_ordinal$indv_indirect[,1])

# Projections GAM
proj_gam_zeta <- projection_gam(data_new, out_ordinal$indv_direct)
proj_gam_delta <- projection_gam(data_new, out_ordinal$indv_indirect)

hist(proj_gam_zeta[,1])
hist(proj_gam_delta[,1])

# Projections Tree
proj_tree_zeta <- projection_tree(data_new, model_y, out_ordinal$indv_direct)
proj_tree_delta <- projection_tree(data_new, model_y, out_ordinal$indv_indirect)

hist(proj_tree_zeta[,1])
hist(proj_tree_delta[,1])

# Actual delta, zeta
hist(out_ordinal$avg_direct)
hist(out_ordinal$avg_indirect)

hist(out_ordinal$indv_direct[,1])
hist(out_ordinal$indv_indirect[,1])


# Waterfall Plots
get_quantiles <- function(samples, probs) {
  n <- ncol(samples)
  quantiles <- matrix(unlist(lapply(1:n, function(i) quantile(samples[,i], probs = probs))), nrow = n, byrow = TRUE)
  return (quantiles)
}

zeta_quantiles <- get_quantiles(out_ordinal$indv_direct, probs = c(0.025, 0.5, 0.975))
delta_quantiles <- get_quantiles(out_ordinal$indv_indirect, probs = c(0.025, 0.5, 0.975))

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


