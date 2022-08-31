## Packages ----

library(SoftBart)
library(glmnet)
library(truncnorm)
library(rpart)
library(gam)
library(tidyverse)
library(rpart)
library(rpart.plot)

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
  opts_mu_m$update_sigma <- FALSE
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
  forest_mu_m$set_sigma(1)
  forest_tau$set_sigma(1)
  
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
    
    mu_m  <- forest_mu_m$do_gibbs(X_mu_m, m_scale - A * tau, X_mu_m, 1)
    tau <- forest_tau$do_gibbs(X_m[A == 1,], m_scale[A == 1] - mu_m[A == 1], X_m, 1)
    
    # Yi(a)
    mu_y  <- forest_mu_y$do_gibbs(X_y, Y_scale - A * zeta - m_scale * d, X_y, 1)
    
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

## projection_gam <- function(data, samples, type) {
##     projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
##     for (i in 1:nrow(samples)) {
##       formula <- samples[i,] ~ ns(age, 5) + race_white + ns(inc, 5) + ns(bmi, 5) +
##         edu + ns(povlev, 5) + phealth
##       gam_fit <- gam(formula, data = data)
##       projections[i,] <- predict(gam_fit, type = 'response')
##     }
  
##   return(projections = projections)
## }

## projection_tree <- function(data, model_y, samples, type) {
##     projections <- matrix(NA, nrow = nrow(samples), ncol = ncol(samples))
##     for (i in 1:nrow(samples)) {
##       X <- model.matrix(model_y, data = data)
##       formula <- samples[i,] ~ .
##       tree_fit <- rpart(formula, data = as.data.frame(X))
##       projections[i,] <- predict(tree_fit)
##     }
  
##   return(projections = projections)
## }

# Model for m and y
model_m <- phealth ~ -1 + age + race_white + inc + bmi + edu + povlev
model_y <- logY ~ -1 + age + race_white + inc + bmi + edu + povlev + phealth

# Estimate clever covariates with BART
meps <- readRDS("Data/meps_logy.rds")
meps$smoke <- ifelse(meps$smoke == 2, 0, 1)

X_m <- model.frame(model_m, data = meps) %>% 
  select(-phealth) %>% preprocess_df() %>% pluck("X")

X_m0 <- X_m[meps$smoke == 0,]
X_m1 <- X_m[meps$smoke == 1,]
m0   <- meps$phealth[meps$smoke == 0]
m1   <- meps$phealth[meps$smoke == 1]

bart_m0 <- softbart(X = X_m0, Y = m0, X_test = X_m)
bart_m1 <- softbart(X = X_m1, Y = m1, X_test = X_m)
## # bart_m0 <- readRDS("bart_m0.rds")
## # bart_m1 <- readRDS("bart_m1.rds")

m0_hat <- bart_m0$y_hat_test_mean
m1_hat <- bart_m1$y_hat_test_mean

# Estimate of propensity score
glm_logit <- glm(smoke ~ age + race_white + inc + bmi + edu + povlev,
                 data = meps, family = binomial)
pi_hat <- predict(glm_logit, type = 'response')

out_var_coef <- bart_mediate(meps, model_m, model_y, pi_hat, m0_hat, 
                             m1_hat, 'phealth', 'logY', 'smoke', 8000, 4000)
# out_var_coef <- readRDS("out_var_coef.rds")

# Traceplots
plot(out_var_coef$avg_direct)
plot(out_var_coef$avg_indirect)
plot(out_var_coef$indv_direct[,1])
plot(out_var_coef$indv_indirect[,1])

## # Projections GAM
## proj_gam_zeta <- projection_gam(data, out_var_coef$indv_direct)
## proj_gam_delta <- projection_gam(data, out_var_coef$indv_indirect)

## hist(proj_gam_zeta[,1])
## hist(proj_gam_delta[,1])

## # Projections Tree
## proj_tree_zeta <- projection_tree(data, model_y, out_var_coef$indv_direct)
## proj_tree_delta <- projection_tree(data, model_y, out_var_coef$indv_indirect)

## hist(proj_tree_zeta[,1])
## hist(proj_tree_delta[,1])

## # Actual delta, zeta
## hist(out_var_coef$avg_direct)
## hist(out_var_coef$avg_indirect)

## hist(out_var_coef$indv_direct[,1])
## hist(out_var_coef$indv_indirect[,1])


# Waterfall Plots
get_quantiles <- function(samples, probs) {
  n <- ncol(samples)
  quantiles <- matrix(unlist(lapply(1:n, 
    function(i) quantile(samples[,i], probs = probs))), nrow = n, byrow = TRUE)
  return (quantiles)
}

zeta_quantiles <- get_quantiles(out_var_coef$indv_direct, probs = c(0.025, 0.5, 0.975))
delta_quantiles <- get_quantiles(out_var_coef$indv_indirect, probs = c(0.025, 0.5, 0.975))

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
  # geom_hline(yintercept = 7030, color = "blue", size = 0.3)

ggplot(df_delta_q) + 
  geom_segment(aes(x = start, y = 1:nrow(zeta_quantiles), xend = end, yend = 1:nrow(zeta_quantiles)), color = "dark gray") +
  geom_point(data = df_delta_q, aes(x = med, y = 1:nrow(zeta_quantiles)), size = 0.2, color = 'red') +
  xlab('Indirect Effect (delta)') +
  ylab('Index') +
  geom_vline(xintercept = 0, size = 1)
  # geom_hline(yintercept = 6250, color = "blue", size = 0.3)

## Posterior projection, summarizing with trees ----

meps_post <- meps %>% 
  mutate(delta_hat = colMeans(out_var_coef$indv_indirect), 
         zeta_hat = colMeans(out_var_coef$indv_direct)) %>% 
  select(age, race_white, loginc, bmi, edu, povlev, delta_hat, zeta_hat)

rpart.plot(rpart(delta_hat ~ . - zeta_hat, data = meps_post))


## Posterior projection, summarizing with a GAM ----

pdf(NULL)
project_gam <- gam(delta_hat ~ s(age, k = 10) + race_white + s(loginc, k = 10) + s(bmi, k = 10) + 
                     s(edu, k = 4) + s(povlev, k = 10), data = meps_post)

smooths <- plot(project_gam)

smooth_age <- list()
idx <- floor(seq(from = 1, to = 4000, length = 200))
for(i in 1:length(idx)) {
  print(i)
  project_gam <- gam(delta_hat ~ s(age, k = 10) + race_white + 
                       s(loginc, k = 10) + s(bmi, k = 10) + 
                       s(edu, k = 4) + s(povlev, k = 10), data = meps_post %>% mutate(delta_hat = out_var_coef$indv_indirect[i,]))
  p <- plot(project_gam)
  smooth_age[[i]] <- data.frame(x = p[[1]]$x, y = p[[1]]$fit, iteration = i)
}
dev.off()

smooth_age_df <- do.call(rbind, smooth_age)

smooth_age_summary <- smooth_age_df %>% 
  group_by(x) %>% summarize(mu = mean(y), LCL = quantile(y, 0.025), 
                            UCL = quantile(y, 0.975))

ggplot(smooth_age_df, aes(x = x, y = y)) + 
  geom_line(aes(group = factor(iteration)), alpha = .1) + 
              stat_summary(geom = "line", color = "chartreuse3", 
                           alpha = 3, size = 2)

ggplot(smooth_age_summary, aes(x = x, y = mu, ymin = LCL, ymax = UCL)) + 
  geom_line(color = "chartreuse3", size = 2) + geom_ribbon(alpha= 0.3)
