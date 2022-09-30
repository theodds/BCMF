## 249
softbart_bcmf <- function(formula_y, formula_m, trt, data, test_data, 
                          hypers_y = NULL, hypers_m = NULL, opts = NULL, 
                          verbose = TRUE) {
  
  dv_y <- dummyVars(formula_y, data)
  dv_m <- dummyVars(formula_m, data)
  terms_y <- attr(dv_y$terms, "term.labels")
  terms_m <- attr(dv_m$terms, "term.labels")
  group_y <- dummy_assign(dv_y)
  group_m <- dummy_assign(dv_m)
  
  treated <- which(trt == 1)
  stopifnot(all(trt == 1 | trt == 0))
  
    suppressWarnings({
      Xy_train <- predict(dv_y, data)
      Xy_test  <- predict(dv_y, test_data)
      Xm_train <- predict(dv_m, data)
      Xm_test  <- predict(dv_m, test_data)
    })
    
    Y_train <- model.response(model.frame(formula_y, data))
    Y_test  <- model.response(model.frame(formula_y, test_data))
    M_train <- model.response(model.frame(formula_m, data))
    M_test  <- model.response(model.frame(formula_m, test_data))
    
    stopifnot(is.numeric(Y_train))
    stopifnot(is.numeric(M_train))
    
    mean_Y <- mean(Y_train)
    mean_M <- mean(M_train)
    sd_Y <- sd(Y_train)
    sd_M <- sd(M_train)
    Y_train <- (Y_train - mean_Y) / sd_Y
    Y_test  <- (Y_test - mean_Y) / sd_Y
    M_train <- (M_train - mean_M) / sd_M
    M_test <- (M_test - mean_M) / sd_M
    
    if(is.null(hypers_y)) {
      hypers_y <- Hypers(X = cbind(Xy_train, M_train, trt), 
                         Y = Y_train, normalize_Y = FALSE)
    }
    hypers_y$sigma_mu <- 3 / hypers_y$k / sqrt(hypers_y$num_tree)
    hypers_y$group <- group_y
    
    if(is.null(hypers_m)) {
      hypers_m <- Hypers(X = cbind(Xm_train, trt), 
                         Y = M_train, normalize_Y = FALSE)
    }
    hypers_m$sigma_mu <- 3 / hypers_m$k / sqrt(hypers_m$num_tree)
    hypers_m$group <- group_m
    
    if(is.null(opts)) {
      opts <- Opts()
      opts$update_s <- FALSE
    }
    opts$num_print <- .Machine$integer.max
    
    ## Normalize!
    
    make_01_norm <- function(x) {
      a <- min(x)
      b <- max(x)
      return(function(y) (y - a) / (b - a))
    }
    
    ecdfs_y   <- list()
    for(i in 1:ncol(Xy_train)) {
      ecdfs_y[[i]] <- ecdf(Xy_train[,i])
      if(length(unique(Xy_train[,i])) == 1) ecdfs_y[[i]] <- identity
      if(length(unique(Xy_train[,i])) == 2) ecdfs_y[[i]] <- make_01_norm(Xy_train[,i])
    }
    for(i in 1:ncol(Xy_train)) {
      Xy_train[,i] <- ecdfs_y[[i]](Xy_train[,i])
      Xy_test[,i] <- ecdfs_y[[i]](Xy_test[,i])
    }

    ecdfs_m   <- list()
    for(i in 1:ncol(Xm_train)) {
      ecdfs_m[[i]] <- ecdf(Xm_train[,i])
      if(length(unique(Xm_train[,i])) == 1) ecdfs_m[[i]] <- identity
      if(length(unique(Xm_train[,i])) == 2) ecdfs_m[[i]] <- make_01_norm(Xm_train[,i])
    }
    for(i in 1:ncol(Xy_train)) {
      Xm_train[,i] <- ecdfs_m[[i]](Xm_train[,i])
      Xm_test[,i] <- ecdfs_m[[i]](Xm_test[,i])
    }
    
    ## Make forests ----
    mu_y_forest  <- MakeForest(hypers_y, opts, FALSE)
    zeta_forest  <- MakeForest(hypers_y, opts, FALSE)
    d_forest     <- MakeForest(hypers_y, opts, FALSE)
    mu_m_forest  <- MakeForest(hypers_m, opts, FALSE)
    tau_m_forest <- MakeForest(hypers_m, opts, FALSE)
    
    ## Output ----
    
    mu_y_train  <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    mu_y_test   <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    zeta_train  <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    zeta_test   <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    d_train     <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    d_test      <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    delta_train <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    delta_test  <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    mu_m_train  <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    mu_m_test   <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    tau_m_train <- matrix(NA, nrow = opts$num_save, ncol = length(Y_train))
    tau_m_test  <- matrix(NA, nrow = opts$num_save, ncol = length(Y_test))
    sigma_y_out <- numeric(opts$num_save)
    sigma_m_out <- numeric(opts$num_save)
    
    ## Preparing to run the chain
    mu_y  <- as.numeric(mu_y_forest$do_predict(Xy_train))
    zeta  <- as.numeric(zeta_forest$do_predict(Xy_train))
    d     <- as.numeric(d_forest$do_predict(Xy_train))
    mu_m  <- as.numeric(mu_m_forest$do_predict(Xm_train))
    tau_m <- as.numeric(tau_m_forest$do_predict(Xm_train))
    
    sigma_y <- mu_y_forest$get_sigma()
    sigma_m <- mu_m_forest$get_sigma()
    
    ## Warmup ----
    
    pb <- progress_bar$new(
      format = "  warming up [:bar] :percent eta: :eta",
      total = opts$num_burn, clear = FALSE, width= 60)
    
    for(i in 1:opts$num_burn) {
      if(verbose) pb$tick()
      
      ## Update mu_y
      R       <- Y_train - zeta * trt - d * M_train
      mu_y    <- mu_y_forest$do_gibbs(Xy_train, R, Xy_train, 1)
      sigma_y <- mu_y_forest$get_sigma()
      
      ## Update zeta
      zeta_forest$set_sigma(sigma_y)
      
      R    <- Y_train - mu_y - d * M_train
      zeta <- zeta_forest$do_gibbs(Xy_train[treated, ], R[treated], Xy_train, 1)
      
      ## Update d
      d_forest$set_sigma(sigma_y)
      
      R <- (Y_train - mu_y - trt * zeta) / M_train
      d <- d_forest$do_gibbs_weighted(Xy_train, R, M_train^2, Xy_train, 1)
      
      ## Update mu_m
      R       <- M_train - tau_m * trt
      mu_m    <- mu_m_forest$do_gibbs(Xm_train, R, Xm_train, 1)
      sigma_m <- mu_m_forest$get_sigma()
      
      ## Update tau_m
      tau_m_forest$set_sigma(sigma_m)
      
      R <- M_train - mu_m
      tau_m <- tau_m_forest$do_gibbs(Xm_train[treated,], R[treated], Xm_train, 1)
    }
  
    ## Save ----
      
    pb <- progress_bar$new(
      format = "  saving [:bar] :percent eta: :eta",
      total = opts$num_save, clear = FALSE, width= 60)
    
    for(i in 1:opts$num_save) {
      for(j in 1:opts$num_thin) {
        if(verbose) pb$tick()
        
        ## Update mu_y
        R       <- Y_train - zeta * trt - d * M_train
        mu_y    <- mu_y_forest$do_gibbs(Xy_train, R, Xy_train, 1)
        sigma_y <- mu_y_forest$get_sigma()
        
        ## Update zeta
        zeta_forest$set_sigma(sigma_y)
        
        R    <- Y_train - mu_y - d * M_train
        zeta <- zeta_forest$do_gibbs(Xy_train[treated, ], R[treated], Xy_train, 1)
        
        ## Update d
        d_forest$set_sigma(sigma_y)
        
        R <- (Y_train - mu_y - trt * zeta) / M_train
        d <- d_forest$do_gibbs_weighted(Xy_train, R, M_train^2, Xy_train, 1)
        
        ## Update mu_m
        R       <- M_train - tau_m * trt
        mu_m    <- mu_m_forest$do_gibbs(Xm_train, R, Xm_train, 1)
        sigma_m <- mu_m_forest$get_sigma()
        
        ## Update tau_m
        tau_m_forest$set_sigma(sigma_m)
        
        R <- M_train - mu_m
        tau_m <- tau_m_forest$do_gibbs(Xm_train[treated,], R[treated], Xm_train, 1)        
      }
      
      
      ## NOT QUITE RIGHT (ALSO: ADD AS.NUMERIC)
      # mu_y_train[i,] <- mu_y * sd_Y + mean_Y - mu_M * sd_Y / sd_M * d
      # mu_y_test[i,]  <- mu_y_forest$do_predict(Xy_train) * sd_Y + mean_Y - mu_M * sd_Y / sd_M * d_forest$do_predict(Xy_train)
      zeta_train[i,] <- zeta * sd_Y
      zeta_test[i,]  <- zeta_forest$do_predict(Xy_train) * sd_Y
      d_train[i,]    <- d * sd_Y / sd_M
      d_test[i,]     <- d_forest$do_predict(Xy_train) * sd_Y / sd_M
      tau_m_train[i,] <- tau_m * sd_M
      tau_m_test[i,] <- tau_m_forest$do_predict(Xm_train) * sd_M
      delta_train[i,] <- d_train[i,] * tau_m_train[i,]
      delta_test[i,] <- d_test[i,] * tau_m_test[i,]
      mu_m_train[i,] <- sd_M * mu_m + mean_M
      mu_m_test[i,] <- mu_m_forest$do_predict(Xm_train) * sd_M + mean_M
      sigma_y_out[i] <- sigma_y * sd_Y
      sigma_m_out[i] <- sigma_m * sd_M
      mu_y_train[i,] <- mu_y * sd_Y + mean_Y - mean_M * d_train[i,]
      mu_y_test[i,] <- mu_y_forest$do_predict(Xy_train) * sd_Y + mean_Y - mean_M * d_test[i,]
    }
    
    out <- list(mu_y_train = mu_y_train, mu_y_test = mu_y_test, 
                zeta_train = zeta_train, zeta_test = zeta_test, 
                d_train = d_train, d_test = d_test, tau_m_train = tau_m_train, 
                tau_m_test = tau_m_test, delta_train = delta_train, 
                delta_test = delta_test, mu_m_train = mu_m_train, 
                mu_m_test = mu_m_test, sigma_y = sigma_y_out, 
                sigma_m = sigma_m_out)
    
    class(out) <- "vc_bcmf"
    
    return(out)
    
}

dummy_assign <- function(dummy) {
  terms <- attr(dummy$terms, "term.labels")
  group <- list()
  j     <- 0
  for(k in terms) {
    if(k %in% dummy$facVars) {
      group[[k]] <- rep(j, length(dummy$lvls[[k]]))
    } else {
      group[[k]] <- j
    }
    j <- j + 1
  }
  return(do.call(c, group))
}
