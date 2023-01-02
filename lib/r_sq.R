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