# Estimates the variances with overlapping blocks
estimate_variances <- function(Y, kMax = length(Y) / 2) {
  TMax <- length(Y)
  
  tmp <- foreach(k = 1:kMax, .combine = "rbind") %do% {
    n_blocks <- TMax - k + 1
    
    block_means <- 
      map(1:n_blocks, 
          ~mean(Y[(.):(. +  k - 1)])) %>% 
      as.numeric
    
    all_mean <- mean(block_means)
    sample_var <- 1 / n_blocks * sum((block_means - all_mean)^2)
    
    return(c(sample_var))
  }
  c(tmp)
}

# Computes the periodogram of a time series
### This function is taken from fracdiff:fdGPH()
estimate_periodogram <- function(x) {
  # error and conversion
  if (NCOL(x) > 1) 
    stop("only implemented for univariate time series")
  x <- as.numeric(na.fail(as.ts(x)))
  if (any(is.na(x))) 
    stop("NAs in x")
  
  # process
  n <- length(x)
  kk <- 1:(n-1)
  x <- x - mean(x)
  var.x <- sum(x^2)/n
  cov.x <- numeric(n - 1L)
  for (k in kk) cov.x[k] <- sum(x[1:(n - k)] * x[(1 + k):n])/n
  
  ## periodogram
  w <- 2 * pi * seq_along(x)/n
  periodogram <- numeric(n)
  for (i in seq_along(w)) { 
    periodogram[i] <- var.x + 2 * sum(cov.x * cos(w[i] * kk))
  }
  periodogram
}


compute_grid <- function(TMax) {
  tibble(
    lBd = 1:(TMax-1),
    uBd = map(lBd + 1, ~(.:TMax))
  ) %>% 
    unnest(uBd) 
}

estimate_slope <- function(variance, lBd, uBd) {
  n <- seq_along(variance)
  n <- n[between(n, lBd, uBd)]
  
  lm.fit(cbind(1, log(n)), log(variance[n])) %>%
    coefficients() %>%
    .[[2]]
}

classify_slope <- function(slope) {
  if_else(slope > (-1), "LRD", "SRD")
}


estimate_d_GPH <- function(periodogram, lBd, uBd) {
  n <- seq_along(periodogram)
  n <- n[between(n, lBd, uBd) & (periodogram > 0)]
  w <- 2 * pi * n/length(periodogram)
  
  y.reg <- log(periodogram[n]/(2 * pi))
  x.reg <- 2 * log(2 * sin(w/2))
  fit <- lm.fit(cbind(1, x.reg), y.reg)
  d.GPH <- coef(fit)[["x.reg"]]
  -d.GPH
}

classify_d_GPH <- function(d) {
  if_else(d > 0, "LRD", "SRD")
}

estimate_var_slope_grid <- function(variance, grid) {
  grid %>% 
    mutate(
      slope = map2_dbl(lBd, uBd, ~estimate_slope(variance, .x, .y)),
      classification = classify_slope(slope)
    )
}

estimate_GPH_slope_grid <- function(periodogram, grid)  {
  grid %>% 
    mutate(
      d_GPH = map2_dbl(lBd, uBd, ~estimate_d_GPH(periodogram, .x, .y)),
      classification = classify_d_GPH(d_GPH)
    )
}

compute_GPH_grid_estimates <- function(periodograms, grid_GPH){
  tmp <- future_map(
    periodograms$periodogram,
    ~estimate_GPH_slope_grid(., grid = grid_GPH)
  )
  
  tibble(
    simu_id = periodograms$simu_id,
    H = periodograms$H,
    grid_estimates = tmp
  )
}