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
  w <- 2 * pi * n / length(periodogram)
  freqs <- 2 * sin(w/2)
  freqs <- freqs[freqs > 0]
  if (is_empty(n)) return(NA)
  
  x.reg <- 2 * log(freqs)
  y.reg <- log(periodogram[n]/(2 * pi))
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

compute_var_grid_estimates_helper <- function(variances, grid){
  tmp <- future_map(
    variances$variance,
    ~estimate_var_slope_grid(., grid = grid)
  )
  
  tibble(
    simu_id = variances$simu_id,
    H = variances$H,
    grid_estimates = tmp
  )
}

compute_var_grid_estimates <- function(variances, grid, TMax) {
  if (TMax <= 100) {
    # If TMax is small enough, then the grid estimates can be computed in one go
    grid_estimates <- compute_var_grid_estimates_helper(variances, grid)
    
    # Save as intermediate result
    write_rds(
      grid_estimates, 
      glue::glue('computed_data/{variance_tag}_grid_estimates_var_{TMax}.rds')
    )
  }
  
  if (TMax > 100) {
    # If TMax is not small enough, then compute grid estimates in chunks
    # Split grid into chunks
    nrow_75 <- compute_grid(75) %>% nrow()
    n_grid_slices <- ceiling(nrow(grid) / nrow_75)
    set.seed(4635635)
    grid_slices <- split(
      grid, 
      sample(1:n_grid_slices, nrow(grid), replace=T)
    )
    
    # Run calculation for each chunk and save it on disk
    for (idx in seq_along(grid_slices)) {
      grid_slice <- grid_slices[[idx]]
      grid_estimates <- compute_var_grid_estimates_helper(variances, grid_slice)
      
      # Save as intermediate result
      write_rds(
        grid_estimates, 
        glue::glue(
          'computed_data/{variance_tag}_grid_estimates_var_{TMax}_{idx}.rds'
        )
      )
      
      rm(grid_estimates)
      if (idx %% 10 == 0) {
        gc()
      }
    }
    
  }
}





safe_estimate_d_GPH <- safely(estimate_d_GPH, otherwise = NA)



estimate_GPH_slope_grid <- function(periodogram, grid)  {
  grid %>% 
    mutate(
      d_GPH = map2_dbl(lBd, uBd, ~safe_estimate_d_GPH(periodogram, .x, .y)$result),
      classification = classify_d_GPH(d_GPH)
    )
}

compute_GPH_grid_estimates <- function(periodograms, grid_GPH, TMax) {
  if (TMax <= 100) {
    # If TMax is small enough, then the grid estimates can be computed in one go
    grid_estimates_GPH <- compute_GPH_grid_estimates_helper(periodograms, grid_GPH)
    
    # Save as intermediate result
    write_rds(
      grid_estimates_GPH, 
      glue::glue('computed_data/{variance_tag}_grid_estimates_GPH_{TMax}.rds')
    )
  }
  
  if (TMax > 100) {
    # If TMax is not small enough, then compute grid estimates in chunks
    # Split grid into chunks
    nrow_75 <- compute_grid(75) %>% nrow()
    n_grid_slices <- ceiling(nrow(grid_GPH) / nrow_75)
    set.seed(4635635)
    grid_GPH_slices <- split(
      grid_GPH, 
      sample(1:n_grid_slices, nrow(grid_GPH), replace=T)
    )
    
    # Run calculation for each chunk and save it on disk
    for (idx in seq_along(grid_GPH_slices)) {
      grid_slice <- grid_GPH_slices[[idx]]
      grid_estimates_GPH <- compute_GPH_grid_estimates_helper(periodograms, grid_slice)
      
      # Save as intermediate result
      write_rds(
        grid_estimates_GPH, 
        glue::glue(
          'computed_data/{variance_tag}_grid_estimates_GPH_{TMax}_{idx}.rds'
        )
      )
      
      rm(list = c('grid_estimates_GPH'))
      if (idx %% 3 == 0) {
        gc()
      }
    }
    
  }
}


compute_GPH_grid_estimates_helper <- function(periodograms, grid_GPH){
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