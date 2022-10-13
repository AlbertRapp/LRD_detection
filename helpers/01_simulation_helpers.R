# Covariance function of fractional Gaussian oise
cov_fGN <- function(k, H, sigma = 1) {
  sigma * 
    (abs(k - 1)^(2 * H) + 
       abs(k + 1)^(2 * H)  - 
       2 * abs(k)^(2 * H)) / 2
}

# Compute covariance matrix of fGN at time points dts
compute_sigma_matrix_fGN <- function(dts, H, sigma = 1) {
  n_dts <- length(dts)
  expand_grid(
    s = dts,
    t = dts
  ) %>% 
    mutate(cov = map2_dbl(s, t, ~cov_fGN(abs(.x - .y), H = H))) %>% 
    pull(cov) %>% 
    matrix(data = ., nrow = n_dts, ncol = n_dts)
}

# Simulate fGN by simulating multivariate Gaussian distribution
# with appropriate covariance matrix
simulate_fGN <- function(dts, n = 1, H = 1 / 2, sigma = 1){
  n_dts <- length(dts)
  mu <- rep(0, n_dts)
  Sigma <- compute_sigma_matrix_fGN(dts, H, sigma)
  matrix <- MASS::mvrnorm(n, mu, Sigma)
  
  if (is.null(dim(matrix))) {
    matrix %>%
      as_tibble() %>%
      mutate(t = dts[-length(dts)], .before = 1) %>% 
      mutate(H = H)
  } else {
    matrix %>%
      as_tibble(
        .name_repair = ~as.character(0:(ncol(matrix) - 1))
      ) %>%
      pivot_longer(cols = everything(), names_to = "t") %>%
      mutate(t = as.numeric(t)) %>%
      mutate(simu_id = rep(1:nrow(matrix), each = ncol(matrix))) %>% 
      mutate(H = H)
  }
}

# Same as simulate_fGN() but for multiple Hurst parameters and
# Parallelized with {furrr} (needs multisession plan)
# Output format is a desired nested tibble
simulate_many_fGN <- function(TMax, Hs, n_seeds) {
  future_map(
    Hs, 
    ~simulate_fGN(1:TMax, n = n_seeds, .),
    .options = furrr_options(seed = 246)
  ) %>% 
    bind_rows() %>% 
    group_by(H, simu_id) %>% 
    nest() %>% 
    ungroup() %>% 
    hoist('data', 'value') %>% 
    rename(ts = value) %>% 
    select(!data)
}


# Converts a time series to a Bernoulli time series with discrete measure mu
convert_to_bernoulli <- function(ts, mu) {
  
  # Compute Y by rowwise means
  Y <-  quantile(ts, mu, names = F) %>% 
    map(., ~(ts > .x)) %>% 
    pmap_dbl(~mean(c(...))) 
  
  Y
}