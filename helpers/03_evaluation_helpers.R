perform_majority_vote_on_grid <- function(collected_grid_estimates) {
  collected_grid_estimates |> 
    group_by(lBd, uBd) |> 
    summarize(
      LRD_votes = sum(slope >= (-1)), 
      SRD_votes = sum(slope < (-1)),
      .groups = 'drop'
    ) |> 
    mutate(classification = if_else(LRD_votes >= SRD_votes, "LRD", 'SRD'))
}

# Computes metrics (Accuracy, Sensitivity, etc.)
evaluate_metrics <- function(grid_estimates, mset, infinite_variance) {
  threshold <- if (infinite_variance) (3 / 4) else (1 / 2)
  
  grouped_and_cleaned_estimates <- grid_estimates %>%
    select(H, grid_estimates) %>% 
    unnest(grid_estimates) %>% 
    mutate(
      mem = if_else(H > threshold, "LRD", "SRD"),
      mem = factor(mem, levels = c("LRD", "SRD")),
      classification = factor(classification, levels = c("LRD", "SRD"))
    ) %>% 
    group_by(lBd, uBd) 
  
  grouped_and_cleaned_estimates %>% 
    mset(truth = mem, estimate = classification) %>% 
    select(-.estimator) %>% 
    rename(metric = .metric, estimate = .estimate) %>% 
    mutate(metric = case_when(
      metric == 'accuracy' ~ 'Accuracy',
      metric == 'sens' ~ 'Sensitivity',
      metric == 'spec' ~ 'Specificity'
    ))
}


collect_metrics <- function(TMax, mset, infinite_variance, dir = 'computed_data/', estimator) {
  if (str_ends(dir, '/', negate = T)) {
    dir <- paste0(dir, '/')
  }
  
  estimator_tag <- if (estimator == 'variance') 'var' else 'GPH'
  
  locations <- list.files(
    path = dir, 
    pattern = glue::glue('{variance_tag}_grid_estimates_{estimator_tag}_{TMax}_\\d+.rds')
  ) %>% 
    paste0(dir, .)
  
  
  metrics_list <- foreach (loc = locations) %do% {
    tmp <- read_rds(loc)
    
    metrics_tmp <- evaluate_metrics(tmp, mset, infinite_variance)
  
    metrics_tmp
  } 
  
  bind_rows(metrics_list)
}

collect_and_evaluate_metrics <- function(TMax, mset, infinite_variance, estimator) {
  estimator_tag <- if (estimator == 'variance') 'var' else 'GPH'
  if (TMax <= 100) {
    grid_estimates <- read_rds(
      glue::glue('computed_data/{variance_tag}_grid_estimates_{estimator_tag}_{TMax}.rds')
    )
    metrics <- evaluate_metrics(
      grid_estimates, my_mset, infinite_variance
    )
  }
  if (TMax > 100) {
    metrics <- collect_metrics(
      TMax = TMax, my_mset, infinite_variance, estimator = estimator
    ) 
  }
  metrics
}


# Helpers for plot generation
font_family <- 'Libre Baskerville'
font_size <- 12
text_color <- 'grey40'
grid_lines_color <- 'grey80'
grid_lines_size <- unit(0.2, 'mm')
plot_metrics <- function(dat, uBd_max) {
  if (uBd_max > 200) {
    tile_border_color <- NA
  } else {
    tile_border_color <- 'grey60'
  }
  if (uBd_max <= 100) {
    legend_pos <- c(0.415, 1.38)
  } else {
    legend_pos <- c(0.42, 1.38)
  }
  dat %>% 
    ggplot(aes(lBd, uBd, fill = estimate)) +
    geom_tile(col = tile_border_color) +
    facet_wrap(vars(metric)) +
    scale_fill_viridis_c(
      limits = c(0.5, 1),
      labels = scales::label_percent()
    ) +
    coord_cartesian(expand = F, clip = 'off') +
    theme_minimal() +
    theme(
      text = element_text(
        family = font_family,
        size = font_size,
        color = text_color
      ),
      legend.position = legend_pos,
      legend.direction = 'horizontal',
      strip.text = element_text(color = text_color, margin = margin(t = 2.5, b = 0.1, unit = 'cm')),
      plot.background = element_rect(fill = 'white', color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = grid_lines_color, linetype = 2, size = grid_lines_size),
      axis.ticks.y = element_line(color = grid_lines_color, size = grid_lines_size),
      axis.ticks.length.y = unit(1, 'mm')
    ) +
    guides(
      fill = guide_colorbar(
        barwidth = unit(12, 'cm'),
        title.position = 'top'
      )
    ) +
    labs(
      x = 'Lower Cutoff',
      y = 'Upper Cutoff',
      fill = glue::glue('Metric value for time series of length {TMax}')
    )
}

