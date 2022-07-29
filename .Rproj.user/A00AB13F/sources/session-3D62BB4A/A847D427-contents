# Computes metrics (Accuracy, Sensitivity, etc.)
evaluate_metrics <- function(grid_estimates, metric_set) {
  grid_estimates %>%
    select(H, grid_estimates) %>% 
    unnest(grid_estimates) %>% 
    mutate(
      mem = if_else(H > 1 / 2, "LRD", "SRD"),
      mem = factor(mem, levels = c("LRD", "SRD")),
      classification = factor(classification, levels = c("LRD", "SRD"))
    ) %>% 
    group_by(lBd, uBd) %>% 
    metric_set(truth = mem, estimate = classification) %>% 
    select(-.estimator) %>% 
    rename(metric = .metric, estimate = .estimate) %>% 
    mutate(metric = case_when(
      metric == 'accuracy' ~ 'Accuracy',
      metric == 'sens' ~ 'Sensitivity',
      metric == 'spec' ~ 'Specificity'
    ))
}

collect_metrics_GPH <- function(dir, TMax) {
  if (str_ends(dir, '/', negate = T)) {
    dir <- paste0(dir, '/')
  }
  list.files(dir) %>%
    str_subset(glue::glue('grid_estimates_GPH_{TMax}')) %>%
    paste0(dir, .) %>% 
    map_dfr(read_rds) %>% 
    group_by(simu_id, H) %>% 
    summarise(grid_estimates = list(bind_rows(grid_estimates)), .groups = 'drop')
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

