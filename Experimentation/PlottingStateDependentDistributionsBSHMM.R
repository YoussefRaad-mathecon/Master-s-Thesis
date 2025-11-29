####################################################################################################################
# BS–HMM: State-dependent densities (12 models) — compact + fixed non-overlapping labels + fitted-model mixture
# (labels moved LEFT of the gutter; extended mask behind text; zero overlap; TITLES CENTERED)
####################################################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(gridExtra)
})

###--------------------------------------------------------------------------------------------------
### Paths & data
###--------------------------------------------------------------------------------------------------
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
load("fitted_params_BSHMM_2019.RData")   # list: models_BSHMM_2019
load("sp500_clean.Rdata")                # expects `spx$Return`

x  <- as.numeric(na.omit(spx$Return))
dt <- 1/252

###--------------------------------------------------------------------------------------------------
### STYLE
###--------------------------------------------------------------------------------------------------
palette_states   <- c("#DEAB7E", "#78C2AF", "#77ACC9", "#D98E96", "#B39AD8")
title_size       <- 18
axis_title_size  <- 13
axis_text_size   <- 10
annot_size       <- 4.2
line_width       <- 0.9
panel_spacing_pt <- 6
hist_alpha       <- 0.55

# Default export sizes
PNG_WIDTH_IN  <- 13.5
PNG_HEIGHT_IN <- 16
PNG_DPI       <- 200

# Gutter (right strip normally reserved for labels)
GUTTER_FRAC   <- 0.10    # 10%

###--------------------------------------------------------------------------------------------------
### Helpers
###--------------------------------------------------------------------------------------------------
unpack_bshmm <- function(fit, N, family = c("mu","sigma","mu_sigma")) {
  family <- match.arg(family); th <- fit$estimate
  if (family == "mu") {
    mu <- th[1:N]; log_sigma <- th[N+1]; eta <- th[(N+2):length(th)]
    sigma <- rep(exp(log_sigma), N)
  } else if (family == "sigma") {
    mu <- rep(th[1], N); log_sigma <- th[2:(N+1)]; eta <- th[(N+2):length(th)]
    sigma <- exp(log_sigma)
  } else {
    mu <- th[1:N]; log_sigma <- th[(N+1):(2*N)]; eta <- th[(2*N+1):length(th)]
    sigma <- exp(log_sigma)
  }
  Gamma <- diag(N); Gamma[!Gamma] <- exp(eta); Gamma <- sweep(Gamma, 1, rowSums(Gamma), "/")
  delta <- solve(t(diag(N) - Gamma + 1), rep(1, N)); delta <- as.numeric(delta / sum(delta))
  list(mu = mu, sigma = sigma, Gamma = Gamma, delta = delta)
}

title_for <- function(N, fam) {
  if (fam == "mu")        return(paste0(N, "-state BS-HMM (", "\u03BC", ")"))
  if (fam == "sigma")     return(paste0(N, "-state BS-HMM (", "\u03C3", ")"))
  if (fam == "mu_sigma")  return(paste0(N, "-state BS-HMM (", "\u03BC", " & ", "\u03C3", ")"))
  stop("Unknown family")
}

###--------------------------------------------------------------------------------------------------
### Grid & data prep
###--------------------------------------------------------------------------------------------------
xrng_all <- quantile(x, c(0.001, 0.999), na.rm = TRUE)
x_grid   <- seq(xrng_all[1], xrng_all[2], length.out = 1000)

specs <- tribble(
  ~obj,               ~N, ~Family,
  "mod2_mu",           2,  "mu",
  "mod3_mu",           3,  "mu",
  "mod4_mu",           4,  "mu",
  "mod5_mu",           5,  "mu",
  "mod2_sigma",        2,  "sigma",
  "mod3_sigma",        3,  "sigma",
  "mod4_sigma",        4,  "sigma",
  "mod5_sigma",        5,  "sigma",
  "mod2_mu_sigma",     2,  "mu_sigma",
  "mod3_mu_sigma",     3,  "mu_sigma",
  "mod4_mu_sigma",     4,  "mu_sigma",
  "mod5_mu_sigma",     5,  "mu_sigma"
) %>%
  mutate(Family = factor(Family, levels = c("mu","sigma","mu_sigma"))) %>%
  arrange(Family, N)

###--------------------------------------------------------------------------------------------------
### Panel builder (hist + state densities + MIXTURE + FIXED LABELS with extended mask)
###   — consistent spacing + header safely above labels
###--------------------------------------------------------------------------------------------------
plot_state_densities <- function(x, fit, N, family) {
  pars  <- unpack_bshmm(fit, N, family)
  means <- (pars$mu - 0.5 * pars$sigma^2) * dt
  sds   <- pars$sigma * sqrt(dt)
  
  xrng       <- xrng_all
  xrange     <- diff(xrng)
  gutter_w   <- GUTTER_FRAC * xrange
  x_data_max <- xrng[2] - gutter_w            # end of histogram/data area
  
  # Anchor near gutter
  x_label_anchor <- xrng[2] - gutter_w - 0.15 * xrange
  # Text (header + labels) moved clearly left inside the masked band
  x_label_text   <- x_label_anchor - 0.04 * xrange
  x_right_edge   <- xrng[2]
  
  # --- DENSITY GRID: FULL RANGE (no truncation) --- #
  xg_data <- x_grid
  
  dens_df <- purrr::map_dfr(seq_len(N), function(i) {
    tibble(
      x       = xg_data,
      density = dnorm(xg_data, mean = means[i], sd = sds[i]),
      i       = i
    )
  })
  
  mix_df <- tibble(
    x       = xg_data,
    density = rowSums(sapply(seq_len(N), function(i) {
      pars$delta[i] * dnorm(xg_data, mean = means[i], sd = sds[i])
    }))
  )
  
  # Histogram
  h  <- hist(x, breaks = 60, plot = FALSE, include.lowest = TRUE, right = TRUE)
  bw <- diff(h$breaks)
  hist_densities    <- h$counts / sum(h$counts) / bw
  hist_density_max  <- max(hist_densities, na.rm = TRUE)
  
  line_density_max  <- max(dens_df$density, na.rm = TRUE)
  mix_density_max   <- max(mix_df$density,  na.rm = TRUE)
  y_max_base        <- max(line_density_max, mix_density_max, hist_density_max)
  
  # EXTRA HEADROOM DEPENDING ON N (more states -> more vertical space)
  if (N <= 2) {
    y_max_factor <- 1.15
  } else if (N == 3) {
    y_max_factor <- 1.20
  } else if (N == 4) {
    y_max_factor <- 1.25
  } else {
    y_max_factor <- 1.30
  }
  y_max <- y_max_base * y_max_factor
  
  cols <- palette_states[seq_len(N)]
  
  # VERTICAL LAYOUT: S1 + equal gaps + header just above
  top_frac <- 0.88      # S1 position as fraction of y_max
  gap_frac <- 0.048     # vertical gap between labels, as fraction of y_max
  
  y_top_labels <- y_max * top_frac
  y_vals       <- y_top_labels - (seq_len(N) - 1L) * (gap_frac * y_max)
  
  labels_df <- tibble(
    x   = x_label_text,
    y   = y_vals,
    lab = paste0("S", seq_len(N), ": ", sprintf("%.1f%%", 100 * pars$delta)),
    i   = seq_len(N)
  )
  
  # Header just above S1 but not cramped
  header_y <- y_top_labels + 0.045 * y_max
  
  p <- ggplot() +
    # Histogram only over the data area (stops before label band)
    geom_histogram(
      data = data.frame(x = x[x <= x_data_max]),
      aes(x = x, y = after_stat(density)),
      bins = 60, color = "#4D5D53", fill = "#F0EAD6", alpha = hist_alpha
    ) +
    # State density lines (FULL x-range; mask will hide label band)
    geom_line(
      data = dens_df,
      aes(x = x, y = density, color = factor(i)),
      linewidth = line_width
    ) +
    # Mixture curve (FULL x-range)
    geom_line(
      data = mix_df,
      aes(x = x, y = density),
      color = "black",
      linewidth = line_width + 0.2,
      linetype = "dotted"
    ) +
    # Extended mask: start just left of the text and cover to the right edge
    geom_rect(
      aes(
        xmin = x_label_text - 0.01 * xrange,
        xmax = x_right_edge,
        ymin = 0,
        ymax = y_max
      ),
      inherit.aes = FALSE, fill = "white", color = NA
    ) +
    # Fixed labels (left-aligned)
    geom_text(
      data = labels_df,
      aes(x = x, y = y, label = lab, color = factor(i)),
      size = annot_size, hjust = 0, vjust = 0.5
    ) +
    # "Fitted model" header
    annotate(
      "text",
      x = x_label_text, y = header_y,
      label = "Fitted model",
      hjust = 0, vjust = 0.5, size = annot_size * 0.9, color = "black"
    ) +
    scale_color_manual(values = cols, guide = "none") +
    theme_bw() +
    labs(title = title_for(N, family), x = "Log-returns", y = "Density") +
    scale_x_continuous(
      breaks = pretty(xrng, n = 5),
      limits = xrng,
      expand = expansion(mult = 0)
    ) +
    scale_y_continuous(
      breaks = NULL,
      limits = c(0, y_max),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(clip = "off") +
    theme(
      plot.title       = element_text(size = title_size, hjust = 0.5),  # CENTERED TITLES
      axis.title       = element_text(size = axis_title_size),
      axis.text        = element_text(size = axis_text_size),
      panel.grid.minor = element_blank(),
      plot.margin      = margin(5, 45, 5, 5, unit = "pt")
    )
  
  p
}

###--------------------------------------------------------------------------------------------------
### Build panels
###--------------------------------------------------------------------------------------------------
plots_list <- pmap(specs, function(obj, N, Family) {
  fit <- models_BSHMM_2019[[obj]]
  stopifnot(!is.null(fit))
  plot_state_densities(x, fit, N, as.character(Family))
})

# order as rows by N (2..5), columns by family (mu, sigma, mu_sigma)
order_mat <- expand.grid(N = c(2,3,4,5), Family = c("mu","sigma","mu_sigma")) %>%
  arrange(N, match(Family, c("mu","sigma","mu_sigma")))
index_vec <- purrr::map_int(seq_len(nrow(order_mat)), function(k) {
  which(specs$N == order_mat$N[k] & specs$Family == order_mat$Family[k])
})

g <- arrangeGrob(
  grobs = plots_list[index_vec], nrow = 4, ncol = 3, top = NULL,
  heights = rep(1, 4), widths = rep(1, 3), padding = grid::unit(panel_spacing_pt, "pt")
)

###--------------------------------------------------------------------------------------------------
### EXPORT OPTION A — big ggsave (inches × dpi)
###--------------------------------------------------------------------------------------------------
ggsave(
  filename = "BSHMM_state_densities_grid_compact_labels_LEFT_OF_GUTTER.png",
  plot     = g,
  width    = 15,
  height   = 20,
  dpi      = 300,
  units    = "in",
  limitsize = FALSE
)

# Also draw to screen
grid::grid.newpage(); grid::grid.draw(g)

###--------------------------------------------------------------------------------------------------
### EXPORT OPTION B — base::png (pixels × res)
###--------------------------------------------------------------------------------------------------
png("BSHMM_state_densities_grid_compact_labels_LEFT_OF_GUTTER_big.png",
    width = 3600, height = 4800, res = 300)
grid::grid.newpage(); grid::grid.draw(g)
dev.off()

###--------------------------------------------------------------------------------------------------
### Moments table (printed to console)
###--------------------------------------------------------------------------------------------------
moments_tbl <- pmap_dfr(specs, function(obj, N, Family) {
  fit  <- models_BSHMM_2019[[obj]]
  pars <- unpack_bshmm(fit, N, as.character(Family))
  tibble(
    Family = recode(as.character(Family),
                    "mu" = "μ only", "sigma" = "σ only", "mu_sigma" = "μ & σ"),
    N = N,
    State = paste0("S", seq_len(N)),
    mu_i = pars$mu,
    sigma_i = pars$sigma,
    mean = (pars$mu - 0.5 * pars$sigma^2) * dt,
    variance = (pars$sigma^2) * dt,
    delta = pars$delta
  )
}) %>%
  arrange(match(Family, c("μ only","σ only","μ & σ")), N, State)

moments_pretty <- moments_tbl %>%
  mutate(
    mu_i = round(mu_i, 6),
    sigma_i = round(sigma_i, 6),
    mean = round(mean, 8),
    variance = round(variance, 8),
    delta = paste0(sprintf("%.1f", 100*delta), "%")
  )

print(moments_pretty, n = nrow(moments_pretty))
