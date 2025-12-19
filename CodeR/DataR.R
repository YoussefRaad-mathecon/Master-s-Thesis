# S&P 500 (Yahoo): Load, clean, enrich, and plot (tick-label sizes matched to BS–HMM figure style)

# ---- Packages ----
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, readr, ggplot2 piping
  library(lubridate)   # dates
  library(ggplot2)     # plots
  library(scales)      # comma labels for y-axis
})

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

# ---- Load CSV ----
csv_all <- "sp500_yahoo_daily_full.csv"
stopifnot(file.exists(csv_all))
spx <- readr::read_csv(csv_all, show_col_types = FALSE)

# ---- Clean & enrich ----
# Parse Date robustly (handles "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD")
d1 <- suppressWarnings(as.Date(ymd_hms(spx$Date, quiet = TRUE)))
d2 <- suppressWarnings(ymd(spx$Date, quiet = TRUE))
spx$Date <- dplyr::coalesce(d1, d2)

# Prefer 'Adj Close'; fallback to 'Close'
if (!("Adj Close" %in% names(spx))) spx$`Adj Close` <- spx$Close
spx$AdjClose <- suppressWarnings(as.numeric(spx$`Adj Close`))

spx <- spx %>%
  arrange(Date) %>%
  distinct(Date, .keep_all = TRUE) %>%
  filter(!is.na(Date), !is.na(AdjClose))

# Continuous time axis (same construction as your other code)
spx$DateCont <- year(spx$Date) + (month(spx$Date) - 1)/12 + (day(spx$Date) - 1)/365
# Alternative:
# spx$DateCont <- year(spx$Date) + (yday(spx$Date) - 1) / 365

# Daily log-returns
spx$Return <- c(NA_real_, diff(log(spx$AdjClose)))

# ---- Style (match BS–HMM plot sizing) ----
title_size      <- 18
axis_title_size <- 13
axis_text_size  <- 15

tick_lw     <- 0.8
tick_len_pt <- 6

theme_spx <- theme_bw() +
  theme(
    plot.title        = element_text(size = title_size, hjust = 0.5),
    axis.title        = element_text(size = axis_title_size),
    axis.text         = element_text(size = axis_text_size),
    axis.ticks        = element_line(linewidth = tick_lw),
    axis.ticks.length = grid::unit(tick_len_pt, "pt"),
    panel.grid.minor  = element_blank()
  )

# ---- Timeseries plot: price levels ----
p_price <- spx %>%
  ggplot(aes(x = DateCont, y = AdjClose)) +
  geom_line(color = "#901a1E") +
  theme_spx +
  xlab("Time") +
  ylab("S&P 500 (Price Levels)") +
  ggtitle("S&P 500 vs. Time") +
  theme(
    axis.text.x = element_text(size = axis_text_size, angle = 0, vjust = 0.7),
    axis.text.y = element_text(size = axis_text_size)
  ) +
  scale_x_continuous(breaks = seq(1930, 2030, by = 10)) +
  scale_y_continuous(labels = comma)

print(p_price)

# ---- ACF plots ----
op <- par(no.readonly = TRUE)
par(cex.axis = 1.15, cex.lab = 1.25, cex.main = 1.35)

acf(na.omit(spx$AdjClose), lag.max = 200,
    main = "S&P 500 ACF (200 Lags)", ci.col = "#901A1E", ci.lwd = 2)

acf(na.omit(spx$Return), lag.max = 200,
    main = "S&P 500 ACF (200 Lags): Log-Returns", ci.col = "#901A1E", ci.lwd = 2)

par(op)

# ---- Timeseries plot: log-returns ----
p_ret <- spx %>%
  ggplot(aes(x = DateCont, y = Return)) +
  geom_line(color = "#901a1E") +
  theme_spx +
  xlab("Time") +
  ylab("S&P 500 Log-Returns") +
  ggtitle("S&P 500 Log-Returns vs. Time") +
  theme(
    axis.text.x = element_text(size = axis_text_size, angle = 0, vjust = 0.7),
    axis.text.y = element_text(size = axis_text_size)
  ) +
  scale_x_continuous(breaks = seq(1930, 2030, by = 10))

print(p_ret)

# ---- Density plots ----
p_dens_price <- ggplot(spx, aes(x = AdjClose)) +
  geom_density(fill = "#901a1E", alpha = 0.6, color = "#901a1E") +
  theme_spx +
  xlab("Closing Price") +
  ylab("Density") +
  ggtitle("S&P 500 Density: Closing Price")

print(p_dens_price)

p_dens_ret <- ggplot(spx, aes(x = Return)) +
  geom_density(fill = "#901a1E", alpha = 0.6, color = "#901a1E") +
  theme_spx +
  xlab("Daily Log-Returns") +
  ylab("Density") +
  ggtitle("S&P 500 Density: Log-Returns")

print(p_dens_ret)

# ---- Save clean data ----
save(spx, file = "sp500_clean.Rdata")

# ---- Sanity checks ----
sanity <- tibble(
  start_date = min(spx$Date, na.rm = TRUE),
  end_date   = max(spx$Date, na.rm = TRUE),
  rows       = nrow(spx)
)
print(sanity)

print(summary(spx$AdjClose))
print(summary(spx$Return))
