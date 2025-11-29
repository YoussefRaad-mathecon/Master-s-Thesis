#------------------------------------- Packages -------------------------------------#
library(tidyverse)   # dplyr, readr, ggplot2 piping
library(lubridate)   # dates
library(ggplot2)     # plots
library(scales)      # comma labels for y-axis
# (gridExtra, latex2exp, etc. not required here unless you need them)

#---------------------------------- Working directory -------------------------------#
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

#-------------------------------------- Load CSV ------------------------------------#
csv_all <- "sp500_yahoo_daily_full.csv"
stopifnot(file.exists(csv_all))
spx <- readr::read_csv(csv_all, show_col_types = FALSE)

#------------------------------------ Clean & enrich --------------------------------#
# Parse Date robustly (handles "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD")
d1 <- suppressWarnings(as.Date(ymd_hms(spx$Date, quiet = TRUE)))
d2 <- suppressWarnings(ymd(spx$Date, quiet = TRUE))
spx$Date <- dplyr::coalesce(d1, d2)

# Prefer 'Adj Close'; fallback to 'Close'
if (!("Adj Close" %in% names(spx))) spx$`Adj Close` <- spx$Close
spx$AdjClose <- as.numeric(spx$`Adj Close`)

spx <- spx %>%
  arrange(Date) %>%
  distinct(Date, .keep_all = TRUE) %>%
  filter(!is.na(Date), !is.na(AdjClose))

# Continuous time axis (same construction as your FRED code)
spx$DateCont <- year(spx$Date) + (month(spx$Date) - 1)/12 + (day(spx$Date) - 1)/365
# Alternative:
# spx$DateCont <- year(spx$Date) + (yday(spx$Date) - 1) / 365

# Daily log-returns
spx$Return <- c(NA_real_, diff(log(spx$AdjClose)))

#------------------------------------- Timeseries plot ------------------------------#
spx %>%
  ggplot(aes(x = DateCont, y = AdjClose)) +
  geom_line(color = "#901a1E") +
  theme_bw() +
  xlab("Time") +
  ylab("S&P 500 (Adj Close)") +
  ggtitle("S&P 500 vs. Time") +
  theme(plot.title  = element_text(size = 17, hjust = 0.5)) +
  theme(axis.title  = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size = 9)) +
  scale_x_continuous(breaks = seq(1930, 2030, by = 10)) +
  scale_y_continuous(labels = comma)

#--------------------------------------- ACF plot -----------------------------------#
acf(na.omit(spx$AdjClose), lag.max = 200,
    main = "S&P 500 ACF (200 Lags)",ci.col = "#901A1E", ci.lwd = 2)

#--------------------------------------- ACF plot: returns -----------------------------------#
acf(na.omit(spx$Return), lag.max = 200,
    main = "S&P 500 ACF (200 Lags): Log-Retuns",ci.col = "#901A1E", ci.lwd = 2)

#----------------------------------- Log-returns plot -------------------------------#
spx %>%
  ggplot(aes(x = DateCont, y = Return)) +
  geom_line(color = "#901a1E") +
  theme_bw() +
  xlab("Time") +
  ylab("S&P 500 Log-Returns") +
  ggtitle("S&P 500 vs. Time") +
  theme(plot.title  = element_text(size = 17, hjust = 0.5)) +
  theme(axis.title  = element_text(size = 13)) +
  theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.7),
        axis.text.y = element_text(size = 9)) +
  scale_x_continuous(breaks = seq(1930, 2030, by = 10))

#----------------------------------- Density: returns -------------------------------#
ggplot(spx, aes(x = AdjClose)) +
  geom_density(fill = "#901a1E", alpha = 0.6, color = "#901a1E") +
  theme_bw() +
  xlab("Daily Returns") +
  ylab("Density") +
  ggtitle("S&P 500 Density") +
  theme(plot.title = element_text(size = 15, hjust = 0.5))

summary(spx$`Adj Close`)
#-------------------------------- Density: log returns ------------------------------#
ggplot(spx, aes(x = Return)) +
  geom_density(fill = "#901a1E", alpha = 0.6, color = "#901a1E") +
  theme_bw() +
  xlab("Daily Log-Returns") +
  ylab("Density") +
  ggtitle("S&P 500 Density: Log-Returns") +
  theme(plot.title = element_text(size = 15, hjust = 0.5))

#---------------------------------- Save clean data ---------------------------------#
save(spx, file = "sp500_clean.RData")

# Quick sanity check
cat("Range:", min(spx$Date), "to", max(spx$Date), "Rows:", nrow(spx), "\n")

summary(spx$`Adj Close`)
