#------------------------------------- Packages -------------------------------------#
library(tidyverse)   # dplyr, readr, ggplot2, etc.
library(lubridate)   # date handling
library(scales)      # axis labels
# Parquet reader (optional, but faster/safer for types)
suppressPackageStartupMessages({
  has_arrow <- requireNamespace("arrow", quietly = TRUE)
  if (has_arrow) library(arrow)
})

#---------------------------------- Working directory -------------------------------#
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")

#---------------------------------- File paths --------------------------------------#
csv_all_gspc <- "sp500_yahoo_daily_full.csv"                     # your GSPC export
csv_div      <- "sp500_dividend_daily_MATCHED_to_gspc.csv"       # from Python
par_div      <- "sp500_dividend_daily_MATCHED_to_gspc.parquet"   # alt: Parquet

stopifnot(file.exists(csv_all_gspc))
stopifnot(file.exists(csv_div) || file.exists(par_div))

#-------------------------------------- Load GSPC -----------------------------------#
spx <- readr::read_csv(csv_all_gspc, show_col_types = FALSE)

# Robust Date parsing (handles "YYYY-MM-DD HH:MM:SS" or "YYYY-MM-DD")
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

#---------------------------------- Load Dividends ----------------------------------#
if (file.exists(par_div) && has_arrow) {
  div <- arrow::read_parquet(par_div) %>% as_tibble()
} else {
  div <- readr::read_csv(csv_div, show_col_types = FALSE)
}

# Normalize column names just in case
names(div) <- sub("^X?Date$", "Date", names(div))  # sometimes readr names first col X1
if (!("Date" %in% names(div))) {
  # assume first column is Date if unnamed
  names(div)[1] <- "Date"
}
if (!("DivLog" %in% names(div))) {
  # if the file only has one data column, rename it to DivLog
  data_cols <- setdiff(names(div), "Date")
  if (length(data_cols) == 1) names(div)[names(div) == data_cols] <- "DivLog"
}

# Parse Date robustly
d1 <- suppressWarnings(as.Date(ymd_hms(div$Date, quiet = TRUE)))
d2 <- suppressWarnings(ymd(div$Date, quiet = TRUE))
div$Date <- dplyr::coalesce(d1, d2)

div <- div %>%
  arrange(Date) %>%
  distinct(Date, .keep_all = TRUE) %>%
  filter(!is.na(Date), !is.na(DivLog))

#---------------------------------- Align & sanity checks ---------------------------#
# Force to the intersection of your requested span:
spx <- spx %>% filter(Date >= as.Date("1927-12-30"), Date <= as.Date("2025-09-05"))
div <- div %>% filter(Date >= as.Date("1927-12-30"), Date <= as.Date("2025-09-05"))

# Check lengths & alignment
stopifnot(nrow(spx) == 24536)
# Reindex DivLog to SPX dates (ensures identical length/order)
div_matched <- tibble(Date = spx$Date) %>%
  left_join(div %>% select(Date, DivLog), by = "Date")

# You should have complete coverage; if any NA remains, set to 0 (conservative)
n_missing <- sum(is.na(div_matched$DivLog))
if (n_missing > 0) {
  warning(sprintf("DivLog had %d missing rows on SPX dates; filling with 0.", n_missing))
  div_matched$DivLog[is.na(div_matched$DivLog)] <- 0
}

#--------------------------------- Combine & enrich ---------------------------------#
spx_div <- spx %>%
  select(Date, AdjClose, Open, High, Low, Volume, Close = `Adj Close`) %>%
  left_join(div_matched, by = "Date") %>%
  mutate(
    # simple dividend return from log-dividend
    DivRet = exp(DivLog) - 1,
    # dividend points (uses prior day price)
    DivPts = DivRet * dplyr::lag(AdjClose),
    # daily log-returns of price index (capital gains)
    Return = c(NA_real_, diff(log(AdjClose))),
    # continuous annualized constant-q from full sample
    q_annual = 252 * mean(DivLog, na.rm = TRUE)
  )

#------------------------------------- Reports --------------------------------------#
cat("SPX range:", min(spx_div$Date), "to", max(spx_div$Date),
    "| Rows:", nrow(spx_div), "\n")
cat("Missing DivLog after align/fill:", sum(is.na(spx_div$DivLog)), "\n")
cat(sprintf("Annualized q (continuous): %.6f\n", unique(spx_div$q_annual)))

# Optional: split estimate
q_pre  <- 252 * mean(spx_div$DivLog[spx_div$Date < as.Date("1988-01-05")], na.rm = TRUE)
q_post <- 252 * mean(spx_div$DivLog[spx_div$Date >= as.Date("1988-01-05")], na.rm = TRUE)
cat(sprintf("Annualized q pre-1988:  %.6f\n", q_pre))
cat(sprintf("Annualized q 1988+:     %.6f\n", q_post))

#------------------------------------ Quick plots -----------------------------------#
# Rolling 252d average dividend yield
spx_div <- spx_div %>%
  mutate(DivLog_roll252 = zoo::rollapply(DivLog, 252, mean, align = "right", fill = NA))

# Dividend yield (rolling) plot
spx_div %>%
  ggplot(aes(x = Date, y = DivLog_roll252 * 252)) +  # annualized
  geom_line(color = "#901a1E") +
  theme_bw() +
  labs(x = "Time", y = "Annualized dividend yield",
       title = "S&P 500 Dividend Yield (rolling 252-day, continuous)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 0.1))

#---------------------------------- Save clean data ---------------------------------#
save(spx_div, file = "sp500_clean_with_divs.RData")
readr::write_csv(spx_div, "sp500_clean_with_divs.csv")

# Final sanity check
stopifnot(nrow(spx_div) == 24536)

load("sp500_clean_with_divs.RData")  # creates object: spx_div
head(spx_div)
x
testdata <- spx_div$DivLog_roll252*252
max(na.omit(testdata))
