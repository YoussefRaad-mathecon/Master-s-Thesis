# S&P 500 Dividend Yield (constant q): point estimate + SE (i.i.d. and HAC)
# Periods:
#   (1) Full: 1927-12-30 to 2025-09-05
#   (2) Sub : 1927-12-30 to 2019-12-31

# ---- Packages ----
library(tidyverse)
library(lubridate)
library(sandwich)

suppressPackageStartupMessages({
  has_arrow <- requireNamespace("arrow", quietly = TRUE)
  if (has_arrow) library(arrow)
})

# ---- Setup ----
setwd("C:/Users/youss/OneDrive - University of Copenhagen/Master's Thesis/Experimentation")
set.seed(112233)

# ---- File paths ----
csv_all_gspc <- "sp500_yahoo_daily_full.csv"
csv_div      <- "sp500_dividend_daily_MATCHED_to_gspc.csv"
par_div      <- "sp500_dividend_daily_MATCHED_to_gspc.parquet"

stopifnot(file.exists(csv_all_gspc))
stopifnot(file.exists(csv_div) || file.exists(par_div))

# ---- Controls ----
ANNUAL_FACTOR <- 252
FILL_MISSING_DIVLOG_WITH_ZERO <- TRUE  # recommended if dividend file is daily-matched (0 on non-div days)

FULL_START <- as.Date("1927-12-30")
FULL_END   <- as.Date("2025-09-05")

SUB_START  <- as.Date("1927-12-30")
SUB_END    <- as.Date("2019-12-31")

# ---- Helpers ----
parse_date_robust <- function(x) {
  d1 <- suppressWarnings(as.Date(lubridate::ymd_hms(x, quiet = TRUE)))
  d2 <- suppressWarnings(lubridate::ymd(x, quiet = TRUE))
  dplyr::coalesce(d1, d2)
}

read_spx_yahoo <- function(path_csv) {
  spx <- readr::read_csv(path_csv, show_col_types = FALSE)
  
  spx$Date <- parse_date_robust(spx$Date)
  
  if (!("Adj Close" %in% names(spx))) spx$`Adj Close` <- spx$Close
  spx$AdjClose <- suppressWarnings(as.numeric(spx$`Adj Close`))
  
  spx %>%
    arrange(Date) %>%
    distinct(Date, .keep_all = TRUE) %>%
    filter(!is.na(Date), !is.na(AdjClose))
}

read_div_series <- function(csv_div, par_div, has_arrow) {
  if (file.exists(par_div) && has_arrow) {
    div <- arrow::read_parquet(par_div) %>% as_tibble()
  } else {
    div <- readr::read_csv(csv_div, show_col_types = FALSE)
  }
  
  if (!("Date" %in% names(div))) names(div)[1] <- "Date"
  div$Date <- parse_date_robust(div$Date)
  
  if (!("DivLog" %in% names(div))) {
    data_cols <- setdiff(names(div), "Date")
    if (length(data_cols) == 1) {
      div <- div %>% rename(DivLog = all_of(data_cols[1]))
    } else {
      stop("Could not identify 'DivLog' column in dividend file.")
    }
  }
  
  div$DivLog <- suppressWarnings(as.numeric(div$DivLog))
  
  div %>%
    arrange(Date) %>%
    distinct(Date, .keep_all = TRUE) %>%
    filter(!is.na(Date), !is.na(DivLog))
}

estimate_q_iid <- function(divlog_vec, annual_factor = 252) {
  x <- divlog_vec[!is.na(divlog_vec)]
  n <- length(x)
  
  q_hat_daily <- mean(x)
  se_daily    <- sd(x) / sqrt(n)
  
  list(
    n            = n,
    q_hat_daily  = q_hat_daily,
    se_daily     = se_daily,
    q_hat_annual = annual_factor * q_hat_daily,
    se_annual    = annual_factor * se_daily
  )
}

estimate_q_hac <- function(divlog_vec, annual_factor = 252, prewhite = FALSE, adjust = TRUE) {
  x <- divlog_vec[!is.na(divlog_vec)]
  if (length(x) <= 5) stop("Too few observations for HAC estimate.")
  
  mod <- lm(x ~ 1)
  bw  <- sandwich::bwNeweyWest(mod)
  lag <- max(0L, as.integer(round(bw)))
  
  vc <- sandwich::NeweyWest(mod, lag = lag, prewhite = prewhite, adjust = adjust)
  
  q_hat_daily <- unname(coef(mod)[1])
  se_daily    <- sqrt(diag(vc))[1]
  
  list(
    n            = length(x),
    bw           = bw,
    lag          = lag,
    q_hat_daily  = q_hat_daily,
    se_daily     = se_daily,
    q_hat_annual = annual_factor * q_hat_daily,
    se_annual    = annual_factor * se_daily
  )
}

run_period <- function(spx_all, div_all, start_date, end_date,
                       annual_factor = 252,
                       fill_missing_with_zero = TRUE,
                       label = "") {
  spx <- spx_all %>% filter(Date >= start_date, Date <= end_date)
  if (nrow(spx) == 0) stop(sprintf("[%s] No SPX rows in requested range.", label))
  
  div <- div_all %>% filter(Date >= start_date, Date <= end_date)
  
  div_matched <- tibble(Date = spx$Date) %>%
    left_join(div %>% dplyr::select(Date, DivLog), by = "Date")
  
  n_missing <- sum(is.na(div_matched$DivLog))
  if (n_missing > 0) {
    if (fill_missing_with_zero) {
      warning(sprintf("[%s] DivLog missing on %d SPX dates; filling with 0.", label, n_missing))
      div_matched$DivLog[is.na(div_matched$DivLog)] <- 0
    } else {
      warning(sprintf("[%s] DivLog missing on %d SPX dates; keeping NA.", label, n_missing))
    }
  }
  
  div_vec <- div_matched$DivLog
  
  iid <- estimate_q_iid(div_vec, annual_factor = annual_factor)
  hac <- estimate_q_hac(div_vec, annual_factor = annual_factor)
  
  tibble(
    Period     = label,
    StartDate  = as.Date(min(spx$Date)),
    EndDate    = as.Date(max(spx$Date)),
    N          = nrow(spx),
    Method     = c("IID", "HAC"),
    q_hat_ann  = c(iid$q_hat_annual, hac$q_hat_annual),
    se_ann     = c(iid$se_annual,    hac$se_annual),
    lag_used   = c(NA_integer_,      hac$lag)
  )
}

# ---- Load data once ----
spx_all <- read_spx_yahoo(csv_all_gspc)
div_all <- read_div_series(csv_div, par_div, has_arrow)

# Hard stop if the price file does not cover FULL_END
if (max(spx_all$Date, na.rm = TRUE) < FULL_END) {
  stop(sprintf("SPX file ends at %s, which is before FULL_END=%s.",
               format(max(spx_all$Date, na.rm = TRUE), "%Y-%m-%d"),
               format(FULL_END, "%Y-%m-%d")))
}

# ---- Run the two periods ----
res_full <- run_period(
  spx_all, div_all,
  start_date = FULL_START,
  end_date   = FULL_END,
  annual_factor = ANNUAL_FACTOR,
  fill_missing_with_zero = FILL_MISSING_DIVLOG_WITH_ZERO,
  label = "FULL PERIOD (1927-12-30 to 2025-09-05)"
)

res_sub <- run_period(
  spx_all, div_all,
  start_date = SUB_START,
  end_date   = SUB_END,
  annual_factor = ANNUAL_FACTOR,
  fill_missing_with_zero = FILL_MISSING_DIVLOG_WITH_ZERO,
  label = "SUB PERIOD (1927-12-30 to 2019-12-31)"
)

q_results <- bind_rows(res_full, res_sub)

# ---- Output ----
print(q_results, n = Inf)

# ---- LaTeX table output ----
# library(knitr); library(kableExtra)
# q_results %>%
#   mutate(
#     q_hat_ann = round(q_hat_ann, 8),
#     se_ann    = round(se_ann, 8)
#   ) %>%
#   select(Period, StartDate, EndDate, N, Method, q_hat_ann, se_ann, lag_used) %>%
#   knitr::kable(
#     "latex",
#     booktabs = TRUE,
#     escape   = FALSE,
#     caption  = "Constant dividend yield estimates (annualised) for two sample periods. IID uses the standard error of the sample mean; HAC uses Newey--West standard errors with data-driven bandwidth (reported as the implied lag).",
#     label    = "tab:qConstIIDHAC"
#   ) %>%
#   kableExtra::kable_styling(position = "center", latex_options = "hold_position")
