"""
Build & align a daily S&P 500 dividend-yield series (DivLog) to EXACTLY match
the ^GSPC trading-day index from 1927-12-30 to 2025-09-05 (24,536 rows).

Outputs (in the script folder):
  - sp500_dividend_daily_MATCHED_to_gspc.csv      (Date, DivLog)
  - sp500_dividend_daily_MATCHED_to_gspc.parquet

Math:
  DivLog_t = ln(TR_t/TR_{t-1}) - ln(PR_t/PR_{t-1})              # 1988+
  Monthly backfill: DivLog_m = ln(1 + (TTM_div/12)/P_month_end) # pre-1988
  We distribute DivLog_m evenly across business days (Mon–Fri) so daily logs
  add to the monthly total; then splice with TR–PR on 1988-01-05.

Notes:
  - We *force* the final index to GSPC's own trading days in [FORCE_START, FORCE_END].
  - If Yahoo's calendar differs slightly, the reindex step guarantees an exact match.
"""

from pathlib import Path
import numpy as np
import pandas as pd

# -------- You can edit these if needed --------
FORCE_START = "1927-12-30"
FORCE_END   = "2025-09-05"
EXPECTED_ROWS = 24536
PR_TICKER = "^GSPC"
TR_TICKER = "^SP500TR"  # use "^SP500NTR" for net-of-withholding
SHILLER_CSV = "https://datahub.io/core/s-and-p-500/r/data.csv"  # Shiller-derived monthly Price & TTM Dividend
# ----------------------------------------------

# yfinance import with a friendly error if missing
try:
    import yfinance as yf
except ModuleNotFoundError as e:
    raise SystemExit("Missing dependency: yfinance.\nInstall with:\n  pip install yfinance --upgrade") from e

def tidy_index(idx: pd.DatetimeIndex) -> pd.DatetimeIndex:
    """Drop timezone, normalize to midnight, sort, and dedupe."""
    try:
        if getattr(idx, "tz", None) is not None:
            try:
                idx = idx.tz_convert(None)
            except Exception:
                idx = idx.tz_localize(None)
    except Exception:
        pass
    idx = idx.normalize()
    idx = idx.sort_values()
    # ensure unique
    _, unique_pos = np.unique(idx.values, return_index=True)
    idx = idx[sorted(unique_pos)]
    return idx

def tidy_df(df: pd.DataFrame) -> pd.DataFrame:
    """Apply tidy_index to the index and drop duplicate rows by date."""
    df = df.copy()
    df.index = tidy_index(df.index)
    df = df[~df.index.duplicated(keep="last")].sort_index()
    return df

def fetch_gspc_index() -> pd.DatetimeIndex:
    """Fetch ^GSPC from Yahoo and return its trading-day index within the forced range."""
    g = yf.Ticker(PR_TICKER).history(period="max", interval="1d", auto_adjust=False)
    if g.empty:
        raise RuntimeError("No data returned for ^GSPC.")
    g = tidy_df(g)
    # clip to the forced range
    mask = (g.index >= pd.Timestamp(FORCE_START)) & (g.index <= pd.Timestamp(FORCE_END))
    g = g.loc[mask]
    idx = g.index
    if len(idx) == 0:
        raise RuntimeError("No ^GSPC dates found within the forced range.")
    return idx

def daily_divlog_from_tr_pr() -> pd.Series:
    """
    Compute daily DivLog = ln(TR_t/TR_{t-1}) - ln(PR_t/PR_{t-1})
    on the intersection of dates where both ^GSPC and ^SP500TR exist.
    """
    pr = yf.Ticker(PR_TICKER).history(period="max", interval="1d", auto_adjust=False)
    tr = yf.Ticker(TR_TICKER).history(period="max", interval="1d", auto_adjust=False)
    if pr.empty or tr.empty:
        raise RuntimeError("No data returned for PR or TR.")
    pr = tidy_df(pr)[["Adj Close"]].rename(columns={"Adj Close": "PR"})
    tr = tidy_df(tr)[["Adj Close"]].rename(columns={"Adj Close": "TR"})
    df = pr.join(tr, how="inner")
    divlog = (np.log(df["TR"]).diff() - np.log(df["PR"]).diff()).dropna()
    divlog.name = "DivLog"
    # de-dupe safety
    divlog = divlog[~divlog.index.duplicated(keep="last")]
    return divlog

def monthly_shiller_divlog() -> pd.DataFrame:
    """
    Load Shiller-derived monthly data and compute:
      DivLog_m = ln(1 + (D_TTM/12) / P_month_end)
    Ensure exactly one row per month via resample('M').
    """
    m = pd.read_csv(SHILLER_CSV, parse_dates=["Date"]).set_index("Date").sort_index()
    expected_cols = {"SP500", "Dividend"}
    if not expected_cols.issubset(m.columns):
        raise RuntimeError(f"Shiller CSV must contain columns {expected_cols}.")
    m = m.rename(columns={"SP500": "P", "Dividend": "D_TTM"}).dropna(subset=["P", "D_TTM"])
    mM = m.resample("M").agg({"P": "last", "D_TTM": "last"}).dropna()
    mM["DivLog_m"] = np.log(1.0 + (mM["D_TTM"] / 12.0) / mM["P"])  # TTM -> per-month
    return mM[["DivLog_m"]]

def expand_monthly_to_daily(mM: pd.DataFrame, until_date: pd.Timestamp) -> pd.Series:
    """
    Spread monthly log dividend yield evenly across business days (Mon–Fri) in each month,
    up to 'until_date' inclusive. Daily logs then sum to the monthly total.
    """
    mM = mM.loc[:until_date]
    pieces = []
    for month_end, row in mM.iterrows():
        period = month_end.to_period("M")
        start = period.start_time.normalize()
        end = period.end_time.normalize()
        days = pd.date_range(start, end, freq="B")  # Mon–Fri business days
        if len(days) == 0:
            continue
        per_day = row["DivLog_m"] / len(days)
        pieces.append(pd.Series(per_day, index=days))
    if not pieces:
        return pd.Series(dtype=float, name="DivLog")
    daily = pd.concat(pieces).sort_index()
    daily.index = tidy_index(daily.index)
    daily = daily[~daily.index.duplicated(keep="last")]
    daily.name = "DivLog"
    return daily

def main():
    # 0) Target index from ^GSPC within the forced range
    gspc_idx = fetch_gspc_index()

    # 1) TR–PR daily series (where available)
    div_trpr = daily_divlog_from_tr_pr()
    splice_date = div_trpr.index.min()  # first day with TR–PR (expected 1988-01-05)

    # 2) Shiller monthly -> daily backfill strictly before TR–PR starts
    mM = monthly_shiller_divlog()
    div_pre = expand_monthly_to_daily(mM, until_date=splice_date - pd.Timedelta(days=1))

    # 3) Stitch: pre-TR days + TR–PR days
    div_all = pd.concat([div_pre[div_pre.index < splice_date], div_trpr], axis=0).sort_index()
    div_all = div_all[~div_all.index.duplicated(keep="last")]
    div_all.name = "DivLog"

    # 4) Align EXACTLY to the ^GSPC trading-day index in the forced range
    div_matched = div_all.reindex(gspc_idx)
    missing = int(div_matched.isna().sum())
    if missing > 0:
        # In practice this should be zero; if not, forward-fill zero for gaps (rare holidays/cal changes).
        # Prefer 0 over ffill to avoid smearing dividends across days.
        div_matched = div_matched.fillna(0.0)

    # 5) Save
    base = Path(__file__).parent.resolve()
    out_csv = base / "sp500_dividend_daily_MATCHED_to_gspc.csv"
    out_par = base / "sp500_dividend_daily_MATCHED_to_gspc.parquet"
    div_matched.to_csv(out_csv, header=["DivLog"])
    div_matched.to_frame().to_parquet(out_par)

    # 6) Reports
    print("Saved:")
    print(f"- {out_csv}")
    print(f"- {out_par}")
    print("Forced GSPC span:", pd.Timestamp(FORCE_START).date(), "->", pd.Timestamp(FORCE_END).date())
    print("Matched rows:", len(div_matched), "(expected", EXPECTED_ROWS, ")")
    print("First date:", div_matched.index.min().date(), "| Last date:", div_matched.index.max().date())
    print("Missing rows before fill (if any):", missing)
    # sanity split
    pre_mask = div_matched.index < splice_date
    q_all  = 252 * div_matched.mean()
    q_pre  = 252 * div_matched.loc[pre_mask].mean()
    q_post = 252 * div_matched.loc[~pre_mask].mean()
    print("\nAnnualized q (continuous):")
    print(f"- Full sample : {q_all:.6f}")
    print(f"- Pre-TR      : {q_pre:.6f}")
    print(f"- TR+         : {q_post:.6f}")
    print("\nHead:")
    print(div_matched.head(5))
    print("\nTail:")
    print(div_matched.tail(5))

if __name__ == "__main__":
    main()
import pandas as pd

gspc = pd.read_csv(r"C:\Users\youss\OneDrive - University of Copenhagen\sp500_yahoo_daily_full.csv",
                   parse_dates=["Date"]).set_index("Date").sort_index()
div  = pd.read_csv(r"C:\Users\youss\OneDrive - University of Copenhagen\sp500_dividend_daily_MATCHED_to_gspc.csv",
                   parse_dates=["Date"]).set_index("Date").sort_index()

out = gspc.join(div, how="left")  # DivLog aligned to your trading days
out["DivRet"] = (out["DivLog"]).pipe(pd.Series.map, lambda x: None if pd.isna(x) else float(np.exp(x)-1))
out["DivPts"] = out["Adj Close"].shift(1) * out["DivRet"]  # dividend “index points” (where price exists)

out.to_csv(r"C:\Users\youss\OneDrive - University of Copenhagen\gspc_with_divs.csv")
print(out[["Adj Close","DivLog","DivRet","DivPts"]].tail())
