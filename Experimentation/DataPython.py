"""
Fetch S&P 500 (^GSPC) daily history from Yahoo Finance and save outputs
(CSV + Parquet) to the same directory as this script.
"""

from pathlib import Path
import pandas as pd
import yfinance as yf

def main():
    # --- Config ---
    sym = "^GSPC"  # S&P 500 price index

    # --- Download ---
    df = yf.Ticker(sym).history(period="max", interval="1d", auto_adjust=False)
    if df.empty:
        raise RuntimeError("No data returned. Check your internet connection or the ticker symbol.")

    # --- Tidy index (remove timezone and dupes) ---
    if getattr(df.index, "tz", None) is not None:
        try:
            df.index = df.index.tz_convert(None)
        except Exception:
            df.index = df.index.tz_localize(None)
    df = df[~df.index.duplicated(keep="last")].sort_index()

    # --- AdjClose only (if needed) ---
    px = df[["Adj Close"]].rename(columns={"Adj Close": "AdjClose"})

    # --- Output folder = script directory ---
    base = Path(__file__).parent.resolve()
    out_all = base / "sp500_yahoo_daily_full.csv"
    out_adj = base / "sp500_adjclose_daily_full.csv"
    out_par = base / "sp500_yahoo_daily_full.parquet"

    # --- Save ---
    df.to_csv(out_all)
    px.to_csv(out_adj)
    df.to_parquet(out_par)

    # --- Log ---
    print("Saved to:")
    print(f"- {out_all}")
    print(f"- {out_adj}")
    print(f"- {out_par}")
    print(f"Date range: {df.index.min().date()} -> {df.index.max().date()} | Rows: {len(df)}")
    print(df.tail(5))

if __name__ == "__main__":
    main()
