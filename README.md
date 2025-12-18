# Master's Thesis: The Black-Scholes Option Pricing Model - A Markov-Switching Extension

## Overview
This Master's thesis explores extensions to the Black-Scholes option pricing model as a frequent critique of said model is it lack of adaptivity to economical turbulent environments. The extensions are via continuous and discrete state-space models often called Markov regime-switching models. The analysis is done using S\&P 500 price data and total return data.

**Practical Information**
- **Author**: Youssef Raad
- **Supervisor**: Rolf Poulsen
- **University**: University of Copenhagen
- **Department**: Mathmematical Science
- **Date**: 22<sup>nd</sup> December 2025

## Repository Structure
```
├── README.md                          
├── CodePython/                              
│   ├── DividendsPython.py
│   └── SP500DataGeneratorPython.py
├── CodeR/                              
│   ├── players.json
│   └── 
├── DataFramesExcel/                                
│   ├── sp500_clean_with_divs.csv
│   ├── sp500_clean_with_divs_1927_2019.csv
│   ├── sp500_clean_with_divs_FULL.csv
│   ├── sp500_dividend_daily_MATCHED_to_gspc.csv
│   ├── sp500_dividend_daily_backfilled.csv
│   ├── sp500_yahoo_daily_full.csv
│   ├── statewise_dividend_yields_4and5state_with_SEs.csv          
│   └── statewise_dividend_yields_with_SEs.csv
├── DataFramesR/                         
│   ├── BSSSM_loglik_comparison.RData        
│   ├── fitted_params_BSHMM_2019.RData 
│   ├── fitted_params_BSSSM_2019.RData          
│   ├── fitted_params_BSSSM_2019_beta.RData      
│   ├── fitted_params_BSSSM_2019_test.RData
│   ├── fitted_params_BSSSM_grid_2019.RData
│   ├── fitted_params_BSSSM_grid_2019_beta.RData     
│   ├── fitted_params_BS_2019.RData
│   ├── sp500_clean_with_divs.RData
│   ├── sp500_clean_with_divs_1927_2019.RData
│   └── sp500_clean_with_divs_FULL.RData
```
